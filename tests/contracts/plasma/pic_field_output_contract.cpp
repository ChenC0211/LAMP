#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/IO/PicFieldOutput.h"
#include "Modules/Plasma/IO/PicFieldSelection.h"
#include "Modules/Plasma/IO/PicFieldFrame.h"
#include "Modules/Plasma/IO/PicProfiles.h"
#include "Modules/Plasma/IO/PicFieldOutputQueue.h"
#include "Modules/Plasma/IO/PicOutputReconstruction.h"
#include "Modules/Plasma/IO/PicSpeciesSeries.h"
#include "Modules/Plasma/IO/PicStepSeries.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicOverrides.h"
#include "Modules/Plasma/IO/PicOutputManifest.h"
#include "Modules/Plasma/IO/PicPhysicsDiagnostics.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/Diagnostics/PicParticleFields.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"
#include "Core/IO/ConsoleControl.h"
#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include "Modules/Plasma/Boundary/PicScheduledGeometry.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Solvers/Poisson/PoissonFVSolverCellCentered.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicFieldSolve.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Poisson/PoissonBC.h"
#include "Core/Field/VertField.h"
#include "Core/Field/CellField.h"
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <type_traits>

namespace {
void check(bool value,const char* why) { if(!value) throw std::runtime_error(why); }
void checkLiveParticleSamples() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<ParticleType> species;
    species.emplace_back("first",2.*PhysicalConstants::eV_J,-1.);
    species.emplace_back("second",4.*PhysicalConstants::eV_J,2.);
    species.emplace_back("missing",2.*PhysicalConstants::eV_J,0.);
    Particle p;p.velocity={1.,2.};p.vz=2.;p.w=0.;p.rngState=7654;
    species[0].particle.push_back(p);p.w=3.;species[0].particle.push_back(p);
    p.w=nan;species[1].particle.push_back(p);p.w=-1.;species[2].particle.push_back(p);
    species[0].number=999; // Container size, not a possibly stale cached count.
    static_assert(std::is_same_v<decltype(countPicMacroParticles(species)),std::uint64_t>);
    check(countPicMacroParticles(species)==4 && countPicMacroParticles({})==0,"integer live macro counts");
    check(PicSpeciesCount::sample(species[0],2.).macro==2 && PicSpeciesCount::sample(species[0],2.).physical==5.,
        "live counts must use positive particle weights and original default fallback");
    check(PicSpeciesCount::sample(species[1],-2.).physical==-2.
        && std::isnan(PicSpeciesCount::sample(species[1],nan).physical),"count fallback must not sanitize legacy values");
    PicDiagnosticWindow window(species.size());window.beginFromParticles(.125,species,{2.});
    check(window.totals().speciesMacroBegin==std::vector<long long>({2,1,1})
        && window.totals().speciesPhysBegin==std::vector<double>({5.,1.,1.}),"window samples positional default weights");
    const std::unordered_map<std::string,std::size_t> index{{"first",1},{"second",0},{"missing",99}};
    std::vector<std::string> order;std::vector<double> weights,energies,charges;
    observeExplicitPicBulkEnergy(species,index,{2.,7.},2.,
        [&](const std::string& name,double charge,double energy,double weight) {
            order.push_back(name);weights.push_back(weight);energies.push_back(energy);charges.push_back(charge);
        });
    check(order==std::vector<std::string>({"first","first","second","missing"})
        && weights==std::vector<double>({7.,3.,2.,1.}) && charges==std::vector<double>({-1.,-1.,2.,0.}),
        "explicit bulk samples must preserve name-index defaults and particle order");
    check(std::abs(energies[0]-36.)<1e-12 && std::abs(energies[2]-72.)<1e-12,
        "bulk samples must include three velocity components and the original scale");
    int calls=0;bool rejected=false;
    try { observeExplicitPicBulkEnergy(species,index,{2.,7.},2.,
        [&](const std::string&,double,double,double) { if(++calls==2)throw std::runtime_error("observer failure"); }); }
    catch(const std::runtime_error&) { rejected=true; }
    check(rejected && calls==2 && species[0].particle.size()==2 && species[0].particle[0].rngState==7654
        && species[0].number==999 && species[0].particle[0].velocity.x==1.,"sampling failure must propagate without mutating particles");
    observeExplicitPicBulkEnergy(species,index,{},2.,{});
    weights.clear();observeExplicitPicBulkEnergy(species,{}, {},2.,
        [&](const std::string&,double,double,double weight) { weights.push_back(weight); });
    check(weights==std::vector<double>({1.,3.,1.,1.}),"missing names must keep unit fallback weights");
    species[0].particle.clear();window.beginFromParticles(3.,species,{9.});
    check(window.startTimeCode()==.125 && window.totals().speciesPhysBegin[0]==5.,"active interval must retain its first snapshot");
    window.finishReport();window.beginFromParticles(3.,species,{9.});
    check(window.startTimeCode()==3. && window.totals().speciesMacroBegin[0]==0
        && window.totals().speciesPhysBegin[0]==0.,"new interval must resample live particles");
    p.w=1e16;species[0].particle={p};p.w=1.;species[0].particle.push_back(p);species[0].particle.push_back(p);
    check(PicSpeciesCount::sample(species[0],2.).physical==1e16,"physical counts must retain sequential floating-point addition");
}
void checkImplicitParticleSummary() {
    const double eV=PhysicalConstants::eV_J;
    std::vector<ParticleType> species;
    species.emplace_back("electron",2.*eV,-2.);species.emplace_back("ion",4.*eV,4.);
    species.emplace_back("neutral",eV,0.);species.emplace_back("tiny",eV,5e-31);
    Particle p;p.velocity={1.,2.};p.vz=2.;p.w=2.;p.rngState=42;p.trackId=51;
    species[0].particle.push_back(p);p.w=0.;species[0].particle.push_back(p);
    p.w=-1.;species[1].particle.push_back(p);p.w=2.;species[2].particle.push_back(p);
    p.w=1.;species[3].particle.push_back(p);
    std::vector<std::string> names;std::vector<double> charges,energies,weights;
    const auto result=summarizeImplicitPicParticles(species,{3.},2.,
        [&](const std::string& name,double q,double energy,double w) {
            names.push_back(name);charges.push_back(q);energies.push_back(energy);weights.push_back(w);
        });
    const auto near=[](double a,double b){return std::abs(a-b)<=1e-13*std::max(std::abs(a),std::abs(b));};
    check(result.macroCharge==-6. && near(result.kinetic.electronJ,180.*eV)
        && near(result.kinetic.ionJ,72.*eV) && near(result.kinetic.neutralJ,54.*eV),
        "implicit charge and kinetic sums retain missing/default weights, 3D speed and near-neutral classification");
    check(names==std::vector<std::string>{"electron","electron","ion","neutral","tiny"}
        && charges==std::vector<double>{-2.,-2.,4.,0.,5e-31} && weights==std::vector<double>{2.,3.,1.,2.,1.},
        "bulk observer retains species/particle order, charge and fallback weight");
    check(energies.size()==5 && near(energies[0],36.) && near(energies[1],36.)
        && near(energies[2],72.) && near(energies[3],18.) && near(energies[4],18.),
        "bulk observer receives unweighted per-particle eV samples");
    const auto silent=summarizeImplicitPicParticles(species,{3.},2.);
    check(silent.macroCharge==result.macroCharge && silent.kinetic.electronJ==result.kinetic.electronJ
        && silent.kinetic.ionJ==result.kinetic.ionJ && silent.kinetic.neutralJ==result.kinetic.neutralJ,
        "disabled observer does not change numerical accumulation");
    int seen=0;bool failed=false;
    try {summarizeImplicitPicParticles(species,{3.},2.,[&](const std::string&,double,double,double) {
        if(++seen==2)throw std::runtime_error("observer failed");});}catch(const std::runtime_error&){failed=true;}
    check(failed && seen==2 && species[0].particle[0].rngState==42 && species[0].particle[0].trackId==51
        && species[0].particle[0].velocity.x==1.,"observer failure stops at its original particle without physical mutation");
    const auto negative=summarizeImplicitPicParticles(species,{-3.},2.);
    check(negative.macroCharge==6. && near(negative.kinetic.electronJ,72.*eV),
        "negative fallback contributes charge but is rejected by the existing kinetic filter");
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite=summarizeImplicitPicParticles(species,{nan},2.);
    check(std::isnan(nonfinite.macroCharge) && near(nonfinite.kinetic.electronJ,72.*eV),
        "nonfinite fallback charge propagates while kinetic filter skips its energy");
    species[0].particle[0].velocity.x=nan;const auto invalidEnergy=summarizeImplicitPicParticles(species,{3.},2.);
    check(invalidEnergy.macroCharge==-6. && near(invalidEnergy.kinetic.electronJ,108.*eV),
        "nonfinite kinetic energy does not suppress charge accounting");
    const auto empty=summarizeImplicitPicParticles({}, {},2.);
    check(empty.macroCharge==0. && empty.kinetic.total()==0.,"empty particle sample starts with zero results");
}
void checkDiagnosticSample(bool implicit) {
    using Diagnostics=ExplicitPIC::Diagnostics::PicDiagnostics;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=64;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    CellField<double> rhs(mesh,0.),phi(mesh,0.),charge(mesh,0.);
    phi.value[0]=-2.;phi.value[1]=5.;phi.value[2]=1.;phi.value[3]=3.;
    UnitsScales scales;scales.L0_m=2.;scales.T0_s=1.;
    auto units=std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));
    NormalizationConfig config;config.slabThickness_m=3.;
    const auto norm=NormalizationContext::MakeFromUnits(units,config);
    std::vector<ParticleType> species;
    species.emplace_back("negative",2.*PhysicalConstants::eV_J,-1.);
    species.emplace_back("Hot-ELECTRONS",2.*PhysicalConstants::eV_J,1.);
    species.emplace_back("neutral",1.,0.);species.emplace_back("empty",1.,1.);
    Particle p;p.velocity={10.,0.};p.w=0.;p.rngState=123;p.trackId=456;
    species[0].particle.push_back(p);
    p.velocity={1.,2.};p.vz=3.;p.w=2.;species[1].particle.push_back(p);
    p.velocity={0.,0.};p.vz=2.;p.w=0.;species[1].particle.push_back(p);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    p.w=nan;species[2].particle.push_back(p);
    const std::vector<double> weights{3.};
    const PicDomainBC bc;const PicEmbeddedBC eb;const PicEmbeddedBCMap byId;
    std::vector<double> kappaSegments(graph.segments().size(),2.),kappaCells(mesh.cellIndexNum(),3.);
    const auto sample=[&](const Solvers::Poisson::CutCellMetrics* metrics=nullptr) {
        auto result=Diagnostics::sampleExplicit(mesh,graph,bc,rhs,phi,charge,4.,species,weights,norm,metrics,eb,byId);
        if(implicit)result.values=Diagnostics::sampleImplicit(mesh,graph,bc,rhs,phi,charge,4.,species,weights,norm,
            metrics,&eb,&byId,kappaSegments,kappaCells);
        return result;
    };
    auto result=sample();const auto& values=result.values;
    check(result.phiMin==-2. && result.phiMax==5.,"explicit potential extrema remain in code units");
    check(values.speciesNames==std::vector<std::string>{"negative","Hot-ELECTRONS","neutral","empty"}
        && values.speciesMacroCount==std::vector<long long>{1,2,1,0}
        && values.speciesPhysCount==std::vector<double>{3.,3.,1.,0.}
        && values.macro_total==4 && values.phys_total==7.,"species order, empty species and missing/default weights");
    check(values.speciesNmean_1pm3==std::vector<double>{.25,.25,1./12.,0.},"density uses physical area and slab thickness");
    check(std::abs(values.electron_mean_energy_eV-128./3.)<1e-13
        && std::abs(values.Te_e_mean_eV-256./9.)<1e-13,"name precedes charge and energy includes weighted vz and velocity units");
    const auto gauss=Diagnostics::compute(mesh,graph,bc,rhs,phi,&charge,4.,nullptr,nullptr,&eb,&byId,
        implicit?&kappaSegments:nullptr,implicit?&kappaCells:nullptr);
    check(values.rG_L2==gauss.rG_L2 && values.rG_Linf==gauss.rG_Linf && values.energyW==gauss.energyW
        && values.Qgrid==gauss.Qgrid && values.dQ==gauss.dQ,"sample retains the assembled Gauss calculation");
    // Controlled non-const fixture: exercise legacy partial/invalid alpha input.
    Solvers::Poisson::CutCellMetrics metrics;
    auto& cells=const_cast<std::vector<Solvers::Poisson::CutCellMetrics::CellMetric>&>(metrics.cellMetricsByIndex());
    cells.resize(4);cells[0].alpha=nan;cells[1].alpha=-2.;cells[2].alpha=2.;cells[3].alpha=.5;
    check(sample(&metrics).values.speciesNmean_1pm3[0]==3./7.5,"nonfinite alpha fallback and finite alpha clipping");
    cells.resize(1);cells[0].alpha=.5;
    check(sample(&metrics).values.speciesNmean_1pm3[0]==3./10.5,"missing metric entries retain full volume");
    cells.resize(4);for(auto& cell:cells)cell.alpha=0.;
    check(sample(&metrics).values.speciesNmean_1pm3[0]==3./1e-300,"closed volume retains density denominator floor");
    species[1]=ParticleType("positive",2.*PhysicalConstants::eV_J,1.);
    check(sample().values.electron_mean_energy_eV==400.,"first negative species is the electron fallback");
    for(const auto& name:{"E","Electron","preElectronPost"}) {
        species[1]=ParticleType(name,1.,1.);
        check(sample().values.electron_mean_energy_eV==0.,"empty named electron still takes precedence");
    }
    check(species[0].particle[0].rngState==123 && species[0].particle[0].trackId==456
        && species[0].particle[0].velocity.x==10. && species[0].particle[0].w==0.,"diagnostic sampling leaves particles unchanged");
    // Controlled stale child slot: retain each sampler's existing volume test.
    auto& leaf=mesh.cells(mesh.cellIndexToLabel[0]);const int child=leaf.son[1];leaf.son[1]=0;
    const double density=sample().values.speciesNmean_1pm3[0];leaf.son[1]=child;
    check(density==(implicit?3./9.:3./12.),"implicit all-child and explicit first-child volume predicates remain distinct");
    species.clear();result=sample();
    check(result.values.speciesNames.empty() && result.values.phys_total==0.
        && result.values.electron_mean_energy_eV==0.,"empty species produce a fresh empty summary");
    phi.value[0]=nan;check(sample().phiMin==1. && sample().phiMax==5.,"NaN potential is ignored by the original min/max scan");
    for(int ci=0;ci<phi.size();++ci)phi.value[ci]=nan;
    result=sample();check(result.phiMin==std::numeric_limits<double>::infinity()
        && result.phiMax==-std::numeric_limits<double>::infinity(),"all-NaN potential retains sentinel extrema");
    Diagnostics::ExplicitSample report{};report.phiMin=-2.;report.phiMax=5.;
    report.values.rG_L2=4.;report.values.rG_Linf=3.;report.values.energyW=7.;
    report.values.Qgrid=999.; // Console uses deposition totals, not this Gauss sum.
    struct FlushBuffer:std::stringbuf {int flushes=0;int sync() override {++flushes;return 0;}} buffer;
    std::ostream stream(&buffer);PicStepSeries::writeExplicitDiagnostic(stream,8,2.,-1.,report,"solved");
    check(buffer.str()=="[ExplicitPIC][diag] step=8 Qgrid=2 Qmacro=-1 dQ=3 rel=-3 rG_L2=4 rG_Linf=3 W=7 phi[min,max]=[-2, 5] poisson: solved\n"
        && stream.precision()==16 && buffer.flushes==1,"explicit console exact format, charge source, precision and flush");
    std::ostringstream zero;PicStepSeries::writeExplicitDiagnostic(zero,0,2.,0.,report,{});
    check(zero.str().find("dQ=2 rel=0")!=std::string::npos && zero.str().find("poisson:")==std::string::npos,
        "zero macro charge and absent solver message");
}
void checkDomainOwnership() {
    static_assert(std::is_same_v<decltype(std::declval<const PicDomain&>().mesh()),const QuadtreeMesh&>);
    static_assert(std::is_same_v<decltype(std::declval<const PicDomain&>().faceGraph()),const FaceSegmentGraph&>);
    static_assert(std::is_same_v<decltype(std::declval<const PicDomain&>().refiner()),const QuadtreeMeshRefiner&>);
    for(auto mode:{PicDomain::Mode::Explicit,PicDomain::Mode::Implicit}) {
        auto makeDomain=[&]() {
            DictionaryMesh meshInput;meshInput.Nx=meshInput.Ny=2;meshInput.maxCellNum=64;
            meshInput.maxLevel_local={1};meshInput.minLevel_neg={0};meshInput.minLevel_pos={0};
            auto policy=meshInput;policy.geometryRefinement.enabled=true;
            DictionaryMesh::GeometryRefinementConfig::RegionRule rule;
            rule.name="upper";rule.shape=DictionaryMesh::GeometryRefinementConfig::RegionShape::kBox;
            rule.min={0.,.5};rule.max={1.,1.};rule.maxLevel=1;
            policy.geometryRefinement.regionRules={rule};
            auto domain=std::make_unique<PicDomain>(meshInput,policy,mode);
            policy.geometryRefinement.enabled=false;policy.geometryRefinement.regionRules.clear();
            return domain;
        };
        auto first=makeDomain(),second=makeDomain();
        auto& mesh=first->mesh();auto& graph=first->faceGraph();auto& refiner=first->refiner();
        check(refiner.attachGeometry(nullptr,nullptr),"refiner owns a stable policy after its input is changed and destroyed");
        CellField<double> target(mesh,0.);refiner.buildTargetLevels(mesh,target);
        bool hasRefinement=false;
        for(int ci=0;ci<target.size();++ci)hasRefinement|=target.value[ci]==1.;
        check(hasRefinement,"retained region policy produces refinement targets");
        check(graph.segments().empty(),"domain construction leaves graph building to its original stage");
        first->rebuildFaceGraph();const auto initial=graph.segments().size();
        check(initial>0 && graph.cellOffsets().size()==static_cast<std::size_t>(mesh.cellIndexNum()+1),"geometry-free domain includes all cells");
        DictionaryMesh::RefineBox box;box.min={0.,.5};box.max={1.,1.};box.level=1;
        mesh.applyRefineBoxes({box});first->rebuildFaceGraph();
        check(mesh.cellIndexNum()==10 && graph.segments().size()>initial && second->mesh().cellIndexNum()==4
            && second->faceGraph().segments().empty(),"topology and derived graph storage are independent between domains");
        check(&mesh==&first->mesh() && &graph==&first->faceGraph() && &refiner==&first->refiner(),"AMR preserves borrowed object identities");
        PicFieldState fields(mesh);PicPoissonFieldSolve poisson(fields,{});PicPoissonFieldSolve::Report report;
        poisson.solveExplicit({graph,PicDomainBC{}},0,report);
        check(report.success && poisson.hasCache(),"field services use the owned refined mesh and face graph");
    }
}
void checkMovingGeometry(const DictionaryMesh& input,
                         const DictionaryGeometry& source, PicMovingGeometry::Mode mode) {
    using Mode=PicMovingGeometry::Mode;
    for(bool erosion:{false,true}) {
        auto dict=source;dict.moving={true};
        for(auto& line:dict.lines)line.levelSetIndex=0;
        PicDomain ready(input,input,mode);auto& mesh=ready.mesh();auto& embedded=ready.embedded();bool eb=true,refine=false;
        check(ready.initializeGeometry(&dict,eb,refine,{}) && ready.geometry(),"moving geometry fixture initializes");
        ready.refineInitial(&dict,nullptr,eb,true);
        auto& geometry=*ready.geometry();
        MovingEBConfig config;config.enabled=true;config.updateEvery=2;config.outputEvery=1;
        const std::string prefix=std::string(mode==Mode::Explicit?"moving_explicit_":"moving_implicit_")+(erosion?"erosion":"rigid");
        config.metricsPath=prefix+"_metrics.csv";config.erosionMetricsPath=prefix+"_regions.csv";config.impactStatsPath=prefix+"_impacts.csv";
        MovingEBRegion region;region.name="wall";region.boundaryId=1;region.levelSetIndex=0;
        region.selector.wholeBoundary=true;region.motionMode=erosion?"normalErosion":"rigidShift";
        region.velocityX_mps=.01;region.normalSpeed_mps=.001;config.regions={region};
        PicMovingGeometry moving(mode);
        check(moving.initialize(mesh,&geometry,&dict,config,nullptr,nullptr,std::nullopt)
            && moving.controller() && !moving.updated(),"module owns an initialized wall controller");
        auto& graph=ready.faceGraph();FaceSegmentGraph allCells;allCells.build(mesh);
        auto& refiner=ready.refiner();check(refiner.attachGeometry(&dict,&geometry),"moving refiner attaches");
        PicFieldState solveFields(mesh);PicPoissonFieldSolve poisson(solveFields,{});ImplicitPicFieldSolve polarization(solveFields,{},PicConfig::PolarizationPoissonConfig{},false,false);
        const std::vector<int> lineMap(4,1);
        const auto cache=[&]() {PicPoissonFieldSolve::Report report;poisson.solveExplicit({allCells,PicDomainBC{}},0,report);check(report.success && poisson.hasCache(),"moving factor cache is populated");};
        const auto refresh=[&](const std::vector<int>& ids) {
            moving.refreshAfterUpdate(2,{mesh,geometry,embedded,graph,ids,true,&refiner,poisson,&polarization});
        };
        moving.beginStep();
        ExplicitPIC::Pusher::EmbeddedImpact impact;impact.pos={.3,.5};impact.cellLabel=mesh.findCellLabel(impact.pos);
        impact.particle.position=impact.pos;impact.particle.velocity={2.,0.};impact.crossing=true;
        moving.ingestImpact(mesh,geometry,nullptr,lineMap,impact,0,6.63e-26,1.,3.,1.,.1);
        moving.ingestImpact(mesh,geometry,nullptr,lineMap,impact,1,9.1e-31,-1.,3.,1.,.1);
        check(moving.controller()->resolvedImpactCount()==1 && moving.controller()->ignoredNonIonImpactCount()==1,
            "lazy impact lookup resolves ions and retains non-ion accounting");
        cache();
        check(!moving.advance(1,.1,.1,mesh,&geometry) && !moving.updated() && poisson.hasCache(),"wall cadence skip preserves factorization");
        moving.beginStep();check(moving.controller()->resolvedImpactCount()==0,"begin step clears impact counters");
        check(moving.advance(2,.2,.1,mesh,&geometry) && moving.updated(),"wall update completes at configured cadence");
        refresh(lineMap);
        check(!poisson.hasCache() && embedded.metricsReady(),"moved geometry publishes metrics and invalidates matrix");
        check(moving.controller()->linesInCurrentGeometryFrame()!=erosion
            && moving.controller()->lastUpdateUsedAdvectedPath()==erosion,"rigid and advected geometry frames remain distinct");
        check(std::abs(geometry.lines[0].points[0].x-(erosion?.3:.302))<1e-12,"rigid lines move by the whole update-window duration");
        moving.beginStep();check(moving.updated(),"begin step preserves last update for existing diagnostic timing");
        cache();check(!moving.advance(3,.3,.1,mesh,&geometry),"next noncadence step reports unchanged geometry");
        refresh(lineMap);check(poisson.hasCache(),"refresh is inactive after an unchanged update");
        const std::filesystem::path saved=prefix+"_cache";moving.writeCache(saved);
        PicMovingGeometry restored(mode);
        check(restored.initialize(mesh,&geometry,&dict,config,nullptr,nullptr,saved)
            && !restored.updated()
            && restored.controller()->linesInCurrentGeometryFrame()==moving.controller()->linesInCurrentGeometryFrame(),"checkpoint restoration retains geometry-frame state without reporting a new update");
        check(std::abs(restored.controller()->boundaryRigidShiftCode(1).x-(erosion?0.:.002))<1e-12,"restored controller retains accumulated rigid displacement");
        PicMovingGeometry missing(mode);
        check(!missing.initialize(mesh,&geometry,&dict,config,nullptr,nullptr,prefix+"_missing"),"required missing MovingEB cache fails initialization");
        config.restart=false;PicMovingGeometry noRestart(mode);
        check(noRestart.initialize(mesh,&geometry,&dict,config,nullptr,nullptr,prefix+"_missing"),"restart disabled ignores missing wall cache");
        noRestart.writeCache(prefix+"_disabled_cache");
        check(!std::filesystem::exists(prefix+"_disabled_cache"),"restart disabled does not create wall checkpoint");
        if(erosion && mode==Mode::Implicit) {
            moving.beginStep();check(moving.advance(4,.4,.1,mesh,&geometry),"erosion failure fixture advances");
            bool rejected=false;
            try { refresh(std::vector<int>(4,0)); }
            catch(const std::runtime_error& e) {rejected=std::string(e.what()).find("invalid cut metrics after erosion update")!=std::string::npos;}
            check(rejected && moving.updated(),"invalid erosion boundary tags retain the original failure-stop behavior");
        }
    }
}
void checkDomainStartup(const DictionaryMesh& input, DictionaryGeometry dict, PicDomain::Mode mode) {
    const bool implicit=mode==PicDomain::Mode::Implicit;
    const std::string dir=implicit?"domain_startup_implicit/":"domain_startup_explicit/";
    std::filesystem::create_directory(dir);
    const std::vector<int> lineMap(4,1);
    const auto fileBytes=[](const std::string& path){std::ifstream f(path,std::ios::binary);
        check(bool(f),"domain cache missing");return std::string(std::istreambuf_iterator<char>(f),{});};
    const auto rejected=[](const auto& action){try{action();}catch(const std::logic_error&){return true;}return false;};
    // Independent reference: the driver's original operation order, including
    // mode-specific phase/composite order and refinement after explicit restore.
    const auto original=[&](PicDomain& domain,bool changes,const PicDomain::Restore& restore,bool skipImplicitRefinement) {
        bool eb=true,refine=true;
        check(domain.initializeGeometry(&dict,eb,refine,restore),"manual startup initializes geometry");
        auto& mesh=domain.mesh();auto& embedded=domain.embedded();auto* geometry=domain.geometry();
        bool active=false;
        const auto attach=[&](){active=domain.refiner().attachGeometry(&dict,geometry);
            if(active && domain.refiner().phaseActive() && geometry)geometry->markPhases(mesh);};
        if(!implicit)attach();
        if(embedded.ready())embedded.rebuildComposite(mesh,geometry);
        if(implicit)attach();
        if((!implicit || !skipImplicitRefinement) && active && geometry && mesh.maxInitialIterations()>0)
            domain.refineInitial(&dict,nullptr,eb,domain.refiner().phaseActive());
        if(active && !changes){domain.refiner().clearGeometry();active=false;}
        domain.rebuildFaceGraph();
        if((!implicit || eb) && embedded.ready())embedded.rebuildMetrics(mesh,geometry,domain.faceGraph(),lineMap,nullptr);
        return PicDomain::InitialState{eb,active};
    };
    const auto same=[&](PicDomain& a,PicDomain& b) {
        for(auto item:{std::make_pair(&a,"actual/"),std::make_pair(&b,"reference/")}) {
            const auto prefix=dir+item.second;std::filesystem::create_directory(prefix);
            item.first->mesh().writeCache(prefix);item.first->writeGeometryCache(prefix);
        }
        for(const auto& name:{"mesh_cache.bin","iGeometry_cache.bin"})
            check(fileBytes(dir+"actual/"+name)==fileBytes(dir+"reference/"+name),"domain startup preserves exact mesh/geometry caches");
        check(a.faceGraph().cellOffsets()==b.faceGraph().cellOffsets()
            && a.faceGraph().segments().size()==b.faceGraph().segments().size(),"startup face graph topology");
        for(std::size_t k=0;k<a.faceGraph().segments().size();++k) {
            const auto& x=a.faceGraph().segments()[k];const auto& y=b.faceGraph().segments()[k];
            check(x.K==y.K && x.L==y.L && x.boundary==y.boundary && x.side==y.side && x.isHorizontal==y.isHorizontal
                && x.edgeLabel==y.edgeLabel && x.Af==y.Af && x.n.x==y.n.x && x.n.y==y.n.y && x.df==y.df
                && x.faceCenter.x==y.faceCenter.x && x.faceCenter.y==y.faceCenter.y && x.edgeCut==y.edgeCut
                && x.edgeLength==y.edgeLength && x.edgeCutLength==y.edgeCutLength,"startup oriented face geometry");
        }
        const auto& x=a.embedded();const auto& y=b.embedded();
        check(x.alphaByCell()==y.alphaByCell() && x.boundaryIdsByCell()==y.boundaryIdsByCell()
            && x.metrics().betaBySegment()==y.metrics().betaBySegment()
            && x.metrics().activeMaskByCellIndex()==y.metrics().activeMaskByCellIndex(),"startup alpha, beta and boundary mapping");
        const auto& xm=x.metrics().cellMetricsByIndex();const auto& ym=y.metrics().cellMetricsByIndex();
        check(xm.size()==ym.size(),"startup cut metric size");
        for(std::size_t k=0;k<xm.size();++k)
            check(xm[k].alpha==ym[k].alpha && xm[k].A_gamma==ym[k].A_gamma && xm[k].d_gamma==ym[k].d_gamma
                && xm[k].n_gamma.x==ym[k].n_gamma.x && xm[k].n_gamma.y==ym[k].n_gamma.y && xm[k].isCut==ym[k].isCut
                && xm[k].gammaBoundaryId==ym[k].gammaBoundaryId && xm[k].gammaLineIndex==ym[k].gammaLineIndex,"startup full cut metric values");
    };
    for(bool changes:{false,true}) {
        PicDomain actual(input,input,mode),reference(input,input,mode);
        check(rejected([&]{actual.initializeForRun(&dict,true,true,changes,lineMap);}) && !actual.geometry(),"run startup requires mesh preparation before mutation");
        check(actual.prepareMesh(nullptr),"fresh mesh preparation");
        check(rejected([&]{actual.prepareMesh(nullptr);}),"mesh preparation cannot repeat");
        const auto* mesh=&actual.mesh();const auto* graph=&actual.faceGraph();
        auto result=actual.initializeForRun(&dict,true,true,changes,lineMap);
        const auto expected=original(reference,changes,{},false);
        check(result && result->embeddedEnabled==expected.embeddedEnabled && result->refinementActive==expected.refinementActive
            && result->refinementActive==changes,"startup returns runtime flags after static refinement removal");
        check(actual.mesh().cellIndexNum()>input.Nx*input.Ny && actual.embedded().metricsReady(),"startup performs real refinement before derived metrics");
        same(actual,reference);
        const auto* geometry=actual.geometry();
        check(rejected([&]{actual.initializeForRun(nullptr,false,false,false,{});})
            && rejected([&]{actual.prepareMesh(nullptr);}) && geometry==actual.geometry()
            && mesh==&actual.mesh() && graph==&actual.faceGraph(),"repeated startup is rejected without replacing live owners");
    }
    // Restore a deliberately coarse checkpoint: explicit still refines it;
    // implicit must keep its restored topology, even with startup AMR requested.
    auto coarseInput=input;coarseInput.maxInitialIterations=0;
    PicDomain coarse(coarseInput,coarseInput,mode);original(coarse,true,{},false);
    const std::string checkpoint=dir+"pic.2/";std::filesystem::create_directory(checkpoint);
    coarse.mesh().writeCache(checkpoint);coarse.writeGeometryCache(checkpoint);
    check(PicCheckpointStore::writeMetadata(checkpoint,2,.2,.1),"startup checkpoint metadata fixture");
    PicCheckpointStore store(dir,{});double dt=.1;
    auto reader=PicCheckpointReader::select(store,2,dt,implicit?"ImplicitPIC":"ExplicitPIC");check(reader.has_value(),"startup checkpoint selected");
    PicDomain restored(input,input,mode),reference(input,input,mode);
    check(restored.prepareMesh(&*reader),"run mesh restores checkpoint");
    auto result=restored.initializeForRun(&dict,true,true,true,lineMap);check(result.has_value(),"run geometry restores checkpoint");
    reference.mesh().readCache(checkpoint);PicDomain::Restore restore;restore.directory=checkpoint;restore.nextStep=2;
    original(reference,true,restore,true);same(restored,reference);
    check(implicit?restored.mesh().cellIndexNum()==input.Nx*input.Ny:restored.mesh().cellIndexNum()>input.Nx*input.Ny,
        "explicit/implicit startup refinement after checkpoint remains distinct");
    std::filesystem::create_directory(dir+"pic.3");check(PicCheckpointStore::writeMetadata(dir+"pic.3",3,.3,.1),"missing-mesh fixture");
    auto missingReader=PicCheckpointReader::select(store,3,dt,"PIC");check(missingReader.has_value(),"explicit checkpoint selection retains legacy component policy");
    PicDomain missing(input,input,mode);
    check(missing.prepareMesh(&*missingReader)==!implicit,"explicit optional and implicit required mesh behavior");
    if(implicit)check(rejected([&]{missing.initializeForRun(&dict,true,true,false,lineMap);})
        && rejected([&]{missing.prepareMesh(nullptr);}) && !missing.geometry(),"failed required restore cannot initialize or retry");
    else check(missing.initializeForRun(&dict,true,true,false,lineMap).has_value(),"optional missing mesh initializes original mesh");
    PicDomain absent(input,input,mode);check(absent.prepareMesh(nullptr),"no geometry mesh preparation");
    auto empty=absent.initializeForRun(nullptr,true,true,false,{});
    check(empty && !empty->embeddedEnabled && !empty->refinementActive && !absent.geometry()
        && !absent.faceGraph().segments().empty(),"missing geometry disables geometry consumers and builds plain graph");
    if(implicit) {
        const auto meshCache=fileBytes(checkpoint+"mesh_cache.bin"),geometryCache=fileBytes(checkpoint+"iGeometry_cache.bin");
        PicDomain snapshot(input,input,mode);check(snapshot.prepareMesh(nullptr,checkpoint),"initial snapshot loads mesh");
        check(snapshot.initializeForRun(&dict,true,true,true,lineMap).has_value() && snapshot.mesh().cellIndexNum()==input.Nx*input.Ny,
            "loaded snapshot skips initial refinement");
        same(snapshot,coarse);
        check(fileBytes(checkpoint+"mesh_cache.bin")==meshCache && fileBytes(checkpoint+"iGeometry_cache.bin")==geometryCache,"loaded snapshot remains untouched");
        PicDomain created(input,input,mode);check(created.prepareMesh(nullptr,dir+"new_snapshot"),"missing snapshot falls back to fresh mesh");
        check(created.initializeForRun(&dict,true,true,false,lineMap).has_value(),"fresh initialization creates requested snapshot");
        check(std::filesystem::exists(dir+"new_snapshot/mesh_cache.bin") && std::filesystem::exists(dir+"new_snapshot/iGeometry_cache.bin"),"both initial snapshot components saved");
        PicDomain preferred(input,input,mode);check(preferred.prepareMesh(&*reader,dir+"new_snapshot"),"checkpoint takes priority over available snapshot");
        check(preferred.initializeForRun(&dict,true,true,true,lineMap).has_value() && preferred.mesh().cellIndexNum()==input.Nx*input.Ny,"checkpoint topology wins over refined snapshot");
        {std::ofstream blocker(dir+"blocked");blocker<<"preserve";}
        PicDomain blocked(input,input,mode);check(blocked.prepareMesh(nullptr,dir+"blocked/snapshot"),"unavailable snapshot falls back");
        check(blocked.initializeForRun(&dict,true,true,false,lineMap).has_value() && fileBytes(dir+"blocked")=="preserve","snapshot create failure retains original nonfatal behavior");
    } else {
        PicDomain explicitOnly(input,input,mode);
        check(rejected([&]{explicitOnly.prepareMesh(nullptr,checkpoint);}) && explicitOnly.prepareMesh(nullptr),"unsupported explicit snapshot rejected before preparation");
    }
}
void checkGeometrySetup() {
    using Mode=PicDomain::Mode;
    static_assert(!std::is_copy_constructible_v<PicDomain> && !std::is_move_constructible_v<PicDomain>);
    static_assert(std::is_same_v<decltype(std::declval<const PicDomain&>().geometry()),const PicDomain::Geometry*>);
    static_assert(std::is_same_v<decltype(std::declval<const PicDomain&>().embedded()),const PicEmbeddedGeometry&>);
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=256;
    input.maxInitialIterations=2;input.maxLevel_local={1};input.minLevel_neg={0};input.minLevel_pos={0};
    input.geometryRefinement.enabled=true;input.geometryRefinement.present=true;
    input.geometryRefinement.boundaryRules.push_back({1,1.,1});
    DictionaryGeometry dict;dict.levelSetName={"solid"};dict.innerPoints={{.5,.5}};dict.moving={false};
    dict.linesOnInterface={{{0,false},{1,false},{2,false},{3,false}}};
    dict.phaseName={"solid","gas"};dict.levelSetOperation={{LevelSetOperation::non},{LevelSetOperation::non}};
    dict.levelSetIndexA={{{0,false}},{{0,true}}};dict.levelSetIndexB=dict.levelSetIndexA;
    dict.boundaryLines={{{0,false},{1,false},{2,false},{3,false}},{{0,true},{1,true},{2,true},{3,true}}};
    const std::array<Vector2D,4> corners{{{.3,.3},{.7,.3},{.7,.7},{.3,.7}}};
    for(int i=0;i<4;++i) {Line line;line.points[0]=corners[i];line.points[1]=corners[(i+1)%4];dict.lines.push_back(line);}
    DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.kind=GeometricBoundary::fixed;
    boundary.levelSetIndex=0;boundary.phase=0;boundary.lineIds={0,1,2,3};dict.boundaries={boundary};dict.boundaryById[1]=boundary;
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        checkDomainStartup(input,dict,mode);
        bool eb=true,refine=true;
        PicDomain missing(input,input,mode);
        check(missing.initializeGeometry(nullptr,eb,refine,{}) && !missing.geometry() && !eb && !refine,"missing dictionary disables both requested geometry consumers");
        PicDomain unused(input,input,mode);
        check(unused.initializeGeometry(&dict,eb,refine,{}) && !unused.geometry(),"unused geometry is not constructed");
        DictionaryGeometry empty;eb=true;refine=true;
        PicDomain noFields(input,input,mode);check(noFields.initializeGeometry(&empty,eb,refine,{}),"empty geometry initializes");
        check(!noFields.embedded().ready() && (mode==Mode::Explicit
            ? !noFields.geometry() && !eb && !refine : noFields.geometry() && eb && refine),"empty level-set ownership and flags retain mode policy");

        eb=refine=true;
        PicDomain::Restore fallback;fallback.directory="geometry_setup_missing_cache";
        PicDomain ready(input,input,mode);auto& mesh=ready.mesh();auto& embedded=ready.embedded();
        check(ready.initializeGeometry(&dict,eb,refine,fallback) && ready.geometry() && embedded.ready() && eb && refine && !ready.hasScheduledMotion(),"missing optional cache initializes ordinary geometry");
        const auto& initial=ready.geometry()->levelSets[0];
        const auto* meshIdentity=&mesh;
        const auto* rawIdentity=ready.geometry();
        const auto* embeddedIdentity=&ready.embedded();
        for(int vi=0;vi<mesh.vertIndexNum();++vi)check(std::isfinite(initial.value[vi]),"initialized geometry has finite values");
        embedded.rebuildComposite(mesh,ready.geometry());
        const int before=mesh.cellIndexNum();
        ready.refineInitial(&dict,nullptr,eb,true);
        ready.rebuildFaceGraph();
        check(&ready.mesh()==meshIdentity && ready.faceGraph().cellOffsets().size()==static_cast<std::size_t>(mesh.cellIndexNum()+1),
            "geometry refinement rebuilds the owned graph without replacing its mesh");
        for(const auto& segment:ready.faceGraph().segments())
            check(embedded.activeCell(mesh,ready.geometry(),mesh.cells(segment.K)),"domain graph includes only active plasma owners");
        check(mesh.cellIndexNum()>before && ready.geometry()->levelSets[0].size()==mesh.vertIndexNum()
            && embedded.levelSet()==&ready.geometry()->levelSets[0],"startup refinement transfers and rebuilds the selected geometry on a changed mesh");
        if(mode==Mode::Explicit)check(!std::filesystem::exists("output/meta/init_geometry_refine.csv"),"explicit refinement does not create implicit CSV");
        else check(std::filesystem::file_size("output/meta/init_geometry_refine.csv")>200,"implicit refinement writes its original iteration table");
        const int refined=mesh.cellIndexNum();
        ready.refineInitial(&dict,nullptr,eb,true);
        check(mesh.cellIndexNum()==refined,"already converged startup mesh stays unchanged");
        const auto* selected=embedded.levelSet();
        const double retained=selected->value[0];
        bool reinitialized=false;
        try { ready.initializeGeometry(nullptr,eb,refine,{}); }
        catch(const std::logic_error&) { reinitialized=true; }
        check(reinitialized && ready.geometry()==rawIdentity && &ready.embedded()==embeddedIdentity
            && embedded.levelSet()==selected && selected->value[0]==retained && eb && refine,
            "reinitialization must not destroy geometry, derived views or caller flags");

        PicAmrGeometry::Options refreshOptions;refreshOptions.mode=PicAmrGeometry::RefreshMode::Skip;
        PicAmrGeometry refresh(mode,refreshOptions);
        refreshOptions.mode=PicAmrGeometry::RefreshMode::Full;
        const auto& report=refresh.report();
        PicAmrGeometry::Adaptation adaptation{3,.2,true,true,true,false};
        refresh.captureBeforeAdapt(ready.geometry());
        for(int vi=0;vi<mesh.vertIndexNum();++vi)ready.geometry()->LS_dx(0).value[vi]=123.;
        refresh.refresh(mesh,ready.geometry(),embedded,nullptr,adaptation);
        check(embedded.levelSet()==&ready.geometry()->levelSets[0],"AMR refresh publishes the selected field after geometry work");
        if(mode==Mode::Implicit) {
            check(report.geometryReinitSkipped && ready.geometry()->LS_dx(0).value[0]==123.,"owned skip policy leaves raw derivatives untouched");
            refresh.captureBeforeAdapt(ready.geometry());adaptation.movingUpdated=true;
            refresh.refresh(mesh,ready.geometry(),embedded,nullptr,adaptation);
            check(!report.geometryReinitSkipped && ready.geometry()->LS_dx(0).value[0]!=123.,"an updated wall forces the original full refresh even under skip policy");
            PicAmrGeometry::Options adaptive;adaptive.mode=PicAmrGeometry::RefreshMode::Adaptive;
            adaptive.amrGeomRefreshChangeRelTh=.2;PicAmrGeometry threshold(mode,adaptive);
            adaptation.movingUpdated=false;
            threshold.captureBeforeAdapt(ready.geometry());threshold.refresh(mesh,ready.geometry(),embedded,nullptr,adaptation);
            check(threshold.report().geometryReinitSkipped,"adaptive threshold includes equality");
            adaptation.cellCountChange=.20001;
            threshold.captureBeforeAdapt(ready.geometry());threshold.refresh(mesh,ready.geometry(),embedded,nullptr,adaptation);
            check(!threshold.report().geometryReinitSkipped,"adaptive threshold refreshes strictly above the limit");
        } else {
            check(!report.geometryReinitSkipped && ready.geometry()->LS_dx(0).value[0]!=123.,"explicit AMR retains full refresh independent of implicit skip policy");
        }
        refresh.captureBeforeAdapt(nullptr);refresh.refresh(mesh,nullptr,embedded,nullptr,adaptation);
        check(report.geometryRefreshS==0. && !report.geometryReinitSkipped && report.amrSelectiveMovingLevelSets.empty(),"unavailable geometry publishes no stale AMR report");

        if(mode==Mode::Implicit) {
            LegacyGridMotionConfig config;config.enabled=true;config.updateEvery=3;
            config.startupHoldSteps=1;config.stage1Steps=4;config.accelShiftStage1_m=.04;
            LegacyGridMotion::Bindings lines;lines.accelUpperLineIds={0,1,2,3};
            auto motion=std::make_unique<LegacyGridMotionEvolution>(config,lines,
                LegacyGridMotionEvolution::Reference{dict.lines,dict.innerPoints},
                LegacyGridMotionEvolution::InnerPointBindings{{0},{}});
            auto& graph=ready.faceGraph();auto& refiner=ready.refiner();
            check(refiner.attachGeometry(&dict,ready.geometry()),"scheduled geometry refiner attaches");
            PicFieldState solveFields(mesh);PicPoissonFieldSolve poisson(solveFields,{});ImplicitPicFieldSolve polarization(solveFields,{},PicConfig::PolarizationPoissonConfig{},false,false);
            const std::vector<int> lineMap(4,1);
            PicScheduledGeometry::Bindings bindings{mesh,ready.geometry(),embedded,graph,lineMap,nullptr,refiner,poisson,polarization};
            std::vector<std::pair<int,std::string>> events;
            PicScheduledGeometry scheduled(std::move(motion),bindings,{true,true,true,
                [&](int step,const std::string& message){events.emplace_back(step,message);}});
            check(!motion && scheduled.hasSchedule(),"schedule ownership transfers to geometry service");
            check(scheduled.advance(0,0.,true) && embedded.metricsReady()
                && scheduled.current().stageName=="stage0_hold","forced initialization rebuilds geometry and cut metrics");
            check(events.front().first==0 && events.front().second=="applyScheduledGeometry begin force=1 stage=stage0_hold changed=1"
                && events.back().second=="applyScheduledGeometry end","motion observer receives original begin/end events synchronously");
            const auto checkShift=[&](double shift) {
                check(std::abs(ready.geometry()->lines[0].points[0].x-(.3+shift))<1e-14
                    && std::abs(ready.geometry()->levelSets[0].innerPoint().x-(.5+shift))<1e-14
                    && std::abs(scheduled.current().accelShift_m-shift)<1e-14,"lines, inner point and published motion agree");
            };
            check(scheduled.advance(1,1.,false),"stage transition advances between cadence hits");checkShift(.01);
            // Populate a real factor cache: skipped updates retain it, whereas
            // changed and forced updates must invalidate the old matrix.
            FaceSegmentGraph allCells;allCells.build(mesh);
            const auto cache=[&]() {PicPoissonFieldSolve::Report report;poisson.solveExplicit({allCells,PicDomainBC{}},0,report);check(report.success && poisson.hasCache(),"scheduled cache fixture factors successfully");};
            cache();
            const auto eventsBeforeSkip=events.size();
            check(!scheduled.advance(2,2.,false) && !scheduled.current().changed && poisson.hasCache(),"cadence skip retains factorization and publishes unchanged state");checkShift(.01);
            check(events.size()==eventsBeforeSkip+1 && events.back().first==2
                && events.back().second=="applyScheduledGeometry begin force=0 stage=stage1_hot_gap changed=0","skipped motion emits only its original begin event");
            check(scheduled.advance(3,3.,false) && !poisson.hasCache(),"cadence update invalidates factorization");checkShift(.03);
            std::vector<double> values;
            for(int vi=0;vi<mesh.vertIndexNum();++vi)values.push_back(ready.geometry()->levelSets[0].value[vi]);
            cache();
            check(scheduled.advance(3,3.,true) && !scheduled.current().changed && !poisson.hasCache(),"force rebuilds even with unchanged accepted motion");checkShift(.03);
            for(int vi=0;vi<mesh.vertIndexNum();++vi)check(ready.geometry()->levelSets[0].value[vi]==values[vi],"repeated absolute motion does not accumulate deformation");
            PicScheduledGeometry absent(nullptr,bindings,{true,true,true});
            check(!absent.hasSchedule() && !absent.advance(4,4.,true)
                && absent.current().stageName=="static","missing schedule publishes static state without rebuilding");
            auto rejectedMotion=std::make_unique<LegacyGridMotionEvolution>(config,lines,
                LegacyGridMotionEvolution::Reference{dict.lines,dict.innerPoints},
                LegacyGridMotionEvolution::InnerPointBindings{{0},{}});
            PicScheduledGeometry rejecting(std::move(rejectedMotion),bindings,{true,true,true,
                [](int,const std::string&){throw std::runtime_error("observer stopped");}});
            const auto previousState=rejecting.current();const double xBefore=ready.geometry()->lines[0].points[0].x;
            bool observerRejected=false;try {rejecting.advance(4,4.,false);}catch(const std::runtime_error& e){observerRejected=std::string(e.what())=="observer stopped";}
            check(observerRejected && rejecting.current().stageName==previousState.stageName
                && ready.geometry()->lines[0].points[0].x==xBefore,"observer exception propagates at its event before deformation or acceptance");
        }

        const std::string cache=mode==Mode::Explicit?"geometry_setup_explicit_cache/":"geometry_setup_implicit_cache/";
        std::filesystem::create_directory(cache);
        for(int vi=0;vi<mesh.vertIndexNum();++vi)ready.geometry()->levelSets[0].value[vi]=42.;
        mesh.writeCache(cache);ready.writeGeometryCache(cache);
        PicDomain::Restore restore;restore.directory=cache;
        PicDomain restored(input,input,mode);restored.mesh().readCache(cache);
        check(restored.initializeGeometry(&dict,eb,refine,restore),"geometry owner restores existing cache");
        for(int vi=0;vi<mesh.vertIndexNum();++vi)check(restored.geometry()->levelSets[0].value[vi]==42.,"existing geometry cache takes precedence over fresh initialized fields");
        ready.geometry()->levelSets[0].value[0]=43.;
        check(restored.geometry()!=ready.geometry() && restored.geometry()->levelSets[0].value[0]==42.,
            "separate geometry owners must not share mutable level-set storage");
        check(ready.geometry()==rawIdentity && &ready.embedded()==embeddedIdentity,
            "AMR and prescribed motion must preserve the owning geometry identities");
        if(mode==Mode::Implicit) {
            LegacyGridMotionConfig schedule;schedule.enabled=true;
            PicDomain pending(input,input,mode);
            check(pending.initializeGeometry(&dict,eb,refine,{},&schedule) && pending.hasScheduledMotion(),
                "geometry owner retains startup motion until the runtime service is bound");
            const auto* original=pending.geometry();
            auto runtimeMotion=pending.takeScheduledMotion();
            check(runtimeMotion && !pending.hasScheduledMotion() && pending.geometry()==original
                && pending.embedded().ready(),"motion handoff must not transfer or destroy geometry");
        }
        checkMovingGeometry(input,dict,mode);
    }
}
void checkEmbeddedGeometry() {
    using Embedded=PicEmbeddedGeometry;
    DictionaryMesh input;input.Nx=8;input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    DictionaryGeometry dict;dict.levelSetName={"leftSolid","rightSolid"};dict.innerPoints={{.1,.5},{.9,.5}};
    dict.moving={false,false};dict.linesOnInterface={{{0,false}},{{1,true}}};
    dict.phaseName={"left","right"};dict.levelSetOperation={{LevelSetOperation::non},{LevelSetOperation::non}};
    dict.levelSetIndexA={{{0,false}},{{1,false}}};dict.levelSetIndexB={{},{}};dict.boundaryLines={{{0,false}},{{1,true}}};
    for(int k=0;k<2;++k) {
        Line line;line.points[0]={k?.7:.3,0.};line.points[1]={k?.7:.3,1.};dict.lines.push_back(line);
        DictionaryGeometry::BoundaryEntity entity;entity.id=k+1;entity.kind=GeometricBoundary::fixed;
        entity.levelSetIndex=k;entity.phase=k;entity.lineIds={k};dict.boundaries.push_back(entity);dict.boundaryById[k+1]=entity;
    }
    ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh);
    const std::vector<int> lineIds{1,2};
    const auto fill=[&](bool positive=false) {
        for(int k=0;k<2;++k) {
            geometry.levelSets[k].initialize(mesh,0.);geometry.LS_dx(k).initialize(mesh,0.);geometry.LS_dy(k).initialize(mesh,0.);
            for(int vi=0;vi<mesh.vertIndexNum();++vi) {
                const double x=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x;
                geometry.levelSets[k].value[vi]=(positive?-1.:1.)*(k?.7-x:x-.3);
                geometry.LS_dx(k).value[vi]=(positive?-1.:1.)*(k?-1.:1.);
            }
        }
    };
    fill();
    Embedded explicitGeometry(Embedded::Mode::Explicit),implicitGeometry(Embedded::Mode::Implicit);
    check(!explicitGeometry.ready() && !explicitGeometry.metricsReady() && !explicitGeometry.levelSet(),"uninitialized embedded state");
    for(auto* state:{&explicitGeometry,&implicitGeometry}) {
        state->resolveSelection(dict,2);state->rebuildComposite(mesh,&geometry);
        check(state->ready() && state->absorbing() && state->insideIsNegative() && state->usesComposite()
            && state->levelSetIndices()==std::vector<int>({0,1}) && state->boundaryIds()==lineIds,"selected composite ownership");
        check(state->compositeLines().size()==2 && state->compositeLines()[0].value==0
            && state->compositeLines()[1].value==1 && state->compositeLines()[1].isNegative,"line union preserves original signed entries");
        for(int vi=0;vi<mesh.vertIndexNum();++vi) {
            const double x=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x;
            check(state->levelSet()->value[vi]==std::min(x-.3,.7-x),"negative-solid union chooses minimum field");
        }
        check(state->boundaryIdAt(mesh,&geometry,lineIds,{.2,.4},-1)==1
            && state->boundaryIdAt(mesh,&geometry,lineIds,{.8,.4},-1)==2,"nearest polyline fallback before cut metrics");
    }
    const Cell& nearSolid=mesh.cells(mesh.findCellLabel({.2,.25}));
    check(!explicitGeometry.activeCell(mesh,&geometry,nearSolid) && implicitGeometry.activeCell(mesh,&geometry,nearSolid),
        "implicit solid margin keeps a shallow all-inside cell active");
    check(explicitGeometry.containsPoint(mesh,&geometry,{.30001,.25},-1)
        && !implicitGeometry.containsPoint(mesh,&geometry,{.30001,.25},-1),"explicit mesh-scaled versus implicit fixed point tolerance");

    // Make the level sets exactly equal beside a cut cell, with conflicting
    // gradients. The first source must still produce the analytic cut normal.
    int tiedVertices=0;
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        if(mesh.verts(mesh.vertIndexToLabel[vi]).coord().x!=.375)continue;
        geometry.levelSets[1].value[vi]=geometry.levelSets[0].value[vi];
        geometry.LS_dx(1).value[vi]=-3.;++tiedVertices;
    }
    check(tiedVertices>0,"exact tie fixture reaches cut-cell vertices");
    for(auto* state:{&explicitGeometry,&implicitGeometry})state->rebuildComposite(mesh,&geometry);

    FaceSegmentGraph graph;graph.build(mesh);
    VertField<double> psi(mesh,0.),dx(mesh,0.),dy(mesh,0.);
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        const double x=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x;
        const bool first=x-.3<=.7-x;psi.value[vi]=first?x-.3:.7-x;dx.value[vi]=first?1.:-1.;
    }
    Embedded::Metrics reference;Embedded::Metrics::BuildOptions options;
    options.levelSetI=-1;options.insideIsNegative=false;options.alphaSampleN=6;
    options.psiOverride=&psi;options.psiDxOverride=&dx;options.psiDyOverride=&dy;
    options.lineToBoundaryId=&lineIds;options.linesOverride=&explicitGeometry.compositeLines();
    reference.build(mesh,geometry,graph,options);
    for(auto* state:{&explicitGeometry,&implicitGeometry}) {
        state->rebuildMetrics(mesh,&geometry,graph,lineIds);
        check(state->metricsReady() && state->alphaByCell().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"owned cut cache built");
        for(int ci=0;ci<mesh.cellIndexNum();++ci) {
            const auto& actual=state->metrics().cellMetricsByIndex()[ci];const auto& expected=reference.cellMetricsByIndex()[ci];
            check(actual.alpha==expected.alpha && actual.A_gamma==expected.A_gamma && actual.d_gamma==expected.d_gamma
                && actual.n_gamma.x==expected.n_gamma.x && actual.n_gamma.y==expected.n_gamma.y
                && actual.gammaBoundaryId==expected.gammaBoundaryId && actual.gammaLineIndex==expected.gammaLineIndex
                && actual.isCut==expected.isCut,"cut metrics match analytic composite field and first-winner gradient");
            check(state->alphaByCell()[ci]==(reference.isActiveCellIndex(ci)?expected.alpha:0.)
                && state->boundaryIdsByCell()[ci]==expected.gammaBoundaryId,"published physical alpha and boundary labels");
        }
        check(state->metrics().betaBySegment()==reference.betaBySegment(),"composite face openings");
        const int staleHint=mesh.findCellLabel({.1,.8});
        check(state->boundaryIdAt(mesh,&geometry,lineIds,{.7,.2},staleHint)==2,"cut-boundary lookup repairs a stale cell hint");
    }
    UnitsScales scales;scales.L0_m=2.;auto units=std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));
    NormalizationConfig config;config.slabThickness_m=3.;auto norm=NormalizationContext::MakeFromUnits(units,config);
    check(explicitGeometry.surfaceAreaM2(&geometry,dict,lineIds,units.get(),&norm)==12.
        && explicitGeometry.surfaceAreaM2(nullptr,dict,lineIds,units.get(),&norm)==12.
        && explicitGeometry.surfaceAreaM2(&geometry,dict,lineIds,units.get(),nullptr)==0.,"surface area units, dictionary fallback and normalization gate");
    const double owned=explicitGeometry.levelSet()->value[0];
    for(int k=0;k<2;++k)for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry.levelSets[k].value[vi]=-1e-6;
    check(explicitGeometry.levelSet()->value[0]==owned,"composite snapshot owns values");
    check(!explicitGeometry.containsPoint(mesh,&geometry,{.45,.25},-1)
        && implicitGeometry.containsPoint(mesh,&geometry,{.45,.25},-1),"explicit cached full-fluid alpha short-circuit remains distinct");
    geometry.levelSets[0].value[0]=std::numeric_limits<double>::quiet_NaN();geometry.levelSets[1].value[0]=std::numeric_limits<double>::infinity();
    explicitGeometry.rebuildComposite(mesh,&geometry);
    check(explicitGeometry.levelSet()->value[0]==0.,"composite all-nonfinite samples retain zero initialization");

    auto positive=dict;positive.levelSetIndexA={{{0,true}},{{1,true}}};fill(true);
    Embedded positiveGeometry(Embedded::Mode::Implicit);positiveGeometry.resolveSelection(positive,2);positiveGeometry.rebuildComposite(mesh,&geometry);
    check(!positiveGeometry.insideIsNegative(),"positive-solid selection inferred from phase references");
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        const double x=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x;
        check(positiveGeometry.levelSet()->value[vi]==std::max(.3-x,x-.7),"positive-solid union chooses maximum");
    }
    auto single=dict;single.boundaries.resize(1);single.boundaryById.erase(2);
    Embedded singleGeometry(Embedded::Mode::Explicit);singleGeometry.resolveSelection(single,2);singleGeometry.rebuildComposite(mesh,&geometry);
    check(!singleGeometry.usesComposite() && singleGeometry.levelSet()==&geometry.levelSets[0],"single level set is borrowed without copying");
    auto invalid=single;invalid.boundaries[0].levelSetIndex=99;invalid.boundaryById[1]=invalid.boundaries[0];
    singleGeometry.resolveSelection(invalid,2);
    check(!singleGeometry.ready() && !singleGeometry.levelSet() && !singleGeometry.metricsReady(),"invalid selection does not fall back over an explicit boundary");

    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={.5,.5};box.level=1;const int oldCells=mesh.cellIndexNum();
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()>oldCells,"embedded rebuild test changes topology");fill();
    explicitGeometry.rebuildComposite(mesh,&geometry);graph.build(mesh);explicitGeometry.rebuildMetrics(mesh,&geometry,graph,lineIds);
    check(explicitGeometry.alphaByCell().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"AMR replaces owned composite and cut-cell buffers");
    explicitGeometry.rebuildComposite(mesh,nullptr);explicitGeometry.rebuildMetrics(mesh,nullptr,graph,lineIds);
    check(!explicitGeometry.levelSet() && !explicitGeometry.usesComposite() && !explicitGeometry.metricsReady()
        && explicitGeometry.alphaByCell().empty(),"unavailable geometry clears the original derived views");
}
void checkChargeBalance() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<ParticleType> species;
    species.emplace_back("ion",1.,2.);species.emplace_back("electron",1.,-1.);species.emplace_back("neutral",1.,0.);
    const std::unordered_map<std::string,std::size_t> indices{{"ion",0},{"electron",1},{"neutral",2}};
    ParticleSpawnBuffer source;Particle p;p.w=0.;source.Push("ion",p);source.Push("missing",p);
    p.w=nan;source.Push("electron",p);p.w=3.;source.Push("electron",p);source.Push("neutral",p);
    const std::vector<std::array<double,4>> walls{{1.,2.,3.,4.},{2.,0.,0.,0.},{8.,0.,0.,0.}};
    const auto row=summarizeExplicitPicCurrentBalance(source,indices,species,{2.},walls,{1.,3.,8.},2.);
    check(row.sourceC==0. && row.domainC==18. && row.embeddedC==-1.,"signed source, domain and embedded totals");
    check(row.ionInjectionA==2. && row.electronInjectionA==2. && row.injectionMagnitudeA==4.
        && row.domainMagnitudeA==11. && row.embeddedMagnitudeA==2.5,"source fallback and four-side physical outflow");
    check(row.magnitudeResidualA==-9.5 && row.magnitudeRelative==-4.75
        && row.signedResidualA==-8.5 && row.signedRelative==-4.25,"current normalization uses ion injection only");
    for(double dt:{0.,-1.,nan}) {
        const auto zero=summarizeExplicitPicCurrentBalance(source,indices,species,{2.},walls,{1.,3.,8.},dt);
        check(zero.domainC==18. && zero.domainSignedA==0. && zero.magnitudeRelative==0.,"invalid dt retains charge but suppresses finite currents");
    }
    ParticleSpawnBuffer electronOnly;p.w=2.;electronOnly.Push("electron",p);
    const std::vector<std::array<double,4>> noOutflow(3);
    const auto floor=summarizeExplicitPicCurrentBalance(electronOnly,indices,species,{},noOutflow,{0.,0.,0.},1.);
    check(floor.magnitudeRelative==2./1e-30 && floor.signedRelative==-2./1e-30,"electron-only source retains denominator floor");
    p.w=0.;species[0].particle={p};p.w=nan;species[1].particle={p};
    check(sumPicMacroCharge(species,{2.})==3.,"macro charge sum retains fallback and no mass/finite filtering");
    std::vector<ParticleType> ordered;ordered.emplace_back("large",1.,1.);ordered.emplace_back("small",1.,1.);
    p.w=1e16;ordered[0].particle={p};p.w=1.;ordered[1].particle={p,p};
    check(sumPicMacroCharge(ordered,{})==1e16,"macro charge stays particle-wide, not per-species grouping");
    p.w=std::numeric_limits<double>::infinity();source.Push("ion",p);
    const auto infinite=summarizeExplicitPicCurrentBalance(source,indices,species,{2.},walls,{1.,3.,8.},0.);
    check(std::isinf(infinite.sourceC) && std::isnan(infinite.ionInjectionA),"nonfinite charge times zero dt retains legacy NaN");

    PicStepParticleBalance balance(10.);const auto& current=balance.current();
    check(std::isnan(current.macroEndC) && std::isnan(current.relative) && std::isnan(current.stageRelative),"unfinished charge diagnostics stay unavailable");
    balance.addStageMacro(-3.);balance.addStageMacro(-2.);balance.finishStage(-4.);
    check(current.stageResidualC==1. && current.stageRelative==-.2,"stage residual uses signed macro denominator");
    balance.addSecondary(-2.);balance.addSecondary(-1.);balance.includeSecondarySource();balance.addSource(7.);
    balance.addDomain(2.);balance.addEmbedded(-1.);balance.finish(14.);
    check(current.secondaryC==-3. && current.sourceC==4. && current.expectedEndC==13.
        && current.residualC==1. && current.relative==1./14.,"SEE subtotal, signed sinks, macro closure and norm");
    PicStepParticleBalance zero(0.);zero.finishStage(2.);zero.finish(0.);
    check(zero.current().stageResidualC==2. && zero.current().stageRelative==0. && zero.current().relative==0.,"zero stage macro and norm floor");
    PicStepParticleBalance transfers(0.);transfers.addSource(10.);transfers.addDomain(8.);transfers.addEmbedded(1.);transfers.finish(2.);
    check(transfers.current().relative==1./19.,"physical norm includes absolute source and sink totals");
    PicStepParticleBalance grouped(0.);grouped.addSecondary(1e16);grouped.addSecondary(-1e16);
    grouped.includeSecondarySource();grouped.addSource(1.);grouped.finish(1.);
    check(grouped.current().relative==0.,"SEE subtotal committed before ordinary sources");
    PicStepParticleBalance resumed(14.);resumed.finish(14.);
    check(resumed.current().sourceC==0. && resumed.current().relative==0.,"new step/restart has no stale source or sink history");
    resumed.addStageMacro(nan);resumed.finishStage(1.);resumed.addSource(nan);resumed.finish(14.);
    check(resumed.current().stageRelative==0. && std::isnan(resumed.current().relative),"stage NaN fallback and physical NaN propagation stay distinct");
}
void checkImplicitSpawnBalance() {
    using Balance=PicStepParticleBalance;
    using Spawn=Balance::Spawn;
    const double nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<ParticleType> species;
    species.emplace_back("ion",2.,2.);species.emplace_back("electron",4.,-1.);
    species.emplace_back("invalid_mass",-2.,0.);
    const std::vector<double> weights{2.};
    const std::unordered_map<std::string,std::size_t> index{{"ion",0},{"electron",1},{"invalid_mass",2}};
    Particle p;p.velocity={1.,2.};p.vz=2.;p.w=0.;species[0].particle.push_back(p);
    p.w=3.;species[0].particle.push_back(p);p.w=nan;species[1].particle.push_back(p);
    Balance balance(species,weights);
    const auto& secondary=balance.spawns(Spawn::Secondary);
    const auto& source=balance.spawns(Spawn::Source);
    const auto& collision=balance.spawns(Spawn::Collision);
    check(balance.current().macroBeginC==9. && balance.begin(0).macro==2
        && balance.begin(0).physical==5. && balance.begin(1).physical==1.,"step snapshots charge and fallback physical weights");
    species[0].particle.clear();
    check(balance.begin(0).physical==5. && PicSpeciesCount::sample(species[0],2.).macro==0,"begin snapshot survives live particle changes");

    ParticleSpawnBuffer half1,half2;
    p.w=1e16;half1.Push("ion",p);p.w=1.;half2.Push("ion",p);half2.Push("ion",p);
    half2.Push("unknown",p);
    balance.addSpawns(Spawn::Collision,half1,species,weights,index,2.);
    balance.addSpawns(Spawn::Collision,half2,species,weights,index,2.);
    check(collision.species[0].macro==3 && collision.species[0].physical==1e16,
        "MCC halves accumulate in particle order without regrouping small weights");
    check(collision.energyJ==0. && balance.injectionEnergyJ()==0. && balance.current().sourceC==0.,
        "collision spawns do not enter physical source charge or injection energy");

    ParticleSpawnBuffer see;
    p.w=-1.;see.Push("ion",p);p.w=3.;see.Push("electron",p);see.Push("unknown",p);
    balance.addSpawns(Spawn::Secondary,see,species,weights,index,2.);
    const auto energyNear=[](double actual,double expected) {
        check(std::abs(actual-expected)<=1e-13*std::max(1.,std::abs(expected)),"weighted 3D injected energy");
    };
    energyNear(secondary.energyJ,288.);energyNear(balance.injectionEnergyJ(),288.);
    check(secondary.species[0].macro==1 && secondary.species[0].physical==2.
        && secondary.species[1].macro==1 && secondary.species[1].physical==3.
        && balance.current().secondaryC==1. && balance.current().sourceC==0.,
        "retained SEE counts and charge precede the source subtotal commit");
    balance.includeSecondarySource();
    ParticleSpawnBuffer injected;
    p.w=nan;injected.Push("electron",p);p.w=5.;injected.Push("ion",p);
    p.w=1.;injected.Push("invalid_mass",p);injected.Push("unknown",p);
    balance.addSpawns(Spawn::Source,injected,species,weights,index,2.);
    energyNear(source.energyJ,252.);energyNear(balance.injectionEnergyJ(),540.);
    check(source.species[0].physical==5. && source.species[1].physical==1.
        && source.species[2].macro==1 && balance.current().sourceC==10.,
        "source uses retained compensated weights, unknown names skipped, invalid mass only suppresses energy");
    balance.finish(19.);check(balance.current().residualC==0.,"source and SEE share physical charge closure");

    Balance fresh(species,weights);ParticleSpawnBuffer dropped;
    fresh.addSpawns(Spawn::Secondary,dropped,species,weights,index,2.);fresh.includeSecondarySource();
    check(fresh.spawns(Spawn::Secondary).species[0].macro==0 && fresh.injectionEnergyJ()==0.
        && fresh.current().sourceC==0.,"new step and empty retained buffers carry no stale counts, charge or energy");
    fresh.addSpawns(Spawn::Source,injected,species,weights,index,0.);
    check(fresh.injectionEnergyJ()==0. && fresh.spawns(Spawn::Source).species[0].physical==5.,
        "invalid velocity scale preserves counts and the original kinetic helper guard");
}
void checkEnergyState() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<ParticleType> species;
    species.emplace_back("negative",2.,-1.);species.emplace_back("tiny_positive",2.,1e-31);
    species.emplace_back("neutral",2.,0.);species.emplace_back("negative_mass",-2.,-2.);
    species.emplace_back("nan_charge",2.,nan);
    Particle p;p.velocity={1.,2.};p.vz=2.;p.w=0.;species[0].particle.push_back(p);
    p.w=3.;species[0].particle.push_back(p);p.w=-1.;species[1].particle.push_back(p);
    p.w=nan;species[2].particle.push_back(p);p.w=1.;species[3].particle.push_back(p);species[4].particle.push_back(p);
    const auto explicitEnergy=summarizeExplicitPicKineticEnergy(species,{2.},2.);
    check(explicitEnergy.electronJ==144. && explicitEnergy.ionJ==36. && explicitEnergy.neutralJ==72.
        && explicitEnergy.total()==252.,"explicit species-first kinetic sum, 3D speed, charge sign and fallback weights");
    PicKineticEnergy implicitEnergy;
    implicitEnergy.addImplicit(p,2.,-1.,2.,2.);implicitEnergy.addImplicit(p,2.,-1.,3.,2.);
    implicitEnergy.addImplicit(p,2.,1e-31,1.,2.);implicitEnergy.addImplicit(p,2.,0.,1.,2.);
    implicitEnergy.addImplicit(p,-2.,-2.,1.,2.);implicitEnergy.addImplicit(p,2.,nan,1.,2.);
    check(implicitEnergy.electronJ==180. && implicitEnergy.ionJ==36. && implicitEnergy.neutralJ==72.,
        "implicit positive mass filter and near-neutral/NaN charge classification");
    implicitEnergy.addImplicit(p,2.,-1.,nan,2.);implicitEnergy.addImplicit(p,2.,-1.,0.,2.);
    p.vz=nan;implicitEnergy.addImplicit(p,2.,-1.,1.,2.);
    check(implicitEnergy.total()==288.,"implicit nonfinite energy and invalid weight skipped");
    species[0].particle.push_back(p);
    check(std::isnan(summarizeExplicitPicKineticEnergy(species,{2.},2.).electronJ),"explicit nonfinite energy is retained");
    std::vector<ParticleType> ordered;ordered.emplace_back("large",2.,-1.);ordered.emplace_back("small",2.,-1.);
    p.velocity={1.,0.};p.vz=0.;p.w=1e16;ordered[0].particle={p};p.w=1.;ordered[1].particle={p,p};
    check(summarizeExplicitPicKineticEnergy(ordered,{},1.).electronJ==1e16+2.,"explicit grouping must not become one particle-wide sum");

    PicPowerWindow window;
    window.add({10.,5.,2.,1.,7.,8.,9.});window.add({-2.,1.,1.,1.,3.,4.,5.});
    auto first=window.sample(2.,12.,8.);
    check(first.balance.totalJ==20. && first.balance.absorbedJ==5. && first.balance.fieldRate==0.
        && first.balance.residualJ==0. && first.potentialEnergyJ==-8. && first.potentialCurrentW==0.,"first explicit report has no derivative or residual");
    check(window.transfers().reactionLossJ==10. && window.transfers().seeEmissionJ==12.
        && window.transfers().ionImpactJ==14.,"all power transfer channels accumulate and remain until report finishes");
    window.finishReport();check(window.transfers().injectionJ==0. && window.transfers().ionImpactJ==0.,"report clears interval only");
    window.add({4.,2.,1.,2.,999.,888.,777.});auto second=window.sample(4.,20.,10.);
    check(second.balance.fieldRate==4. && second.balance.residualJ==7. && second.potentialCurrentW==-2.,
        "explicit actual report interval, history retained, and original closure excludes diagnostic-only channels");
    window.finishReport();auto zero=window.sample(4.,25.,12.);
    check(zero.balance.fieldRate==0. && zero.balance.residualJ==7.,"zero interval suppresses derivative but retains balance");
    window.sample(3.,30.,14.);auto afterBackwards=window.sample(5.,34.,16.);
    check(afterBackwards.balance.fieldRate==2. && afterBackwards.balance.residualJ==6.,"explicit history advances even for nonpositive interval");
    window.sample(6.,nan,1.);auto afterNan=window.sample(7.,40.,10.);
    check(std::isnan(afterNan.balance.fieldRate) && std::isnan(afterNan.balance.residualJ),"explicit valid flag does not filter NaN history");
    PicPowerWindow restarted;restarted.add({4.,2.,1.,2.});auto resume=restarted.sample(100.,40.,10.);
    check(resume.balance.fieldRate==0. && resume.balance.residualJ==0.,"resumed power window starts without history");

    PicStepEnergyHistory history;
    auto initial=history.sample(.5,12.,8.,8.,6.,3.,2.);
    check(initial.fieldRate==0. && initial.kineticRate==0. && initial.residualJ==0.,"first implicit step has no history");
    auto step=history.sample(.25,20.,10.,4.,2.,1.,2.);
    check(step.fieldRate==32. && step.kineticRate==8. && step.totalRate==40. && step.residualJ==7.,"implicit rates divide by physical dt and add separate derivatives");
    auto invalidDt=history.sample(0.,25.,12.,0.,0.,0.,0.);
    check(invalidDt.fieldRate==0. && invalidDt.kineticRate==0. && invalidDt.residualJ==7.,"implicit zero dt still updates balance/history");
    history.sample(-1.,30.,14.,0.,0.,0.,0.);step=history.sample(2.,34.,16.,0.,0.,0.,0.);
    check(step.fieldRate==2. && step.kineticRate==1.,"implicit nonpositive dt history retained");
    history.sample(1.,nan,1.,0.,0.,0.,0.);step=history.sample(1.,40.,10.,0.,0.,0.,0.);
    check(step.fieldRate==0. && step.kineticRate==0. && step.residualJ==0.,"implicit nonfinite prior pair invalidates both derivatives");
}
void checkDiagnosticWindow() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    PicDiagnosticWindow window(2);
    std::vector<long long> beginMacro{10,5};std::vector<double> beginPhysical{20.,50.};
    window.begin(.25,beginMacro,beginPhysical);beginMacro[0]=999;beginPhysical[0]=999.;
    window.begin(99.,{},{});
    check(window.startTimeCode()==.25 && window.totals().speciesMacroBegin[0]==10
        && window.totals().speciesPhysBegin[0]==20.,"window owns first snapshot and preserves active interval");
    const std::vector<std::array<std::uint64_t,4>> domainMacro{{1,2,0,1},{0,0,0,0}};
    const std::vector<std::array<double,4>> domainPhysical{{2.,4.,0.,2.},{0.,0.,0.,0.}};
    const std::vector<std::array<double,4>> domainEnergy{{1.,2.,0.,3.},{0.,0.,0.,0.}};
    const std::vector<double> embeddedPhysical{2.,0.},embeddedEnergy{4.,0.},managementPhysical{2.};
    const std::vector<long long> embeddedMacro{1,0},managementMacro{1};
    const std::unordered_map<std::string,std::size_t> indices{{"electron",0},{"ion",1}};
    const std::vector<double> weights{3.};
    ParticleSpawnBuffer source,mcc;Particle particle;
    particle.w=2.;source.Push("electron",particle);
    particle.w=nan;source.Push("electron",particle);
    particle.w=-1.;source.Push("ion",particle);source.Push("unknown",particle);
    particle.w=.5;mcc.Push("electron",particle);particle.w=4.;mcc.Push("ion",particle);
    std::vector<Particle> secondary(2);secondary[0].w=0.;secondary[1].w=2.;
    for(int step=0;step<2;++step) {
        window.addEnergy({1.,2.,3.,4.,5.,6.,7.});
        window.addBoundary({domainMacro,domainPhysical,domainEnergy,embeddedPhysical,embeddedMacro,
            embeddedEnergy,managementMacro,managementPhysical});
        window.addSpawns(PicDiagnosticWindow::Spawn::Source,source,indices,weights);
        window.addSpawns(PicDiagnosticWindow::Spawn::Collision,mcc,indices,weights);
        window.addSecondary(0,secondary,weights);window.addReactionEnergy({{"ionization",2.},{"excitation",1.}});
    }
    source.Clear();mcc.Clear();secondary.clear();
    const auto first=window.budget(0,{7,13.});
    check(first.source.macro==4 && first.source.physical==10. && first.mccSpawn.macro==2 && first.mccSpawn.physical==1.
        && first.seeSource.macro==4 && first.seeSource.physical==10.,"independent spawn channels and legacy weight fallbacks");
    check(first.outDomain.macro==8 && first.outDomain.physical==16. && first.outEmbedded.macro==2 && first.outEmbedded.physical==4.
        && first.particleManagement.macro==2 && first.particleManagement.physical==4.,"two-step domain/embedded/management accumulation");
    check(first.residual.macro==3 && first.residual.physical==8. && first.reactionNet.macro==1 && first.reactionNet.physical==4.,
        "species closure distinguishes particle management from reaction net");
    const auto second=window.budget(1,{8,60.});
    check(second.source.macro==2 && second.source.physical==2. && second.mccSpawn.physical==8.
        && second.particleManagement.macro==0 && second.particleManagement.physical==0.
        && second.residual.macro==1 && second.residual.physical==0.,"short weights and management arrays retain zero/default fallback");
    const auto& totals=window.totals();
    check(totals.energy.driveJ==2. && totals.energy.injectionJ==4. && totals.energy.domainAbsorptionJ==6.
        && totals.energy.embeddedAbsorptionJ==8. && totals.energy.reactionLossJ==10.
        && totals.energy.seeEmissionJ==12. && totals.energy.ionImpactJ==14.
        && totals.wallAbsEnergyJBySpecies[0][3]==6. && totals.ebAbsEnergyJBySpecies[0]==8.
        && totals.reactionEnergyJByChannel.at("ionization")==4.,"energy and channel sums use the same owned window");
    window.finishReport();check(window.totals().sourceMacroBySpecies[0]==4,"report close defers clearing");
    window.begin(7.,{7,8},{13.,60.});
    const auto reset=window.budget(0,{7,13.});
    check(window.startTimeCode()==7. && reset.residual.macro==0 && reset.residual.physical==0.
        && reset.source.macro==0 && reset.mccSpawn.macro==0 && reset.seeSource.macro==0 && reset.outDomain.macro==0
        && reset.outEmbedded.macro==0 && reset.particleManagement.macro==0 && window.totals().energy.ionImpactJ==0.
        && window.totals().reactionEnergyJByChannel.empty() && window.totals().wallAbsEnergyJBySpecies[0][3]==0.,
        "next window clears every cumulative channel and captures new begin state");
    PicDiagnosticWindow resumed(1);resumed.begin(3.,{7},{13.});
    check(resumed.budget(0,{7,13.}).residual.physical==0.,"resumed diagnostic starts with an empty window");
    resumed.addEnergy({1e16});resumed.addEnergy({-1e16});resumed.addEnergy({1.});
    check(resumed.totals().energy.driveJ==1.,"energy accumulation keeps step order");
    particle.w=std::numeric_limits<double>::infinity();source.Push("electron",particle);
    resumed.addSpawns(PicDiagnosticWindow::Spawn::Source,source,indices,weights);
    check(std::isinf(resumed.totals().sourcePhysBySpecies[0]),"positive infinite weights keep legacy propagation");
    const std::array<std::uint64_t,4> noMacro{};const std::array<double,4> noPhysical{};
    const auto ordered=summarizePicSpeciesBudget({{0,1e16},{0,-1e16},{0,1.},{0,0.},noMacro,noPhysical,{0,0.},{2,0.},{-1,0.}});
    check(ordered.residual.macro==-2 && ordered.reactionNet.macro==-1 && ordered.residual.physical==1.,
        "shared closure keeps signed residuals and left-to-right physical arithmetic");
    PicDiagnosticWindow empty(0);empty.begin(0.,{},{});empty.addBoundary({{},{},{},{},{},{},{},{}});empty.finishReport();
}
void checkBoundaryWork() {
    std::vector<ParticleType> species;species.emplace_back("positive",1.,2.);species.emplace_back("negative",1.,-1.);
    const std::vector<std::array<double,4>> outflow{{3.,2.,100.,100.},{1.,1.,100.,100.}};
    const double nan=std::numeric_limits<double>::quiet_NaN(),eps=PhysicalConstants::epsilon0_Fpm;
    for(bool trackSum:{false,true}) {
        PicBoundaryWork work(trackSum);
        work.sampleDirichlet({true,true,false,true});work.sampleVoltagesSI({3.,7.,11.,13.});
        for(std::size_t side=0;side<4;++side)work.sampleNormalFieldSI(side,10.*(side+1));
        std::array<double,4> source{1.,.25,0.,0.},secondary{.5,.75,0.,0.};
        const PicBoundaryWork::Injection injection{source,secondary};
        work.resetStep(.5);work.addDomain({2.,4.,5.,0.},species,outflow,trackSum?nullptr:&injection);
        check(work.current().powerW==(trackSum?-72.:-49.),"raw versus net outflow and side/area filtering");
        source[0]=999.;secondary[1]=999.; // The completed reduction must not retain input references.
        const std::unordered_map<int,double> embedded{{10,2.},{20,-1.},{30,0.},{40,nan},{50,4.}};
        std::vector<int> visited,expected;
        for(const auto& item:embedded)if(std::abs(item.second)>0.)expected.push_back(item.first);
        work.addEmbedded(embedded,[&](int id,double& voltage) {
            visited.push_back(id);voltage=id==10?5.:2.;return id!=50;
        });
        check(visited==expected,"embedded voltage resolution order and zero/NaN charge filtering");
        work.finishStep();
        check(work.current().powerW==(trackSum?-88.:-65.) && work.current().energyJ==(trackSum?-44.:-32.5)
            && work.current().potentialCurrentW==(trackSum?88.:0.),"embedded conduction, power sign and time integration");
        work.resetStep(1.);work.sampleDirichlet({false,false,true,true});
        work.sampleNormalFieldSI(2,37.);work.sampleNormalFieldSI(3,49.);
        const std::vector<std::array<double,4>> zero(2);
        work.addDomain({2.,4.,5.,2.},species,zero);work.finishStep();
        check(std::abs(work.current().powerW+619.*eps)<1e-22,
            "all sides acquire history, including previously non-Dirichlet/zero-area sides");
        check(work.current().energyJ==work.current().powerW,"unit-time work integration");
        if(trackSum)check(std::abs(work.current().potentialCurrentW-619.*eps)<1e-22,"raw phi-I sum retained separately");
    }
    PicBoundaryWork history(false);history.sampleDirichlet({true,false,false,false});history.sampleVoltagesSI({1.,0.,0.,0.});
    history.sampleNormalFieldSI(0,10.);history.resetStep(1.);history.addDomain({1.,0.,0.,0.},{},{});history.finishStep();
    check(history.current().energyJ==0.,"first/new or resumed work step has no displacement history");
    history.sampleNormalFieldSI(0,100.);history.resetStep(0.);
    history.addDomain({1.,0.,0.,0.},{},{});
    bool resolved=false;history.addEmbedded({{1,1.}},[&](int,double&) {resolved=true;return true;});history.finishStep();
    check(!resolved && history.current().energyJ==0. && history.current().powerW==0.,"nonpositive dt skips work and embedded resolution");
    history.sampleNormalFieldSI(0,14.);history.resetStep(2.);history.addDomain({1.,0.,0.,0.},{},{});history.finishStep();
    check(history.current().powerW==-2.*eps && history.current().energyJ==-4.*eps,"zero-dt work step must not advance history");
    history.resetStep(1.);history.addDomain({1.,0.,0.,0.},{},{});history.finishStep();
    check(history.current().powerW==0.,"reset clears totals and preserves field snapshots/history");
    history.sampleDirichlet({true,true,true,true});history.resetStep(1.);
    history.addDomain({0.,-1.,nan,0.},species,{});history.finishStep();
    check(history.current().powerW==0.,"invalid areas do not read outflow entries");
    history.resetStep(1.);history.sampleNormalFieldSI(0,nan);history.addDomain({1.,0.,0.,0.},{},{});history.finishStep();
    check(std::isnan(history.current().energyJ),"nonfinite normal fields propagate into work");
}
void checkBoundaryFieldSampling(const QuadtreeMesh& mesh) {
    UnitsScales scales;scales.L0_m=2.;scales.Phi0_V=10.;
    const auto units=UnitsContext::MakeSI(scales);
    PicDomainBC bc;const std::array<PicBoundarySide*,4> sides{&bc.left,&bc.right,&bc.bottom,&bc.top};
    for(int i=0;i<4;++i){auto& side=*sides[i];side.value=100.+i;side.timeDirichlet.enabled=true;
        side.timeDirichlet.dc=3.+i;side.timeDirichlet.ac=2.;side.timeDirichlet.freq_hz=1.;}
    bc.right.type=PicBCType::Neumann;bc.top.timeDirichlet.enabled=false;
    FaceSegmentGraph graph;graph.build(mesh);std::vector<double> normal(graph.segments().size());
    for(std::size_t i=0;i<normal.size();++i)normal[i]=2.+i;
    const std::array<const char*,5> selected{"left","right","bottom","top","invalid"};
    for(int index=0;index<5;++index)for(bool useUnits:{false,true})for(double fallback:{-1.,0.,7.}) {
        PicElectrodeCurrent actual(selected[index],bc),expected(selected[index],bc);
        for(int step=0;step<2;++step) {
            auto fields=normal;for(auto& value:fields)value+=step;
            actual.sampleNormalFieldCode(graph,fields);expected.sampleNormalFieldCode(graph,fields);
            const double time=.125*step;const std::optional<double> scale=fallback<0.?std::nullopt:std::optional<double>(fallback);
            actual.sampleFieldSI(bc,time,useUnits?&units:nullptr,scale);
            double voltage=0.;if(index<4)voltage=PicFieldBoundaryState::voltage(*sides[index],time);
            expected.sampleFieldSI(useUnits?units.ElectricFieldCodeToSI(expected.normalFieldCode()):
                    (scale?*scale*expected.normalFieldCode():expected.normalFieldCode()),
                useUnits?units.PotentialCodeToSI(voltage):voltage);
            actual.finishStep(.5);expected.finishStep(.5);const auto& a=actual.current();const auto& b=expected.current();
            check(a.ready==b.ready && a.voltageV==b.voltageV && (index==4 ||
                (a.displacement==b.displacement && a.conduction==b.conduction && a.total==b.total)),
                "waveform field/voltage sampling preserves selection, time, units, zero/fallback scale and history");
        }
    }
    PicEmbeddedBC embedded;embedded.enabled=true;embedded.type=PicBCType::Dirichlet;embedded.g=19.;
    embedded.timeDirichlet.enabled=true;embedded.timeDirichlet.dc=3.;embedded.timeDirichlet.ac=2.;embedded.timeDirichlet.freq_hz=1.;
    PicEmbeddedBCMap byId;byId[7]=embedded;byId[7].g=23.;byId[7].timeDirichlet.dc=8.;
    byId[8]=embedded;byId[8].type=PicBCType::Neumann;byId[9]=embedded;byId[9].enabled=false;
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const std::unordered_map<int,double> charge{{7,2.},{8,4.},{9,-1.},{10,3.},{99,0.},{101,nan}};
    std::vector<ParticleType> species;species.emplace_back("ion",1.,2.);
    const std::vector<std::array<double,4>> outflow{{2.,3.,4.,5.}};
    const std::array<double,4> areas{2.,3.,4.,5.},source{.1,.2,.3,.4},secondary{.4,.3,.2,.1};
    const PicBoundaryWork::Injection injection{source,secondary};
    for(bool stored:{false,true}) {
        PicFieldBoundaryState state(bc,embedded,byId);PicBoundaryWork actual(stored),expected(stored);
        for(int step=0;step<3;++step) {
            const double time=.125*step,dt=step==2?0.:.5;
            if(step==1)state.update(time);
            const auto& domain=stored?state.domain():bc;
            const std::array<const PicBoundarySide*,4> selectedSides{&domain.left,&domain.right,&domain.bottom,&domain.top};
            auto fields=normal;for(auto& value:fields)value+=step;
            actual.sampleNormalFieldsSI(graph,fields,units);
            const auto averages=averagePicBoundaryNormalFields(graph.segments(),fields);
            for(std::size_t s=0;s<4;++s)expected.sampleNormalFieldSI(s,units.ElectricFieldCodeToSI(averages[s]));
            actual.sampleVoltagesSI(domain,units,stored?std::nullopt:std::optional<double>(time));
            std::array<bool,4> dirichlet{};std::array<double,4> voltages{};
            for(std::size_t s=0;s<4;++s){dirichlet[s]=selectedSides[s]->type==PicBCType::Dirichlet;
                voltages[s]=dirichlet[s]?units.PotentialCodeToSI(stored?selectedSides[s]->value:
                    PicFieldBoundaryState::voltage(*selectedSides[s],time)):0.;}
            expected.sampleDirichlet(dirichlet);expected.sampleVoltagesSI(voltages);
            actual.resetStep(dt);expected.resetStep(dt);
            actual.addDomain(areas,species,outflow,stored?nullptr:&injection);
            expected.addDomain(areas,species,outflow,stored?nullptr:&injection);
            if(stored)actual.addEmbedded(charge,state,units);
            else actual.addEmbedded(charge,embedded,byId,units,time);
            expected.addEmbedded(charge,[&](int id,double& voltage){double code=0.;
                if(!(stored?state.embeddedVoltage(id,code):PicFieldBoundaryState::embeddedVoltage(embedded,byId,id,time,code)))return false;
                voltage=units.PotentialCodeToSI(code);return true;});
            actual.finishStep();expected.finishStep();
            check(actual.current().energyJ==expected.current().energyJ && actual.current().powerW==expected.current().powerW
                && actual.current().potentialCurrentW==expected.current().potentialCurrentW,
                "work sampling preserves raw/evaluated voltages, boundaries, units, injection, map fallback and displacement history");
            check(step==2 || actual.current().powerW!=0.,"boundary sampling fixture has nonzero work");
        }
    }
    check(bc.left.value==100. && embedded.g==19. && byId[7].g==23.,"sampling does not overwrite prepared boundary inputs");
}
void checkElectrodeCurrent(const QuadtreeMesh& mesh) {
    PicDomainBC bc;
    check(PicElectrodeCurrent("auto",bc).boundary()==Boundary::west,"first Dirichlet electrode");
    bc.right.timeDirichlet.enabled=true;bc.bottom.timeDirichlet.enabled=true;
    PicElectrodeCurrent selected("AUTO",bc);
    check(selected.sideIndex()==1,"driven electrode priority before static left");
    bc.right.type=PicBCType::Neumann;
    check(selected.sideIndex()==1 && PicElectrodeCurrent("",bc).sideIndex()==2,"selection owns its identity and ignores driven non-Dirichlet");
    const std::array<const char*,8> names{"LEFT","west","Right","east","Bottom","south","Top","north"};
    for(int i=0;i<8;++i)check(PicElectrodeCurrent(names[i],bc).sideIndex()==i/2,"electrode side aliases/order");
    bc.left.type=bc.right.type=bc.bottom.type=bc.top.type=PicBCType::Neumann;
    PicElectrodeCurrent invalid("auto",bc);invalid.finishStep(.5);
    check(!invalid.current().ready && std::isnan(invalid.current().conduction) && std::isnan(invalid.current().total),
        "unresolved electrode current is unavailable");
    check(PicElectrodeCurrent(" left ",bc).sideIndex()==-1,"legacy side tokens are not whitespace-trimmed");
    PicElectrodeCurrent current("west",bc);FaceSegmentGraph graph;graph.build(mesh);
    current.sampleNormalFieldCode(graph,std::vector<double>(graph.segments().size(),12.));
    check(current.normalFieldCode()==12.,"waveform owns sampled code-unit field");
    std::vector<ParticleType> species;species.emplace_back("positive",1.,2.);species.emplace_back("negative",1.,-1.);
    const std::vector<std::array<double,4>> outflow{{3.,100.,100.,100.},{1.,100.,100.,100.}};
    current.sampleConductionSI(species,outflow,2.,.5);current.sampleFieldSI(10.,7.);current.finishStep(.5);
    check(current.current().ready && current.current().voltageV==7. && current.current().conduction==5.
        && current.current().displacement==0. && current.current().total==5.,"signed outflow/area/time and first-step displacement");
    const double eps=PhysicalConstants::epsilon0_Fpm;
    current.sampleFieldSI(14.,8.);current.finishStep(.5);
    check(current.current().displacement==eps*8.,"consecutive sampled field difference");
    current.sampleFieldSI(22.,9.);current.finishStep(2.);
    check(current.current().displacement==eps*4.,"field history advances on steps without CSV writes");
    current.sampleConductionSI(species,{},0.,1.);current.sampleFieldSI(30.,9.);current.finishStep(0.);
    check(current.current().conduction==0. && current.current().displacement==0.,"invalid area/dt guards");
    current.sampleFieldSI(34.,9.);current.finishStep(1.);
    check(current.current().displacement==eps*4.,"zero-dt step still updates legacy history");
    current.sampleFieldSI(std::numeric_limits<double>::quiet_NaN(),0.);current.finishStep(1.);
    check(std::isnan(current.current().displacement),"undefined field must propagate to current");
    current.sampleFieldSI(1.,0.);current.finishStep(1.);
    check(std::isnan(current.current().displacement),"undefined previous field must not be silently reset");
    current.sampleFieldSI(2.,0.);current.finishStep(1.);
    check(current.current().displacement==eps,"finite history recovers after a full finite interval");
    PicElectrodeCurrent resumed("west",bc);resumed.sampleFieldSI(100.,0.);resumed.finishStep(1.);
    check(resumed.current().displacement==0.,"new/resumed diagnostic has no previous field");
    PicElectrodeCurrent noCsv("west",bc);noCsv.sampleConductionSI(species,{},0.,1.);noCsv.finishStep(1.);
    check(noCsv.current().ready && noCsv.current().total==0.,"diagnostic availability is independent of disabled CSV");
}
void checkBoundaryDiagnostics(const QuadtreeMesh& refined) {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<FaceSegment> segments(8);
    const std::array<Boundary,8> sides{Boundary::west,Boundary::west,Boundary::east,Boundary::north,
        Boundary::non,Boundary::west,Boundary::south,Boundary::north};
    const std::array<double,8> areas{1.,3.,2.,4.,100.,-1.,0.,nan};
    for(std::size_t i=0;i<segments.size();++i) {segments[i].boundary=sides[i];segments[i].Af=areas[i];}
    std::vector<double> fields{2.,6.,-4.,9.,1000.,1000.,1000.,1000.};
    check(averagePicBoundaryNormalFields(segments,fields)==std::array<double,4>{5.,-4.,0.,9.},
        "boundary averages: geometric area, side order and invalid areas");
    check(averagePicBoundaryNormalField(segments,fields,Boundary::west)==5.,"selected electrode average");
    check(averagePicBoundaryNormalField(segments,{2.},Boundary::west)==2.
        && averagePicBoundaryNormalField(segments,{},Boundary::west)==0.,"short/empty segment field array");
    fields[0]=nan;
    check(std::isnan(averagePicBoundaryNormalField(segments,fields,Boundary::west))
        && std::isnan(averagePicBoundaryNormalFields(segments,fields)[0]),"diagnostics must retain undefined fields");

    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);CellField<double> phi(mesh,0.);
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();phi.value[ci]=p.x+2.*p.y;
    }
    using Policy=PicPeriodicMatePolicy;
    for(auto policy:{Policy::AnyMate,Policy::CleanMate}) {
        const auto jumps=measurePicPeriodicJumps(mesh,graph,phi,true,true,nullptr,policy,true);
        check(jumps.maxJumpX==.5 && jumps.maxJumpY==1.,"analytic periodic cell-center jumps");
        check(jumps.debugX.find("jump=0.5")!=std::string::npos && jumps.debugY.find("jump=1")!=std::string::npos,
            "periodic debug owns both measured maxima");
        const auto disabled=measurePicPeriodicJumps(mesh,graph,phi,false,false,nullptr,policy,true);
        check(disabled.maxJumpX==0. && disabled.maxJumpY==0. && disabled.debugX.empty() && disabled.debugY.empty(),
            "disabled periodic axes");
        Solvers::Poisson::CutCellMetrics wrongSize;
        const auto plain=measurePicPeriodicJumps(mesh,graph,phi,true,false,&wrongSize,policy);
        check(plain.maxJumpX==.5 && plain.maxJumpY==0. && plain.debugX.empty(),"mismatched cut metrics and disabled debug");
    }
    FaceSegmentGraph coarseFine;coarseFine.build(refined);CellField<double> asymmetric(refined,0.);
    for(int ci=0;ci<refined.cellIndexNum();++ci)
        asymmetric.value[ci]=refined.cells(refined.cellIndexToLabel[ci]).centerCoord().y;
    const auto any=measurePicPeriodicJumps(refined,coarseFine,asymmetric,false,true,nullptr,Policy::AnyMate,true);
    const auto clean=measurePicPeriodicJumps(refined,coarseFine,asymmetric,false,true,nullptr,Policy::CleanMate,true);
    check(any.maxJumpY==.625 && !any.debugY.empty() && clean.maxJumpY==0. && clean.debugY.empty(),
        "retain explicit/implicit mate policy on actual asymmetric refinement");

    DictionaryGeometry geometryInput;Line wall;wall.points[0]={0.,.625};wall.points[1]={1.,.625};wall.setType(LineType::straight);
    geometryInput.lines={wall};geometryInput.levelSetName={"wall"};geometryInput.innerPoints={{.5,.25}};
    geometryInput.linesOnInterface={{{0,false}}};geometryInput.moving={false};
    ImmersedGeometry<QuadtreeMesh> geometry(geometryInput,mesh);
    VertField<double> psi(mesh,0.),dx(mesh,0.),dy(mesh,1.);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)psi.value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().y-.625;
    Solvers::Poisson::CutCellMetrics cut;Solvers::Poisson::CutCellMetrics::BuildOptions options;
    options.levelSetI=-1;options.psiOverride=&psi;options.psiDxOverride=&dx;options.psiDyOverride=&dy;
    cut.build(mesh,geometry,graph,options);
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();phi.value[ci]=p.x*(p.y<.5?1.:10.);
    }
    const auto unmasked=measurePicPeriodicJumps(mesh,graph,phi,true,true,nullptr,Policy::CleanMate);
    const auto masked=measurePicPeriodicJumps(mesh,graph,phi,true,true,&cut,Policy::CleanMate,true);
    check(unmasked.maxJumpX==5. && masked.maxJumpX==.5 && masked.maxJumpY==0.,
        "periodic cut filtering must exclude the cut row and all vertical pairs");
    check(masked.debugX.find("beta=1 alphaK=1 alphaL=1 cutK=0 cutL=0")!=std::string::npos,
        "periodic diagnostics retain open-cell metric details");
}
void checkOutputCoupling(const QuadtreeMesh& mesh) {
    ParticleType particles("sample",1.,-1.);particles.particle.resize(2);
    std::vector<Vector2D> positions;std::vector<int> hosts;
    for(int i=0;i<2;++i) {
        auto& p=particles.particle[i];p.position={i==0?.03:.97,.6};p.inCellL=mesh.findCellLabel(p.position);
        p.w=2.+i;positions.push_back(p.position);hosts.push_back(p.inCellL);
    }
    const std::vector<int> ignoredHosts(2,-1);std::vector<double> alpha(mesh.cellIndexNum(),1.);alpha[0]=.5;
    CouplingWeightsCache::BuildOptions options;options.periodicX=true;options.periodicY=true;
    options.levelSetI=-1;options.enableCompositeBand=true;options.enableMomentExactCut=false;
    options.alphaByCellIndex=&alpha;options.cutByAlpha=true;options.hostCellLabels=&ignoredHosts;
    options.dumpWeightsOnce=false;
    auto originalOptions=options;originalOptions.hostCellLabels=&hosts;
    const std::vector<double> halfWidths{.4,.2};CouplingWeightsCache expected,actual;
    expected.build(mesh,positions,halfWidths,originalOptions);
    PicCurrentWeights::buildDiagnostic(mesh,particles,halfWidths,options,actual);
    check(actual.offsets()==expected.offsets() && actual.cellIndices()==expected.cellIndices()
        && actual.weights()==expected.weights(),"diagnostic coupling must retain particle order, hosts and supplied periodic/cut options");
    check(options.hostCellLabels==&ignoredHosts && particles.particle[0].position.x==.03
        && particles.particle[1].inCellL==hosts[1],"diagnostic coupling modified caller options or particles");
    const auto saved=actual.weights();particles.particle[0].position={.8,.8};
    check(actual.weights()==saved,"diagnostic cache must own generated weights");
    ParticleType empty("empty",1.,0.);PicCurrentWeights::buildDiagnostic(mesh,empty,halfWidths,options,actual);
    check(actual.offsets()==std::vector<int>({0}) && actual.cellIndices().empty() && actual.weights().empty(),
        "empty diagnostic species must clear a reused cache");
}
void checkFieldFrame() {
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    CellField<double> phi(mesh,1.),raw(mesh,0.);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)
        raw.value[ci]=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x<.5?-3.:1.;
    std::vector<unsigned char> valid(mesh.cellIndexNum(),1);
    PicFieldSelection explicitSelection(PicFieldMode::Explicit,{"n","mean_energy","te","w_mean","w_sum","w_max","np"});
    PicFieldSelection implicitSelection(PicFieldMode::Implicit,{"te","n","mean_energy"});
    PicFieldFrame explicitFrame(mesh,explicitSelection,valid,phi);
    PicFieldFrame implicitFrame(mesh,implicitSelection,valid,phi);
    explicitFrame.addScalar("clamp",raw,[](double v) {return v;},true);
    implicitFrame.addScalar("clamp",raw,[](double v) {return v;},true);
    valid[0]=0;raw.value[0]=99.;
    check(explicitFrame.validCells[0]==1 && implicitFrame.validCells[0]==1,"frame must own its mask");
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[vi]).coord();
        if(p.x==.5 && p.y==.5)
            check(explicitFrame.vertices[0].value[vi]==0. && implicitFrame.vertices[0].value[vi]==.5,
                "frame must preserve clamp/interpolation order and own prepared values");
    }
    valid[0]=1;
    auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(UnitsScales{}));
    NormalizationConfig config;config.slabThickness_m=1.;
    const auto norm=NormalizationContext::MakeFromUnits(units,config);
    std::vector<ParticleType> species;
    species.emplace_back("negative_ion",2.*PhysicalConstants::eV_J,-1.);
    species.emplace_back("hot_ELECTRON",2.*PhysicalConstants::eV_J,1.);
    species.emplace_back("empty",1.,0.);
    for(int si=0;si<2;++si) {
        species[si].particle.resize(1);auto& p=species[si].particle[0];
        p.w=si==0?2.:4.;p.velocity={si==0?1.:3.,0.};p.inCellL=mesh.cellIndexToLabel[0];
        p.position=mesh.cells(p.inCellL).centerCoord();
    }
    std::vector<int> calls(3,0);
    auto build=[&](std::size_t si,CouplingWeightsCache& cache) {
        ++calls[si];cache.buildFromCSR(1,{0,1},{1},{1.});
    };
    explicitFrame.addParticles(species,{},norm,nullptr,build);
    check(calls==std::vector<int>({1,2,0}),"explicit energy must rebuild the named electron cache after all densities");
    const std::vector<std::string> expected={"clamp","n_negative_ion_SI_1pm3","n_hot_electron_SI_1pm3",
        "n_empty_SI_1pm3","electron_mean_energy_eV","Te_e_eV","w_mean","w_sum","w_max","n_macro"};
    explicitFrame.finish(true,true);
    check(explicitFrame.names==expected && explicitFrame.cells.empty(),"explicit fixed order and vertex-only frame");
    calls.assign(3,0);
    PicFieldFrame implicitCells(mesh,implicitSelection,valid,phi,true);
    implicitCells.addParticles(species,{},norm,nullptr,build);
    check(calls==std::vector<int>({1,1,0}),"implicit density and energy must share one cache per species");
    implicitCells.finish(true,true);
    check(implicitCells.names==std::vector<std::string>({"Te_e_eV","n_negative_ion_SI_1pm3",
        "n_hot_electron_SI_1pm3","n_empty_SI_1pm3","electron_mean_energy_eV"}),"implicit sorted particle columns");
    check(implicitCells.cellConstant() && implicitCells.vertices.empty(),"implicit cell frame layout");
    check(std::abs(implicitCells.cells[0].value[1]-6.)<1e-12 && implicitCells.cells[0].value[0]==0.,
        "named electron takes precedence over the first negative species; CSR energy location");
    check(implicitCells.cells[2].value[1]==16. && implicitCells.cells[3].value[1]==0.,"density amount/volume and empty species");
    calls.assign(3,0);
    std::vector<CellField<double>> amounts;amounts.emplace_back(mesh,7.);amounts.emplace_back();
    PicFieldFrame fast(mesh,implicitSelection,valid,phi,true);
    fast.addParticles(species,{},norm,nullptr,build,true,&amounts);fast.finish(false,false);
    check(calls==std::vector<int>({0,1,0}),"fast density uses supplied amount, wrong-sized amount uses CSR once");
    check(fast.cells[1].value[0]==28. && fast.cells[2].value[1]==16.,"fast density and CSR fallback values");
    check(std::abs(fast.cells[0].value[0]-6.)<1e-12 && fast.cells[0].value[1]==0.,
        "fast energy retains host-cell sampling after a density CSR fallback");
    amounts[0].value[0]=100.;species[1].particle[0].w=100.;
    check(fast.cells[1].value[0]==28. && implicitCells.cells[2].value[1]==16.,"frame borrows no deposited/particle value storage");
    calls.assign(3,0);
    PicFieldSelection unknown(PicFieldMode::Implicit,{"unknown"});
    PicFieldFrame empty(mesh,unknown,valid,phi);
    empty.addParticles(species,{},norm,nullptr,build);empty.finish(false,false);
    check(empty.names.empty() && calls==std::vector<int>({0,0,0}),"unselected particle fields must not build caches");
}
void checkGeometryFields(const QuadtreeMesh& mesh) {
    const auto near=[](double actual,double expected) {
        if (!std::isfinite(actual) || std::abs(actual-expected)>=1e-12*std::max(1e-10,std::abs(expected)))
            std::cerr << "geometry/electric actual=" << actual << " expected=" << expected << '\n';
        check(std::isfinite(actual) && std::abs(actual-expected)<1e-12*std::max(1e-10,std::abs(expected)),
            "geometry/electric output differs from analytic value");
    };
    UnitsScales scales;scales.L0_m=2.;scales.Phi0_V=10.;
    auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(scales));
    const auto norm=NormalizationContext::MakeFromUnits(units,NormalizationConfig{});
    CellField<double> phi(mesh,2.),rho(mesh,3.);
    CellField<Vector2D> electric(mesh,Vector2D{4.,-2.});
    VertField<double> composite(mesh,-.5);
    std::vector<LevelSet> levelSets;levelSets.emplace_back(mesh);levelSets.emplace_back(mesh);
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.5;
        levelSets[1].value[vi]=-.25;
    }
    std::vector<unsigned char> valid(mesh.cellIndexNum(),1);
    PicFieldSelection selection(PicFieldMode::Implicit,{"phi","rho","E","psi_all","psi_composite","muCell","muFace"});
    PicFieldFrame frame(mesh,selection,valid,phi,true);
    frame.addElectrostatic(phi,rho,*units,norm);frame.addElectric(electric,*units);
    frame.addLevelSets(*units,&composite,&levelSets,{-1,1,0,1,99},{"Solid+",""});
    check(frame.names==std::vector<std::string>({"phi_SI_V","rho_SI_Cpm3","Ex_SI_Vpm","Ey_SI_Vpm",
        "psi_unnamed_SI_m","psi_solid_SI_m","psi_unnamed_SI_m","psi_composite_SI_m"}),
        "level-set input order, invalid-index skip, empty-name spelling and duplicate retention");
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        near(frame.cells[0].value[ci],20.);near(frame.cells[1].value[ci],7.5*PhysicalConstants::epsilon0_Fpm);
        near(frame.cells[2].value[ci],20.);near(frame.cells[3].value[ci],-10.);
        near(frame.cells[4].value[ci],-.5);near(frame.cells[6].value[ci],-.5);near(frame.cells[7].value[ci],-1.);
        near(frame.cells[5].value[ci],2.*(mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x-.5));
    }
    PicFieldFrame missing(mesh,selection,valid,phi,true);
    missing.addLevelSets(*units,nullptr,&levelSets,{1});
    check(missing.names==std::vector<std::string>({"psi_ls_2_SI_m","psi_composite_SI_m"}),
        "missing level-set name uses one-based index; missing composite still outputs zero");
    for(int ci=0;ci<mesh.cellIndexNum();++ci)near(missing.cells[1].value[ci],0.);
    PicFieldSelection explicitSelection(PicFieldMode::Explicit,{"phi","rho","E","psi"});
    PicFieldFrame explicitFrame(mesh,explicitSelection,valid,phi);
    explicitFrame.addElectrostatic(phi,rho,*units,norm);explicitFrame.addLevelSets(*units,nullptr);
    explicitFrame.addElectric(electric,*units);
    check(explicitFrame.names==std::vector<std::string>({"phi_SI_V","rho_SI_Cpm3","psi_SI_m","Ex_SI_Vpm","Ey_SI_Vpm"}),
        "explicit electric fields stay after scalar/geometry fields");
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        near(explicitFrame.vertices[0].value[vi],20.);near(explicitFrame.vertices[2].value[vi],0.);
        near(explicitFrame.vertices[3].value[vi],20.);near(explicitFrame.vertices[4].value[vi],-10.);
    }
    FaceSegmentGraph graph;graph.build(mesh);
    std::size_t chosen=0;
    for(std::size_t sid=0;sid<graph.segments().size();++sid) {
        const auto& face=graph.segments()[sid];
        if(face.K>=0) {
            const auto range=graph.segmentRangeForCellIndex(mesh.cells(face.K).index());
            if(range.second-range.first>4) {chosen=sid;break;}
        }
    }
    std::vector<double> segments(chosen+1,1.);segments[chosen]=9.;
    frame.addSusceptibility(graph,{7.},{5.},segments);
    const auto& changed=graph.segments()[chosen];bool coarse=false;
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const auto range=graph.segmentRangeForCellIndex(ci);
        // Each interior face has an owner-oriented entry from both cells;
        // the legacy proxy counts both K and L visits, boundaries only once.
        int count=0;
        for(std::size_t sid=range.first;sid<range.second;++sid)
            count+=graph.segments()[sid].L>=0?2:1;
        const int label=mesh.cellIndexToLabel[ci];const bool touched=label==changed.K || label==changed.L;
        if(touched && range.second-range.first>4)coarse=true;
        near(frame.cells[8].value[ci],ci==0?7.:0.);
        near(frame.cells[9].value[ci],touched?8./count*PhysicalConstants::epsilon0_Fpm:0.);
    }
    check(coarse,"susceptibility test must exercise an incident coarse/fine cell");
    FaceSegmentGraph noFaces;noFaces.build(mesh,[](const QuadtreeMesh&,const Cell&) {return false;});
    PicFieldFrame fallback(mesh,selection,valid,phi,true);fallback.addSusceptibility(noFaces,{}, {5.},{});
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        near(fallback.cells[0].value[ci],0.);
        near(fallback.cells[1].value[ci],ci==0?4.*PhysicalConstants::epsilon0_Fpm:0.);
    }
    rho.value[0]=99.;electric.value[0]={99.,99.};composite.value[0]=99.;
    near(frame.cells[1].value[0],7.5*PhysicalConstants::epsilon0_Fpm);
    near(frame.cells[2].value[0],20.);near(frame.cells[7].value[0],-1.);
}
void checkFieldSelection(const QuadtreeMesh& mesh) {
    PicFieldSelection explicitDefault(PicFieldMode::Explicit,{});
    const auto& defaults=explicitDefault.fields();
    check(explicitDefault.configured() && defaults.phi && defaults.rho && defaults.psi && defaults.ex && defaults.ey
        && defaults.n && defaults.te && defaults.macroCount && !defaults.meanEnergy && !defaults.psiAll
        && !defaults.muCell && !defaults.weightMean,"explicit empty-list defaults");
    PicFieldSelection implicitEmpty(PicFieldMode::Implicit,{});
    check(!implicitEmpty.configured() && !implicitEmpty.fields().n && !implicitEmpty.fields().phi,
        "implicit empty selection must stay empty");
    PicFieldSelection explicitUnknown(PicFieldMode::Explicit,{"unrecognized"});
    PicFieldSelection implicitUnknown(PicFieldMode::Implicit,{"unrecognized"});
    check(explicitUnknown.fields().n && !explicitUnknown.fields().te && !explicitUnknown.fields().phi,
        "explicit densities are mandatory even for unknown-only selection");
    check(implicitUnknown.configured() && !implicitUnknown.fields().n && !implicitUnknown.fields().phi,
        "implicit configured unknown token must not enable defaults");
    std::vector<std::string> requested={"TEMPERATURE","n","PHI","electricfield","RHO","levelset_all",
        "psi_composite","MUCell","MUFace","wmean","wphys","wmax","np","mean_energy"};
    PicFieldSelection implicitAliases(PicFieldMode::Implicit,requested),explicitAliases(PicFieldMode::Explicit,requested);
    requested.clear();
    const auto& i=implicitAliases.fields();const auto& e=explicitAliases.fields();
    check(i.te && i.n && i.phi && i.ex && i.ey && i.rho && i.psiAll && i.psiComposite && !i.psi
        && i.muCell && i.muFace && i.meanEnergy && !i.weightMean && !i.macroCount,"implicit supported alias selection");
    check(e.te && e.n && e.phi && e.ex && e.ey && e.rho && !e.psiAll && !e.psiComposite && !e.muCell
        && e.weightMean && e.weightSum && e.weightMax && e.macroCount && e.meanEnergy,"explicit aliases and unsupported implicit fields");
    PicFieldSelection psiAlias(PicFieldMode::Implicit,{"LEVEL-SET"});
    check(psiAlias.fields().psi && psiAlias.fields().psiComposite && !psiAlias.fields().psiAll,"psi alias implies composite field only");
    check(explicitAliases.key(" __AR+*-__ ")=="arplusstarminus" && implicitAliases.key(" __AR+*-__ ")=="ar",
        "preserve distinct explicit/implicit symbol spelling");
    check(explicitAliases.key("  A__ B  ")=="a_b" && implicitAliases.key("!")=="unnamed",
        "output key underscore collapse and empty fallback");
    auto explicitKeys=explicitAliases.keys(4);auto implicitKeys=implicitAliases.keys(4);
    check(explicitKeys.next("Ar+")=="arplus" && explicitKeys.next("ARPLUS")=="arplus_2"
        && explicitKeys.next("Ar-")=="arminus","explicit repeated base naming");
    check(implicitKeys.next("Ar+")=="ar" && implicitKeys.next("Ar-")=="ar_2"
        && implicitKeys.next("ar_2")=="ar_2","retain per-base suffix behavior including an independently suffixed base");
    auto nextFrame=implicitAliases.keys(1);check(nextFrame.next("Ar-")=="ar","key counters must not leak across output frames");

    std::vector<std::string> names={"phi_SI_V","Ey_SI_Vpm","n_e_SI_1pm3","Te_e_eV","Ex_SI_Vpm","n_i_SI_1pm3",
        "rho_SI_Cpm3","psi_gate_SI_m","psi_composite_SI_m","muCell_SI_Fpm","muFaceProxy_SI_Fpm","opaque"};
    const auto originalNames=names;
    std::vector<CellField<double>> cells;
    for(std::size_t index=0;index<names.size();++index)cells.emplace_back(mesh,static_cast<double>(index));
    explicitAliases.order(names,cells);
    check(names==originalNames,"explicit output order is fixed, not user token order");
    implicitAliases.order(names,cells);
    const std::vector<int> expected={3,2,5,0,1,4,6,7,8,9,10,11};
    for(std::size_t index=0;index<expected.size();++index) {
        check(names[index]==originalNames[expected[index]],"implicit field order or stable species/component tie");
        for(int ci=0;ci<mesh.cellIndexNum();++ci)check(cells[index].value[ci]==expected[index],"field/name identity lost during ordering");
    }
    PicFieldSelection duplicates(PicFieldMode::Implicit,{"n","TE_E","n","PHI","Te"});
    std::vector<std::string> vertexNames={"opaque","phi_SI_V","Te_e_eV","n_e_SI_1pm3","rho_SI_Cpm3"};
    const auto oldNames=vertexNames;std::vector<VertField<double>> vertices;
    for(std::size_t index=0;index<vertexNames.size();++index)vertices.emplace_back(mesh,static_cast<double>(index));
    duplicates.order(vertexNames,vertices);
    const std::vector<int> expectedVertices={3,2,1,4,0};
    for(std::size_t index=0;index<expectedVertices.size();++index) {
        check(vertexNames[index]==oldNames[expectedVertices[index]],"first alias occurrence and known/unknown priority");
        check(vertices[index].value[0]==expectedVertices[index],"vertex/name identity after reordering");
    }
    std::vector<std::string> mismatch={"phi_SI_V","n_e_SI_1pm3"};std::vector<int> one{17};
    const auto unchanged=mismatch;duplicates.order(mismatch,one);
    check(mismatch==unchanged && one==std::vector<int>{17},"mismatched field metadata must remain for serializer validation");
}
void checkOutputPreparation(const QuadtreeMesh& refined) {
    using namespace PicOutputReconstruction;
    FaceSegmentGraph graph;graph.build(refined);
    auto valid=buildActiveCellMask(refined,graph,nullptr);
    check(valid.size()==static_cast<std::size_t>(refined.cellIndexNum()),"output mask uses current cell indices");
    for(int ci=0;ci<refined.cellIndexNum();++ci) {
        const auto& cell=refined.cells(refined.cellIndexToLabel[ci]);
        check(valid[ci]==1 && cell.son[0]==-1 && cell.son[1]==-1 && cell.son[2]==-1 && cell.son[3]==-1,
            "valid refined quadtree leaf mask");
    }
    std::vector<double> alpha(valid.size(),1.);
    alpha[0]=-1.;alpha[1]=0.;alpha[2]=.25;alpha[3]=2.;
    alpha[4]=std::numeric_limits<double>::quiet_NaN();alpha[5]=std::numeric_limits<double>::infinity();
    auto filtered=buildActiveCellMask(refined,graph,nullptr,&alpha);
    for(int ci=0;ci<refined.cellIndexNum();++ci)check(filtered[ci]==(ci>=2),"alpha mask clamp/non-finite fallback");
    std::vector<double> wrongSize{0.};
    check(buildActiveCellMask(refined,graph,nullptr,&wrongSize)==valid,"wrong-size alpha must be ignored");
    FaceSegmentGraph empty;empty.build(refined,[](const QuadtreeMesh&,const Cell&) {return false;});
    auto noSegments=buildActiveCellMask(refined,empty,nullptr);
    for(auto value:noSegments)check(value==0,"no face segments means no output without cut metrics");
    Solvers::Poisson::CutCellMetrics inactive;
    auto inactiveMask=buildActiveCellMask(refined,graph,&inactive);
    for(auto value:inactiveMask)check(value==0,"cut activity must take priority over face graph");

    CellField<double> potential(refined,0.),prepared(refined,-17.);
    auto allMarked=markZeroPotentialCells(refined,potential,valid);
    for(bool value:allMarked)check(value,"boundary seeds must propagate across the zero cluster");
    fillMarkedCells(refined,valid,allMarked,prepared);
    for(int ci=0;ci<refined.cellIndexNum();++ci)check(prepared.value[ci]==-17.,"unsupported zero cluster must retain its prepared values");
    auto unmarked=markZeroPotentialCells(refined,potential,noSegments);
    for(bool value:unmarked)check(!value,"invalid cells must not seed the zero cluster");
    potential.value[0]=1e-300;
    auto withDonor=markZeroPotentialCells(refined,potential,valid);
    check(!withDonor[0],"tiny nonzero potential must remain a donor");
    prepared.value[0]=11.;fillMarkedCells(refined,valid,withDonor,prepared);
    for(int ci=0;ci<refined.cellIndexNum();++ci)check(prepared.value[ci]==11.,"constant donor must propagate across coarse/fine output cells");
    check(potential.value[0]==1e-300 && potential.value[1]==0.,"display filling modified physical potential");

    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph missing;missing.build(mesh,[](const QuadtreeMesh&,const Cell&) {return false;});
    DictionaryGeometry geometryInput;
    Line wall;wall.points[0]={.5,0.};wall.points[1]={.5,1.};wall.setType(LineType::straight);
    geometryInput.lines={wall};geometryInput.levelSetName={"wall"};geometryInput.innerPoints={{.25,.5}};
    geometryInput.linesOnInterface={{{0,false}}};geometryInput.moving={false};
    ImmersedGeometry<QuadtreeMesh> geometry(geometryInput,mesh);
    VertField<double> psi(mesh,1.),dx(mesh,0.),dy(mesh,0.);
    Solvers::Poisson::CutCellMetrics cut;Solvers::Poisson::CutCellMetrics::BuildOptions options;
    options.levelSetI=-1;options.psiOverride=&psi;options.psiDxOverride=&dx;options.psiDyOverride=&dy;
    options.isActive=[](const QuadtreeMesh&,const Cell& c) {return c.index()!=1;};
    cut.build(mesh,geometry,missing,options);
    auto activeWithoutFaces=buildActiveCellMask(mesh,missing,&cut);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)
        check(activeWithoutFaces[ci]==(ci!=1),"active cut cells must survive missing face segments");

    auto index=[&](double x,double y) {return mesh.cells(mesh.findCellLabel({x,y})).index();};
    CellField<double> neighbours(mesh,0.);neighbours.value[index(.75,.25)]=1.;
    neighbours.value[index(.25,.75)]=3.;neighbours.value[index(.75,.75)]=9.;
    std::vector<unsigned char> allValid(4,1);std::vector<bool> onlyCorner(4,false);onlyCorner[index(.25,.25)]=true;
    fillMarkedCells(mesh,allValid,onlyCorner,neighbours);
    check(neighbours.value[index(.25,.25)]==17./5.,"retain repeated edge-neighbour contributions through shared vertices");
    CellField<double> raw(mesh,0.),converted(mesh,0.);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)raw.value[ci]=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x<.5?-3.:1.;
    check(convertScalarValues(raw,converted,[](double v) {return v;},true)==2,"cell clamp count");
    auto before=buildMaskedVertField(mesh,converted,&allValid);
    auto after=buildMaskedVertField(mesh,raw,&allValid);
    check(convertScalarValues(*after,*after,[](double v) {return v;},true)>0,"in-place vertex clamp count");
    int center=0;
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[vi]).coord();
        if(p.x==.5 && p.y==.5) {check(before->value[vi]==.5 && after->value[vi]==0.,"cell-first and vertex-first output orders must remain distinct");++center;}
    }
    check(center==1 && raw.value[index(.25,.25)]==-3.,"conversion preserves physical source");

    DictionaryMesh longInput;longInput.Nx=1;longInput.Ny=24;longInput.maxCellNum=128;longInput.maxInitialIterations=0;
    QuadtreeMesh strip(longInput);CellField<double> longValues(strip,-1.);
    const int n=strip.cellIndexNum();check(n==24,"long propagation test mesh");
    std::vector<unsigned char> stripValid(n,1);std::vector<bool> stripMarked(n,true);
    stripMarked[n-1]=false;longValues.value[n-1]=7.;
    fillMarkedCells(strip,stripValid,stripMarked,longValues);
    for(int ci=0;ci<n;++ci)check(longValues.value[ci]==(ci<3?-1.:7.),"original 20-pass fill limit and traversal order");
}
void checkParticleFields(const QuadtreeMesh& mesh) {
    using namespace PicParticleFields;
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double inf=std::numeric_limits<double>::infinity();
    auto near=[](double actual,double expected) {
        check(std::isfinite(actual) && std::abs(actual-expected)<1e-11*(1.+std::abs(expected)),
            "particle display field differs from analytic amount/energy");
    };
    UnitsScales scales;scales.L0_m=2.;scales.T0_s=.5;
    auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(scales));
    NormalizationConfig config;config.slabThickness_m=3.;
    const auto norm=NormalizationContext::MakeFromUnits(units,config);
    CellField<double> amount(mesh,24.),density(mesh,0.);
    std::vector<double> alpha(mesh.cellIndexNum(),1.);
    alpha[0]=0.;alpha[1]=-1.;alpha[2]=.5;alpha[3]=2.;alpha[4]=inf;alpha[5]=nan;
    numberDensitySI(mesh,norm,amount,&alpha,density);
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const double area=mesh.cells(mesh.cellIndexToLabel[ci]).area();
        near(density.value[ci],ci<2?0.:(ci==2?4.:2.)/area);
        check(amount.value[ci]==24.,"density conversion modified deposited amount");
    }
    std::vector<double> wrongSize{0.};
    numberDensitySI(mesh,norm,amount,&wrongSize,density);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)
        near(density.value[ci],2./mesh.cells(mesh.cellIndexToLabel[ci]).area());

    ParticleType electrons("electron",2.*PhysicalConstants::eV_J,-PhysicalConstants::eV_J);
    electrons.particle.resize(3);
    const int first=mesh.cellIndexToLabel[0],second=mesh.cellIndexToLabel[1];
    auto& a=electrons.particle[0];a.w=2.;a.velocity={1.,2.};a.vz=2.;a.inCellL=first;
    a.position=mesh.cells(first).centerCoord();a.rngState=42;
    auto& b=electrons.particle[1];b.w=-1.;b.velocity={0.,0.};b.vz=1.;b.inCellL=-1;
    b.position=mesh.cells(second).centerCoord();
    auto& c=electrons.particle[2];c.w=3.;c.velocity={nan,0.};c.inCellL=first;
    c.position=mesh.cells(first).centerCoord();
    CouplingWeightsCache cache;
    cache.buildFromCSR(3,{0,2,3,4},{0,1,1,0},{.25,.75,1.,1.});
    CellField<double> number(mesh,0.),mean(mesh,0.),temperature(mesh,0.);
    depositNumberAmount(electrons,4.,cache,number);
    near(number.value[0],3.5);near(number.value[1],5.5);
    check(electronEnergyCSR(mesh,electrons,4.,norm,cache,ElectronWeightPolicy::FinitePositive,
        mean,temperature)==1,"explicit invalid-particle diagnostic count");
    near(mean.value[0],144.);near(mean.value[1],280./5.5);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)near(temperature.value[ci],(2./3.)*mean.value[ci]);
    CellField<double> hostMean(mesh,0.),hostTemperature(mesh,0.);
    electronEnergyHostCells(mesh,electrons,4.,norm,hostMean,hostTemperature);
    near(hostMean.value[0],144.);near(hostMean.value[1],16.);
    check(b.inCellL==-1 && a.rngState==42 && a.w==2. && a.velocity.x==1.,
        "display sampling changed particles, host identity or random state");

    CellField<double> strictMean(mesh,0.),strictTemperature(mesh,0.);
    check(electronEnergyCSR(mesh,electrons,inf,norm,cache,ElectronWeightPolicy::FinitePositive,
        strictMean,strictTemperature)==2,"finite policy must reject non-finite fallback weight");
    near(strictMean.value[1],144.);
    CellField<double> positiveMean(mesh,0.),positiveTemperature(mesh,0.);
    electronEnergyCSR(mesh,electrons,inf,norm,cache,ElectronWeightPolicy::Positive,
        positiveMean,positiveTemperature);
    near(positiveMean.value[0],144.);near(positiveMean.value[1],0.);
    ParticleType invalidMass("invalid",0.,-1.);invalidMass.particle=electrons.particle;
    CellField<double> invalidMean(mesh,0.),invalidTemperature(mesh,0.);
    check(electronEnergyCSR(mesh,invalidMass,4.,norm,cache,ElectronWeightPolicy::FinitePositive,
        invalidMean,invalidTemperature)==3,"invalid mass must reject all electron energy samples");
    for(int ci=0;ci<mesh.cellIndexNum();++ci)near(invalidMean.value[ci],0.);
    ParticleType empty("empty",1.,1.);CouplingWeightsCache emptyCache;
    emptyCache.buildFromCSR(0,{0},{},{});CellField<double> emptyAmount(mesh,0.);
    depositNumberAmount(empty,4.,emptyCache,emptyAmount);
    for(int ci=0;ci<mesh.cellIndexNum();++ci)near(emptyAmount.value[ci],0.);

    std::vector<ParticleType> all{electrons,empty};
    all[1].particle.resize(1);all[1].particle[0].w=-1.;all[1].particle[0].inCellL=first;
    CellField<double> sum(mesh,0.),average(mesh,0.),maximum(mesh,0.),count(mesh,0.);
    particleWeights(mesh,all,{4.},sum,average,maximum,count);
    near(sum.value[0],6.);near(average.value[0],2.);near(maximum.value[0],3.);near(count.value[0],3.);
    for(int ci=1;ci<mesh.cellIndexNum();++ci)near(count.value[ci],0.);
}
void checkReconstruction(const QuadtreeMesh& mesh) {
    using namespace PicOutputReconstruction;
    CellField<double> constant(mesh,17.);
    auto reconstructed=buildMaskedVertField(mesh,constant);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)
        check(std::abs(reconstructed->value[vi]-17.)<1e-12,"AMR reconstruction must preserve constants");
    std::vector<unsigned char> mask(mesh.cellIndexNum(),0);
    auto empty=buildMaskedVertField(mesh,constant,&mask);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)
        check(empty->value[vi]==0.,"unsupported vertices must be zero");
    CellField<double> masked(mesh,std::numeric_limits<double>::quiet_NaN());
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
        if(p.x<.5) {mask[ci]=1;masked.value[ci]=7.+p.x;}
    }
    const CellField<double> original(masked);
    auto bounded=buildMaskedVertField(mesh,masked,&mask);
    check(std::memcmp(original.value.get(),masked.value.get(),sizeof(double)*mesh.cellIndexNum())==0,
        "output reconstruction changed the physical input field");
    for(int vi=0;vi<mesh.vertIndexNum();++vi)
        check(std::isfinite(bounded->value[vi]) && bounded->value[vi]>=0. && bounded->value[vi]<=8.,
            "masked solid data leaked into reconstruction or fit exceeded bounds");
    for(int ci=0;ci<mesh.cellIndexNum();++ci)if(mask[ci])
        for(int label:mesh.cells(mesh.cellIndexToLabel[ci]).cellToVert) {
            const int vi=mesh.verts(label).index();
            check(bounded->value[vi]>=7.,"active cell corner lost its field at a coarse/fine junction");
        }
    std::vector<unsigned char> wrongSize{0};
    auto ignored=buildMaskedVertField(mesh,constant,&wrongSize);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)
        check(ignored->value[vi]==reconstructed->value[vi],"wrong-size reconstruction mask must be ignored");
    // One contributing cell exercises averaging/fallback without a fitted stencil.
    mask.assign(mask.size(),0);mask[0]=1;
    auto single=buildMaskedVertField(mesh,constant,&mask);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)
        check(single->value[vi]==0. || single->value[vi]==17.,"single-cell fallback changed its value");

    std::set<std::pair<double,double>> points;
    std::vector<VertField<double>> fields;fields.emplace_back(mesh,0.);fields.emplace_back(mesh,0.);
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[vi]).coord();points.insert({p.x,p.y});
        fields[0].value[vi]=10.+2.*p.x+4.*p.y;fields[1].value[vi]=p.x*p.y;
    }
    int unmatched=0;
    for(int flags=0;flags<4;++flags) {
        auto paired=fields;
        canonicalizePeriodicBoundaryVertices(mesh,(flags&1)!=0,(flags&2)!=0,paired);
        for(int vi=0;vi<mesh.vertIndexNum();++vi) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[vi]).coord();double x=p.x,y=p.y;
            if((flags&2) && (y==0. || y==1.)) {
                if(points.count({x,1.-y}))y=.5;
                else ++unmatched;
            }
            if((flags&1) && (x==0. || x==1.) && points.count({1.-x,p.y}))x=.5;
            check(paired[0].value[vi]==10.+2.*x+4.*y && paired[1].value[vi]==x*y,
                "periodic pairing, corner averaging or unmatched vertex changed");
        }
    }
    check(unmatched>0,"AMR test must include an unmatched periodic vertex");
    std::vector<VertField<double>> none;canonicalizePeriodicBoundaryVertices(mesh,true,true,none);

    DictionaryMesh regularInput;regularInput.Nx=regularInput.Ny=2;
    regularInput.maxCellNum=16;regularInput.maxInitialIterations=0;QuadtreeMesh regular(regularInput);
    CellField<double> affine(regular,0.);
    for(int ci=0;ci<regular.cellIndexNum();++ci) {
        const auto p=regular.cells(regular.cellIndexToLabel[ci]).centerCoord();
        affine.value[ci]=10.+2.*p.x+4.*p.y;
    }
    auto fitted=buildMaskedVertField(regular,affine);int checked=0;
    for(int vi=0;vi<regular.vertIndexNum();++vi) {
        const auto p=regular.verts(regular.vertIndexToLabel[vi]).coord();
        if(p.x==.5 && p.y==.5) {check(std::abs(fitted->value[vi]-13.)<1e-12,"interior affine reconstruction");++checked;}
        if(p.x==0. && p.y==0.) {check(fitted->value[vi]==11.5,"corner extrapolation must be bounded");++checked;}
    }
    check(checked==2,"analytic reconstruction sample vertices missing");
}
std::string bytes(const std::filesystem::path& path) {
    std::ifstream file(path,std::ios::binary);return {std::istreambuf_iterator<char>(file),{}};
}
void checkOutputDirectory() {
    using Directory=PicOutputDirectory;
    auto put=[](const std::filesystem::path& path,const char* text="keep") {
        std::filesystem::create_directories(path.parent_path());std::ofstream file(path);file<<text;
    };
    for(auto mode:{Directory::Mode::Explicit,Directory::Mode::Implicit}) {
        std::filesystem::path base=mode==Directory::Mode::Explicit?"output_directory_explicit":"output_directory_implicit";
        std::filesystem::path diagnosticBase="custom diagnostics/nested";
        const Directory output(mode,diagnosticBase,base);diagnosticBase="changed";base="changed";
        check(!std::filesystem::exists(output.root) && output.diagnostics==output.series/"custom diagnostics/nested",
            "output directory owns paths and construction has no filesystem side effects");
        output.create();output.create();
        check(std::filesystem::is_directory(output.fields) && std::filesystem::is_directory(output.series)
            && std::filesystem::is_directory(output.meta) && !std::filesystem::exists(output.diagnostics),
            "output creates only the three standard directories, idempotently");
        put(output.meta/"keep.txt");put(output.series/"keep.csv");
        put(output.fields/"frame_4.dat");put(output.fields/"frame_5.dat");put(output.fields/"nested/frame_6.dat");
        put(output.fields/"nested/arbitrary_5.txt");put(output.fields/"nested/unindexed.txt");
        output.prepareFields(5);
        check(std::filesystem::exists(output.fields/"frame_4.dat") && !std::filesystem::exists(output.fields/"frame_5.dat")
            && !std::filesystem::exists(output.fields/"nested/frame_6.dat")
            && !std::filesystem::exists(output.fields/"nested/arbitrary_5.txt")
            && bytes(output.fields/"nested/unindexed.txt")=="keep","restore pruning uses the supplied inclusive step threshold recursively, independent of extension");
        output.prepareFields(0);
        check(std::filesystem::is_directory(output.fields) && std::filesystem::is_empty(output.fields)
            && bytes(output.meta/"keep.txt")=="keep" && bytes(output.series/"keep.csv")=="keep",
            "fresh fields reset removes its full subtree and preserves other output categories");
        put(output.series/"current_balance.csv");put(output.series/"impact_report.csv");
        put(output.series/"coupling_fallbacks.csv");put(output.root/"reports/see_report.txt");
        put(output.meta/"see_report.txt");put(output.series/"impact_step.csv");
        put(output.series/"periodic_diagnostics.csv/child.txt");
        std::filesystem::create_directory(output.series/"current_waveform.csv");
        output.removeObsoleteReports();
        const bool explicitMode=mode==Directory::Mode::Explicit;
        check(!std::filesystem::exists(output.series/"current_balance.csv")
            && !std::filesystem::exists(output.series/"current_waveform.csv"),"common retired files and empty legacy directories removed");
        check(std::filesystem::exists(output.series/"impact_report.csv")==!explicitMode
            && std::filesystem::exists(output.root/"reports/see_report.txt")==!explicitMode
            && std::filesystem::exists(output.series/"coupling_fallbacks.csv")==explicitMode,
            "each integrator removes only its own retired names");
        check(bytes(output.series/"periodic_diagnostics.csv/child.txt")=="keep"
            && bytes(output.meta/"see_report.txt")=="keep" && bytes(output.series/"impact_step.csv")=="keep"
            && bytes(output.series/"keep.csv")=="keep","retired nonempty directories and current/unlisted files are preserved");
        const auto failedRoot=output.root/"failed";put(failedRoot/"meta");
        Directory failed(mode,"diagnostics",failedRoot);bool rejected=false;
        try {failed.create();}catch(const std::filesystem::filesystem_error& e){rejected=e.path1()==failed.meta;}
        check(rejected && std::filesystem::is_directory(failed.fields) && std::filesystem::is_directory(failed.series)
            && bytes(failed.meta)=="keep","directory creation order and failure preserve prior stages without rollback");
        const auto blockedRoot=output.root/"blocked";put(blockedRoot);
        Directory blocked(mode,"diagnostics",blockedRoot);rejected=false;
        try {blocked.create();}catch(const std::filesystem::filesystem_error& e){rejected=e.path1()==(explicitMode?blocked.root:blocked.fields);}
        check(rejected && bytes(blockedRoot)=="keep","explicit root creation and implicit fields-first errors remain distinct");
    }
    const auto absolute=std::filesystem::absolute("diagnostics_absolute");
    Directory absoluteDiagnostics(Directory::Mode::Explicit,absolute,"output_absolute");
    check(absoluteDiagnostics.diagnostics==absolute,"absolute diagnostics base retains existing filesystem join semantics");
}
void checkOutputManifest() {
    const std::filesystem::path root="manifest_output";
    std::filesystem::create_directories(root/"meta");
    const auto run=root/"meta/run_manifest.json",schema=root/"meta/schema_manifest.json";
    PicConfig cfg;
    cfg.output.every=2;cfg.output.v2.seriesEvery=3;cfg.output.diagnostics.every=5;
    cfg.output.v2.amountEvery=4;
    const PicOutputCadence explicitCadence(PicOutputCadence::Mode::Explicit,cfg,5);
    check(explicitCadence.fields(0) && explicitCadence.fields(2) && !explicitCadence.fields(5)
        && explicitCadence.explicitFinalFields(),"explicit final frame retains its separate stage");
    check(explicitCadence.series(0) && explicitCadence.series(3) && explicitCadence.series(5)
        && !explicitCadence.series(2) && !explicitCadence.periodicSeries(5),"series end and periodic-only cadence differ");
    check(explicitCadence.amounts(4) && explicitCadence.amounts(5) && !explicitCadence.amounts(3)
        && explicitCadence.explicitDetailFlush(50) && explicitCadence.explicitDetailFlush(5)
        && !explicitCadence.explicitDetailFlush(4),"amount and legacy buffer cadence");
    cfg.output.every=0;cfg.output.v2.seriesEvery=0;cfg.output.v2.amountEvery=0;
    const PicOutputCadence disabled(PicOutputCadence::Mode::Explicit,cfg,5);
    const PicOutputCadence implicitFallback(PicOutputCadence::Mode::Implicit,cfg,5);
    check(!disabled.fields(0) && !disabled.fields(5) && !disabled.explicitFinalFields()
        && implicitFallback.fields(1) && implicitFallback.fieldsEvery()==1
        && !implicitFallback.explicitFinalFields() && disabled.series(1) && !disabled.amounts(5),
        "original nonpositive frequency policies changed");
    check(explicitCadence.fieldsEvery()==2 && explicitCadence.seriesEvery()==3
        && explicitCadence.diagnosticsEvery()==5 && !explicitCadence.fields(-1),"cadence retained mutable config");
    cfg.output.every=2;cfg.output.v2.seriesEvery=3;
    const PicOutputCadence implicitCadence(PicOutputCadence::Mode::Implicit,cfg,5);
    check(implicitCadence.fields(4) && !implicitCadence.fields(5) && !implicitCadence.explicitFinalFields(),
        "implicit fields gained a forced final frame");
    cfg.species.resize(1);cfg.species[0].name="ion\"\\\n\r\t\x01";
    cfg.species[0].q=2.;cfg.species[0].m=3.;cfg.species[0].weight=4.;
    cfg.output.fields={"phi","E"};
    auto& electrode=cfg.output.diagnostics.electrodeCurrent;
    electrode.enabled=true;electrode.perBoundary=false;electrode.includeDomain=true;
    electrode.groups.push_back({"left\"electrode",{100,101}});
    const std::string dict="system\\pic\"special.dict";
    PicOutputManifest::write({cfg,root,dict,nullptr,nullptr},explicitCadence);
    auto text=bytes(run);const auto explicitSchema=bytes(schema);
    check(text.find("\"module\": \"ExplicitPIC\"")!=std::string::npos &&
        text.find("\"fieldsEvery\": 2")!=std::string::npos &&
        text.find("\"seriesEvery\": 3")!=std::string::npos &&
        text.find("\"diagnosticsEvery\": 5")!=std::string::npos,"explicit manifest cadence");
    check(text.find("ion\\\"\\\\\\n\\r\\t\\u0001")!=std::string::npos &&
        text.find("system\\\\pic\\\"special.dict")!=std::string::npos &&
        text.find("left\\\"electrode")!=std::string::npos &&
        text.find("\"boundaryIds\":[100,101]")!=std::string::npos,"manifest JSON strings/groups");
    check(text.find("\"q_C\":2,\"m_kg\":3,\"weight\":4")!=std::string::npos &&
        text.find("\"units\"")==std::string::npos && text.find("\"normalization\"")==std::string::npos,
        "manifest optional reference sections");
    check(explicitSchema.find("rG_L2_per_active")!=std::string::npos &&
        explicitSchema.find("overview_detail_step.csv")!=std::string::npos &&
        explicitSchema.find("physics_step.csv")==std::string::npos,"explicit schema identity");
    UnitsScales scales;scales.L0_m=2.;scales.T0_s=.5;scales.Phi0_V=3.;
    const auto units=UnitsContext::MakeSI(scales);
    NormalizationConfig normalization;normalization.dimension=2;normalization.slabThickness_m=.25;
    cfg.species.clear();cfg.output.fields.clear();electrode.groups.clear();
    PicOutputManifest::write({cfg,root,dict,&units,&normalization},implicitCadence);
    text=bytes(run);const auto implicitSchema=bytes(schema);
    check(text.find("\"module\": \"ImplicitPIC\"")!=std::string::npos &&
        text.find("\"perfEvery\": 3")!=std::string::npos && text.find("\"physicsEvery\": 3")!=std::string::npos &&
        text.find("\"speciesEvery\": 3")!=std::string::npos && text.find("\"amountEvery\": 0")!=std::string::npos &&
        text.find("\"boundaryEvery\": 3")!=std::string::npos && text.find("\"flushEvery\": 3")!=std::string::npos,
        "implicit manifest must use the actual shared cadence");
    check(text.find("\"L0_m\": 2")!=std::string::npos && text.find("\"T0_s\": 0.5")!=std::string::npos &&
        text.find("\"slabThickness_m\": 0.25")!=std::string::npos && text.find("ion\\\"")==std::string::npos &&
        text.find("\"selectedFieldTokens\": []")!=std::string::npos,"manifest must not retain earlier input");
    check(implicitSchema.find("rG_L2_per_active")==std::string::npos &&
        implicitSchema.find("physics_step.csv")!=std::string::npos &&
        implicitSchema.find("\"name\":\"keep_prob\"")!=std::string::npos,"implicit schema identity");
    // One metadata file's open failure must not suppress the other file.
    struct Capture {
        std::ostringstream stream;std::streambuf* previous=std::cerr.rdbuf(stream.rdbuf());
        ~Capture(){std::cerr.rdbuf(previous);}
    } capture;
    std::filesystem::remove(run);std::filesystem::create_directory(run);
    PicOutputManifest::write({cfg,root,dict,nullptr,nullptr},explicitCadence);
    check(bytes(schema)==explicitSchema,"run-manifest failure suppressed schema writer");
    check(capture.stream.str().find("[ExplicitPIC] WARNING: failed to create run manifest:")!=std::string::npos,
        "manifest open warning must identify integrator and file");
    std::filesystem::remove(run);std::filesystem::remove(schema);std::filesystem::create_directory(schema);
    PicOutputManifest::write({cfg,root,dict,&units,&normalization},implicitCadence);
    check(bytes(run)==text && capture.stream.str().find("[ImplicitPIC] WARNING: failed to create schema manifest:")!=std::string::npos,
        "schema failure changed completed run manifest or suppressed warning");
}
std::vector<std::vector<std::string>> csvRows(const std::filesystem::path& path) {
    std::ifstream file(path);check(bool(file),"missing species CSV");
    std::vector<std::vector<std::string>> rows;std::string line;
    while(std::getline(file,line)) {
        std::istringstream stream(line);std::vector<std::string> row;std::string field;
        while(std::getline(stream,field,','))row.push_back(field);
        rows.push_back(std::move(row));
    }
    return rows;
}
void checkCutcellStateOutput() {
    using View=PicCutcellStateView;using Metrics=Solvers::Poisson::CutCellMetrics;
    struct Logs {std::ostringstream out,err;std::streambuf* o=std::cout.rdbuf(out.rdbuf());
        std::streambuf* e=std::cerr.rdbuf(err.rdbuf());~Logs(){std::cout.rdbuf(o);std::cerr.rdbuf(e);}} logs;
    DictionaryMesh settings;settings.Nx=settings.Ny=2;settings.maxCellNum=64;settings.maxInitialIterations=0;
    QuadtreeMesh mesh(settings);CellField<double> rho(mesh,2.),phi(mesh,6.);CellField<Vector2D> electric(mesh,{4.,-2.});
    std::vector<double> kappa{3.};
    // Controlled metrics make the three boundary closures independently analytic.
    Metrics metrics;auto& cells=const_cast<std::vector<Metrics::CellMetric>&>(metrics.cellMetricsByIndex());
    auto& active=const_cast<std::vector<unsigned char>&>(metrics.activeMaskByCellIndex());
    cells.resize(4);active={1,1,1,0};
    for(int i=0;i<4;++i){auto& m=cells[i];m.alpha=.5;m.A_gamma=.5;m.d_gamma=.25;
        m.n_gamma={1.,0.};m.isCut=true;m.gammaBoundaryId=i+1;m.gammaLineIndex=10+i;}
    PicEmbeddedBC fallback;fallback.enabled=true;fallback.type=PicBCType::Robin;fallback.a=2.;fallback.b=.5;fallback.g=1.;
    PicEmbeddedBCMap byId;byId[1]=fallback;byId[1].type=PicBCType::Dirichlet;byId[1].g=2.;
    byId[2]=fallback;byId[2].type=PicBCType::Neumann;byId[2].g=3.;byId[3]=fallback;byId[3].enabled=false;
    std::vector<LevelSet> levels(2);for(auto& ls:levels)ls.initialize(mesh,0.);
    for(int vi=0;vi<mesh.vertIndexNum();++vi){const auto p=mesh.verts(mesh.vertIndexToLabel[vi]).coord();
        levels[0].value[vi]=p.x+2.*p.y;levels[1].value[vi]=3.*p.x-p.y;}
    UnitsScales scales;scales.L0_m=2.;scales.Phi0_V=6.;
    const auto units=std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));
    const auto norm=NormalizationContext::MakeFromUnits(units,{});
    View::Input input{mesh,rho,kappa,phi,electric,&metrics,&fallback,&byId,&levels,units.get(),&norm};
    {
        const View view(input);check(view.size()==4 && view.levelSetCount()==2 && !view.active(3),"cut-cell view layout/mask");
        for(int i=0;i<4;++i){const auto r=view.sample(i);const auto center=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
            check(r.ci==i && r.cLabel==mesh.cellIndexToLabel[i] && r.x_m==2.*center.x && r.y_m==2.*center.y
                && r.areaCode==.25 && r.area_m2==1. && r.aGamma_m==1. && r.dGamma_m==.5,"cut-cell SI geometry");
            check(r.rhoCode==2. && std::abs(r.rhoSI-2.*8.8541878128e-12*6./4.)<1e-25
                && r.phiV==36. && r.exVpm==12. && r.eyVpm==-6. && r.kappa==(i==0?3.:1.),"cut-cell physical scales/default coefficient");
            const double expected[]={16.,-3.,11.,0.};
            check(r.enGammaCode==expected[i] && r.enGammaVpm==3.*expected[i],"Dirichlet/Neumann/Robin normal-field closures and inactive cell");
            check(std::string(r.ebSource)==(i<2?"map":"fallback") && r.eb==(i<2?&byId.at(i+1):&fallback),"enabled map priority/disabled entry fallback");
            check(r.psi==std::vector<double>{center.x+2.*center.y,3.*center.x-center.y},"bilinear level-set sampling at actual cell center");
        }
    }
    const std::filesystem::path root="cutcell_output";
    PicFieldOutput::writeCutcellState(root,{},7,.125,false,input);
    const auto path=root/"cutcell_state_step_000007.csv";const auto all=csvRows(path);
    const std::string header="step,time_s,cellIndex,cellLabel,level,boundary,center_x_code,center_y_code,center_x_m,center_y_m,area_code2,area_m2,active,isCut,alpha,A_gamma_code,A_gamma_m,d_gamma_code,d_gamma_m,n_gamma_x,n_gamma_y,gammaBoundaryId,gammaLineIndex,rho_code,rho_Cpm3,kappa_cell,phi_code,phi_V,Ex_code,Ey_code,Ex_Vpm,Ey_Vpm,EnGamma_code,EnGamma_Vpm,eb_bc_source,eb_bc_type,eb_bc_a,eb_bc_b,eb_bc_g,psi_ls01_code,psi_ls02_code";
    check(bytes(path).rfind(header,0)==0 && all.size()==5,"complete cut-cell header and all-cell rows");
    for(int i=0;i<4;++i){const auto& row=all[i+1];check(row.size()==41 && row[0]=="7" && row[1]=="0.1250000000000000"
        && row[2]==std::to_string(i) && row[27]=="36.0000000000000000" && row[30]=="12.0000000000000000"
        && row[34]==(i<2?"map":"fallback") && row[35]==(i==0?"Dirichlet":(i==1?"Neumann":"Robin")),"cut-cell column order, precision and boundary names");}
    PicFieldOutput::writeCutcellState(root,"active.csv",7,.125,true,input);
    const auto selected=csvRows(root/"active.csv");check(selected.size()==4,"active-only mask filters one cell");
    for(int i=0;i<4;++i)check(selected[i]==all[i],"active filtering preserves selected rows");
    const std::vector<std::filesystem::path> paths{"./dot.csv","cutcell_nested/nested.csv",std::filesystem::absolute("cutcell absolute/space.csv")};
    for(const auto& p:paths){PicFieldOutput::writeCutcellState(root,p,7,.125,false,input);
        const auto resolved=p.parent_path()=="."?root/p:p;check(bytes(resolved)==bytes(path),"legacy relative/absolute dump path resolution");}
    std::ofstream(root/"active.csv")<<"old tail must disappear";
    PicFieldOutput::writeCutcellState(root,"active.csv",7,.125,true,input);
    check(csvRows(root/"active.csv")==selected,"dump truncates and flushes synchronously");
    {
        View::Input plain{mesh,rho,kappa,phi,electric};const View view(plain);const auto r=view.sample(3);
        check(r.active && r.cm.alpha==1. && !r.eb && std::string(r.ebSource)=="none" && r.psi.empty()
            && r.phiV==6. && r.exVpm==4. && r.rhoSI==2. && r.enGammaCode==0.,"null metrics/geometry/units legacy defaults");
        CellField<double> empty;CellField<Vector2D> emptyE;std::vector<double> emptyK;
        const auto shortRow=View({mesh,empty,emptyK,empty,emptyE}).sample(0);
        check(shortRow.phiCode==0. && shortRow.rhoCode==0. && shortRow.eCode.x==0. && shortRow.kappa==1.,"short field arrays use legacy defaults");
        for(auto type:{PicBCType::Periodic,static_cast<PicBCType>(99)}){
            PicEmbeddedBC bc=fallback;bc.type=type;plain.embedded=&bc;
            PicFieldOutput::writeCutcellState(root,"names.csv",1,0.,false,plain);
            check(csvRows(root/"names.csv")[1][35]==(type==PicBCType::Periodic?"Periodic":"Unknown"),"remaining boundary type labels");}
    }
    // Missing mapped boundary fails before truncation; skip keeps the diagnostic.
    auto missing=input;missing.embedded=nullptr;const auto before=bytes(path);bool failed=false;
    try{PicFieldOutput::writeCutcellState(root,{},7,.125,false,missing);}catch(const std::runtime_error&){failed=true;}
    check(failed && bytes(path)==before,"normal-field failure preserves existing dump");
    missing.skipMissingEmbeddedBC=true;check(View(missing).sample(2).enGammaCode==0.,"missing-boundary skip retains zero normal field");
    PicFieldOutput::writeCutcellState(root,"skip.csv",7,.125,false,missing);
    check(csvRows(root/"skip.csv")[3][34]=="none","skip serializes absent boundary source");
    std::filesystem::create_directories(root/"blocked.csv");
    PicFieldOutput::writeCutcellState(root,"blocked.csv",7,.125,false,input);
    check(logs.err.str().find("failed to create cutcell state dump")!=std::string::npos,"open failure retains warning and nonfatal return");
    std::ofstream("cutcell_parent_file")<<"file";failed=false;
    try{PicFieldOutput::writeCutcellState(root,"cutcell_parent_file/child.csv",7,.125,false,input);}
    catch(const std::filesystem::filesystem_error&){failed=true;}
    check(failed && bytes("cutcell_parent_file")=="file","parent creation errors still propagate");
    check(logs.out.str().find("at step=7")!=std::string::npos && phi.value[0]==6. && rho.value[0]==2.
        && electric.value[0].x==4. && kappa==std::vector<double>{3.} && cells[0].d_gamma==.25,"diagnostic reports original stage and does not modify inputs");
}
void checkAmrSeries() {
    using Series=PicStepSeries;
    const std::string dir="amr_series/";
    ExplicitPicRefinement::Report indicator{.125,.25,.5,.75,7,9,2,3,1e12};
    ImplicitPicRefinement::AdaptationReport adaptation;
    adaptation.cellsBefore=4;adaptation.cellsAfter=16;adaptation.useFastRemap=true;
    adaptation.meshAdaptS=4.125;adaptation.movingEbMeshChangedS=5.125;
    adaptation.faceGraphS=14.125;adaptation.embeddedMetricsS=15.125;adaptation.fieldResetS=16.125;
    adaptation.relabelParticlesS=17.125;adaptation.fieldRebuildS=18.125;
    adaptation.coordRefreshS=19.125;adaptation.seconds=20.125;
    PicAmrGeometry::Report report;report.geometryRefreshS=6.125;report.geometryReinitSkipped=true;
    report.geomInitializeLevelSetsS=8.125;report.geomPdeReinitS=9.125;report.geomRefreshDerivedS=10.125;
    report.geomCheckCutInfoS=11.125;report.geomRebuildCompositeS=12.125;report.geomMarkPhasesS=13.125;
    report.amrRawActiveMovingLevelSets={-1,0,0,8};report.amrSelectiveMovingLevelSets={1,8};
    report.amrHybridSensitiveInactiveLevelSets={0,0,-1};
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryGeometry dict;
    dict.levelSetName={"fixed","moving"};dict.innerPoints={{.1,.5},{.9,.5}};
    dict.moving={false,true};dict.linesOnInterface={{},{}};
    PicAmrGeometry::Geometry geometry(dict,mesh);
    auto& diagnostics=const_cast<std::vector<PicAmrGeometry::Geometry::LevelSetUpdateDiagnostics>&>(geometry.levelSetUpdateDiagnostics());
    diagnostics.resize(2);
    diagnostics[0]={-2.,4.,5,6,7,8,.125,.25,true};
    diagnostics[1]={-3.,9.,10,11,12,13,.5,.75,false};
    const auto open=[&](Series& e,Series& i,int step){e.openAmrIndicators(dir+"indicator.csv",step);
        i.openAmrAdaptation(dir+"adapt.csv",step);i.openAmrLevelSets(dir+"levelsets.csv",step);};
    const auto write=[&](Series& e,Series& i,int step){e.writeAmrIndicators(step,.0625,indicator,true,false);
        i.writeAmrAdaptation(step,adaptation,report);i.writeAmrLevelSets(step,&geometry,report);};
    {
        Series e(Series::Mode::Explicit),i(Series::Mode::Implicit);
        check(!e.amrIndicatorsReady() && !i.amrAdaptationReady() && !i.amrLevelSetsReady(),"unopened AMR streams unavailable");
        write(e,i,0);check(!std::filesystem::exists(dir),"disabled AMR writes do not create files");
        open(e,i,0);check(e.amrIndicatorsReady() && i.amrAdaptationReady() && i.amrLevelSetsReady(),"AMR streams opened");
        check(csvRows(dir+"adapt.csv").size()==1 && csvRows(dir+"levelsets.csv").size()==1,"implicit AMR headers immediately flushed");
        i.writeAmrLevelSets(0,nullptr,report);check(csvRows(dir+"levelsets.csv").size()==1,"absent geometry writes no rows");
        for(int step:{1,3,5})write(e,i,step);
        check(csvRows(dir+"indicator.csv").size()==4 && csvRows(dir+"adapt.csv").size()==4
            && csvRows(dir+"levelsets.csv").size()==7,"all AMR rows immediately visible");
        e.closeAmrIndicators();const auto closed=bytes(dir+"indicator.csv");
        e.writeAmrIndicators(6,1.,indicator,false,true);
        check(!e.amrIndicatorsReady() && bytes(dir+"indicator.csv")==closed,"closing indicators disables subsequent output");
    }
    auto indicators=csvRows(dir+"indicator.csv"),adaptations=csvRows(dir+"adapt.csv"),levels=csvRows(dir+"levelsets.csv");
    check(indicators[0]==std::vector<std::string>{"step","time_s","lambdaD_min_m","h_min_m","ratio_max","ratio_mean",
        "refine_cells","coarsen_cells","pending","adapted"},"explicit AMR schema");
    check(indicators[1]==std::vector<std::string>{"1","0.0625","0.125","0.25","0.5","0.75","7","9","1","0"},"explicit AMR defaultfloat values and flags");
    check(adaptations[0]==std::vector<std::string>{"step","cells_before","cells_after","use_fast_remap","mesh_adapt_s",
        "movingeb_onMeshChanged_s","geometry_refresh_s","geometry_reinit_skipped","geom_initialize_levelsets_s","geom_pde_reinit_s",
        "geom_refresh_derived_s","geom_check_cutinfo_s","geom_rebuild_composite_s","geom_mark_phases_s","face_graph_s",
        "embedded_metrics_s","field_buffer_reset_s","relabel_particles_s","field_rebuild_s","coord_refresh_s","total_adapt_s"},"implicit AMR adaptation schema");
    check(adaptations[1].size()==21 && adaptations[1][1]=="4" && adaptations[1][2]=="16"
        && adaptations[1][3]=="1" && adaptations[1][7]=="1","implicit AMR integer counts and flags");
    for(int k=4;k<21;++k)if(k!=7)
        check(adaptations[1][k]==std::to_string(k)+".1250000000000000","implicit AMR timing order and fixed precision");
    check(levels[0]==std::vector<std::string>{"step","level_set_index","moving","active_moving_levelset","selective_reinit_target",
        "hybrid_sensitive_inactive_target","reinitialized","reinit_iterations","reinit_fine_residual","reinit_coarse_residual",
        "min_value","max_value","near_zero_count","nonfinite_count","narrow_band_count"},"AMR level-set schema");
    check(levels[1]==std::vector<std::string>{"1","1","0","1","0","1","1","8","0.1250000000000000",
        "0.2500000000000000","-2.0000000000000000","4.0000000000000000","5","6","7"},"fixed level-set diagnostics and deduplicated valid masks");
    check(levels[2]==std::vector<std::string>{"1","2","1","0","1","0","0","13","0.5000000000000000",
        "0.7500000000000000","-3.0000000000000000","9.0000000000000000","10","11","12"},"moving level-set diagnostics use one-based index");
    {
        Series e(Series::Mode::Explicit),i(Series::Mode::Implicit);open(e,i,3);write(e,i,3);
        check(csvRows(dir+"indicator.csv").size()==3 && csvRows(dir+"adapt.csv").size()==3
            && csvRows(dir+"levelsets.csv").size()==5,"AMR restart trims all rows at/after resume without duplicating header");
        diagnostics.resize(1);i.writeAmrLevelSets(4,&geometry,report);
        check(csvRows(dir+"levelsets.csv").size()==6,"short diagnostic vector clamps level-set loop");
        diagnostics.resize(3);i.writeAmrLevelSets(5,&geometry,report);
        check(csvRows(dir+"levelsets.csv").size()==8,"extra diagnostics do not fabricate level sets");
    }
    {Series e(Series::Mode::Explicit),i(Series::Mode::Implicit);open(e,i,0);}
    for(const auto& file:{"indicator.csv","adapt.csv","levelsets.csv"})
        check(csvRows(dir+file).size()==1,"fresh AMR session replaces old history");
    {std::ofstream blocker(dir+"blocked");blocker<<"preserve";}
    Series failed(Series::Mode::Explicit);failed.openAmrIndicators(dir+"blocked/indicator.csv",0);
    failed.openAmrAdaptation(dir+"blocked/adapt.csv",0);failed.openAmrLevelSets(dir+"blocked/levelsets.csv",0);
    check(!failed.amrIndicatorsReady() && !failed.amrAdaptationReady() && !failed.amrLevelSetsReady()
        && bytes(dir+"blocked")=="preserve","AMR open failure leaves unavailable streams and preserves unrelated files");
    struct Flushed : std::stringbuf {int count=0;int sync() override {++count;return 0;}} buffer;
    std::ostream console(&buffer);
    ExplicitPicRefinement::AdaptationReport explicitAdaptation{0.,{4,8,3},{1,12}};
    Series::writeAmrConsole(console,5,13,.125,explicitAdaptation);
    check(buffer.str()=="[AMR] step=5 cells=13 add=4 del=6 hpMinLeaf=1.250000e-01 levels={0:1,1:12}\n"
        && buffer.count==1 && console.precision()==6 && (console.flags()&std::ios::floatfield)==0,"explicit console level deltas, flush and flags");
    buffer.str("");Series::writeAmrConsole(console,6,0,.25,{});
    check(buffer.str()=="[AMR] step=6 cells=0 add=0 del=0 hpMinLeaf=2.500000e-01 levels={}\n","empty explicit level counts");
    buffer.str("");Series::writeAmrConsole(console,7,adaptation,report);
    check(buffer.str()=="[ImplicitPIC][Event][AMR] step=7 cells 4->16 adapt_s=20.125 fieldRebuild=remap+evalE geom[pde=9.125,derived=10.125,cut=11.125,comp=12.125,mark=13.125]\n"
        && buffer.count==3 && console.precision()==6 && (console.flags()&std::ios::floatfield)==0,"implicit console format, flush and flags");
    buffer.str("");adaptation.useFastRemap=false;Series::writeAmrConsole(console,8,adaptation,report);
    check(buffer.str().find("fieldRebuild=deposit+solve")!=std::string::npos,"implicit full-rebuild console mode");
    buffer.str("");Series::writeAmrDebyeWarning(console,{});check(buffer.str().empty() && buffer.count==4,"no empty Debye warning or flush");
    Series::writeAmrDebyeWarning(console,indicator);
    check(buffer.str()=="[AMR][Debye] lowNe=2 badLambdaD=3 (neFloor=1.000e+12)\n"
        && buffer.count==5 && console.precision()==6 && (console.flags()&std::ios::floatfield)==0,"Debye warning format, flush and flags");
}
void checkGeometryReports() {
    using Series=PicStepSeries;const std::filesystem::path base="geometry_reports",series=base/"series";
    const auto metrics=series/"geometry_step.csv";
    LegacyGridMotionState state;state.stageName="stage1_hot_gap";state.stage1Alpha=.25;state.stage2Alpha=.5;
    state.gsa_m=.001;state.gad_m=.002;state.da_m=.003;state.accelShift_m=.004;
    state.decelShift_m=.005;state.accelRecession_m=.006;state.screenDownstreamX_m=.007;
    state.accelUpstreamX_m=.008;state.accelDownstreamLocalX_m=.009;state.decelUpstreamX_m=.010;
    {Series writer(Series::Mode::Implicit);
     writer.writeGeometryMetrics(0,0.,state,true);
     writer.openGeometryMetrics(series,{},0,3);
     check(csvRows(metrics).size()==1,"geometry metrics header is immediately visible");
     writer.writeGeometryMetrics(0,0.,state,true);writer.writeGeometryMetrics(1,.25,state,false);
     check(csvRows(metrics).size()==2,"unchanged noncadence geometry sample skipped");
     state.changed=true;writer.writeGeometryMetrics(2,.5,state,false);state.changed=false;
     writer.writeGeometryMetrics(3,.75,state,false);writer.writeGeometryMetrics(4,1.,state,false);
     writer.writeGeometryMetrics(5,1.25,state,true);
     check(csvRows(metrics).size()==5,"geometry event, cadence and caller force rows flush immediately");}
    const auto data=csvRows(metrics);
    check(data[0].size()==15 && data[2][0]=="2" && data[3][0]=="3" && data[4][0]=="5"
        && data[2][1]=="0.5000000000000000" && data[2][2]=="stage1_hot_gap"
        && std::stod(data[2][3])==.25 && std::stod(data[2][4])==.5,"geometry schema, state and fixed precision");
    for(int i=5;i<15;++i)check(std::abs(std::stod(data[2][i])-(i-4.))<1e-13,"geometry SI-to-mm column order");
    {Series writer(Series::Mode::Implicit);writer.openGeometryMetrics(series,{},3,3);writer.writeGeometryMetrics(3,.75,state,true);}
    const auto resumed=csvRows(metrics);
    check(resumed.size()==4 && resumed[1][0]=="0" && resumed[2][0]=="2" && resumed[3][0]=="3","geometry CSV resume threshold is supplied by integrator");
    const std::vector<std::filesystem::path> requested={"named.csv","./dot.csv",base/"nested/metrics.csv",std::filesystem::absolute(base/"absolute.csv")};
    const std::vector<std::filesystem::path> expected={series/"named.csv",series/"dot.csv",base/"nested/metrics.csv",std::filesystem::absolute(base/"absolute.csv")};
    for(std::size_t i=0;i<requested.size();++i) {
        Series writer(Series::Mode::Implicit);writer.openGeometryMetrics(series,requested[i],0,0);
        writer.writeGeometryMetrics(1,.25,state,false);
        check(csvRows(expected[i]).size()==2,"geometry basename, dot, nested and absolute path resolution; cadence clamp");
    }
    {Series writer(Series::Mode::Implicit);writer.openGeometryMetrics(series,{},0,3);}
    check(csvRows(metrics).size()==1,"fresh geometry metrics replaces prior history");
    auto text=[](const std::filesystem::path& path) {std::ifstream file(path);return std::string(std::istreambuf_iterator<char>(file),{});};
    const auto trace=base/"trace.log";{std::ofstream file(trace);file<<"seed\n";}
    std::ostringstream console;auto* previousOut=std::cout.rdbuf(console.rdbuf());
    {Series writer(Series::Mode::Implicit);writer.openGeometryTrace(trace,3,false);
     writer.writeGeometryTrace(0,"skip");writer.writeGeometryTrace(1,"first");writer.writeGeometryTrace(5,"last");writer.writeGeometryTrace(6,"skip");}
    std::cout.rdbuf(previousOut);
    const std::string window="[moving-trace] step=1 first\n[moving-trace] step=5 last\n";
    check(console.str()==window && text(trace)=="seed\n"+window,"trace window includes both endpoints and appends without pruning");
    console.str("");previousOut=std::cout.rdbuf(console.rdbuf());
    {Series writer(Series::Mode::Implicit);writer.openGeometryTrace(trace,-1,true);writer.writeGeometryTrace(100,"all");}
    {Series writer(Series::Mode::Implicit);writer.openGeometryTrace({},-1,true);writer.writeGeometryTrace(7,"console only");}
    {Series writer(Series::Mode::Implicit);writer.openGeometryTrace(base/"disabled.log",-1,false);writer.writeGeometryTrace(7,"disabled");}
    std::cout.rdbuf(previousOut);
    check(console.str()=="[moving-trace] step=100 all\n[moving-trace] step=7 console only\n"
        && text(trace)=="seed\n"+window+"[moving-trace] step=100 all\n"
        && std::filesystem::exists(base/"disabled.log") && bytes(base/"disabled.log").empty(),"all-step, console-only and disabled trace file policies");
    std::filesystem::create_directories(base/"occupied/child");
    std::ostringstream warning;auto* previousError=std::cerr.rdbuf(warning.rdbuf());
    {Series writer(Series::Mode::Implicit);writer.openGeometryMetrics(series,base/"occupied",0,1);writer.writeGeometryMetrics(0,0.,state,true);}
    std::cerr.rdbuf(previousError);
    check(warning.str().find("[ImplicitPIC] WARNING: failed to create geometry metrics csv: ")==0,"geometry open failure retains warning");
    console.str("");previousOut=std::cout.rdbuf(console.rdbuf());
    {Series writer(Series::Mode::Implicit);writer.openGeometryTrace(base/"occupied",-1,true);writer.writeGeometryTrace(4,"open failed");}
    std::cout.rdbuf(previousOut);
    check(console.str()=="[moving-trace] step=4 open failed\n","failed trace file still emits console diagnostics");
    {std::ofstream file(base/"blocked");file<<"keep";}
    bool metricsThrow=false,traceThrow=false;
    try {Series writer(Series::Mode::Implicit);writer.openGeometryMetrics(series,base/"blocked/file.csv",0,1);}catch(const std::filesystem::filesystem_error&){metricsThrow=true;}
    try {Series writer(Series::Mode::Implicit);writer.openGeometryTrace(base/"blocked/file.log",-1,false);}catch(const std::filesystem::filesystem_error&){traceThrow=true;}
    check(metricsThrow && traceThrow && bytes(base/"blocked")=="keep","geometry and trace directory errors propagate without changing blocker");
}
void checkExplicitReports() {
    DictionaryMesh settings;settings.Nx=settings.Ny=2;settings.maxCellNum=64;settings.maxInitialIterations=0;
    QuadtreeMesh mesh(settings);ParticleType species("weighted",1.,-1.);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const std::array<double,6> weights{1.,2.,0.,-1.,nan,10.};
    const std::array<int,6> hosts{0,0,1,1,-1,2};
    for(std::size_t i=0;i<weights.size();++i) {
        Particle p;p.w=weights[i];p.inCellL=hosts[i]<0?-1:mesh.cellIndexToLabel[hosts[i]];
        p.rngState=100+i;p.trackId=200+i;species.particle.push_back(p);
    }
    const auto sample=summarizePicParticleWeights(mesh,species,3.);
    check(sample.macro==6 && sample.physical==22. && sample.minimum==1. && sample.maximum==10.
        && sample.mean==22./6. && sample.percentiles==std::array<double,3>{3.,10.,10.},
        "weights use positive actual values, fallback and nearest-rank percentiles");
    check(sample.cellMacroMaximum==2 && sample.cellMacroIndex==0 && sample.cellPhysicalMaximum==10.
        && sample.cellPhysicalIndex==2,"occupancy counts valid hosts only and keeps first maximum on ties");
    for(std::size_t i=0;i<weights.size();++i)
        check(species.particle[i].rngState==100+i && species.particle[i].trackId==200+i
            && (i==4?std::isnan(species.particle[i].w):species.particle[i].w==weights[i]),"weight sampling preserves particle order, weights and RNG");
    const int parent=mesh.cellIndexToLabel[0];DictionaryMesh::RefineBox box;
    box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
    species.particle.resize(1);species.particle[0].w=0.;species.particle[0].inCellL=parent;
    const int savedIndex=mesh.cells(parent).index();mesh.cells(parent).setIndex(-1);
    const auto invalid=summarizePicParticleWeights(mesh,species,-2.);
    check(invalid.macro==1 && invalid.physical==-2. && invalid.minimum==-2. && invalid.maximum==0.
        && invalid.cellMacroMaximum==0 && invalid.cellMacroIndex==-1 && invalid.cellPhysicalIndex==-1,
        "legacy negative fallback and invalid host index are preserved without repair");
    mesh.cells(parent).setIndex(mesh.cellIndexNum());
    check(summarizePicParticleWeights(mesh,species,1.).cellMacroIndex==-1,"upper out-of-range host index is ignored");
    mesh.cells(parent).setIndex(0);
    check(summarizePicParticleWeights(mesh,species,1.).cellMacroIndex==0,"report accepts an in-range stored index without leaf repair");
    mesh.cells(parent).setIndex(savedIndex);
    species.particle[0].w=7.;species.particle[0].inCellL=mesh.cellIndexToLabel[mesh.cellIndexNum()-1];
    check(summarizePicParticleWeights(mesh,species,1.).cellPhysicalIndex==mesh.cellIndexNum()-1,"occupancy uses current AMR indices");
    species.particle.clear();const auto empty=summarizePicParticleWeights(mesh,species,1.);
    check(empty.macro==0 && empty.physical==0. && empty.minimum==0. && empty.cellMacroIndex==-1
        && empty.cellPhysicalIndex==-1,"empty species report has zero values and absent maxima");
    const double ev=PhysicalConstants::eV_J;
    auto energy=summarizePicEmissionEnergy(4.,12.*ev,40.*ev*ev);
    check(std::abs(energy.meanEV-3.)<1e-14 && std::abs(energy.stdEV-1.)<1e-14,"physical-weighted SEE mean and standard deviation");
    check(summarizePicEmissionEnergy(1.,2.*ev,ev*ev).stdEV==0.
        && summarizePicEmissionEnergy(0.,nan,nan).meanEV==0.
        && summarizePicEmissionEnergy(-1.,nan,nan).stdEV==0.,"SEE variance clamp and nonpositive count gates");
    auto nonfinite=summarizePicEmissionEnergy(1.,nan,nan);
    check(std::isnan(nonfinite.meanEV) && nonfinite.stdEV==0.,"SEE legacy nonfinite arithmetic retained");

    PicDomainBC domain;PicEmbeddedBC embedded;
    domain.left.type=PicBCType::Neumann;
    domain.left.timeDirichlet={true,"sine",2.,3.,1.,0.,0.};
    domain.right.timeDirichlet={false,"sine",9.,1.,1.,0.,0.};
    domain.bottom.timeDirichlet={true,"sine",9.,1.,0.,0.,0.};
    domain.top.timeDirichlet={true,"custom",4.,2.,1.,0.,0.};
    embedded.enabled=true;embedded.timeDirichlet={true,"constant",6.,10.,1.,2.,0.};
    UnitsScales scales;scales.Phi0_V=10.;const auto units=UnitsContext::MakeSI(scales);
    const auto phase=summarizePicRfPhases(domain,embedded,.25,nullptr);
    const auto siPhase=summarizePicRfPhases(domain,embedded,.25,&units);
    check(phase.size()==3 && phase[0].target=="left" && phase[1].target=="top" && phase[2].target=="embedded",
        "RF report retains side order, default embedded scope and waveform gates independent of BC type");
    check(phase[0].potential==5. && phase[1].potential==6. && phase[2].potential==6.
        && phase[0].degrees==90. && phase[2].radians==0. && siPhase[0].potential==50.,
        "RF shares the waveform evaluator and converts only the reported potential");
    domain.left.timeDirichlet.type="changed";embedded.enabled=false;
    check(phase[0].type=="sine" && summarizePicRfPhases(domain,embedded,.25,nullptr).size()==2,"RF results own strings and respect embedded enable");

    using Step=PicStepSeries;using Species=PicSpeciesSeries;
    const std::string dir="explicit_reports/";std::filesystem::create_directory(dir);
    PicBoundaryLosses::Totals losses;losses.domainMacroE={1,2,3,4};losses.unknownE=5;
    losses.domainMacroI={6,7,8,9};losses.unknownI=9007199254740993ULL;
    PicSecondaryEmission::Stats emission;emission.ionImpactMacroCathode=1;emission.ionImpactPhysCathode=2.;
    emission.ionImpactEnergyJCathode=3.;emission.emittedMacro=4;emission.emittedPhys=5.;emission.emittedEnergyJ=6.;
    auto open=[&](Step& step,Species& sp,int resume) {
        step.openSecondaryEmission(dir+"see.txt");step.openImpact(dir+"impact.csv",resume);
        step.openRfPhase(dir+"phase.csv",resume,true);sp.openWeightStatistics(dir+"weights.csv",resume);
    };
    auto write=[&](Step& step,Species& sp,int n) {
        step.writeSecondaryEmission(n,.25,emission,{7.,8.},false);step.writeImpact(n,.25,losses);
        step.writeRfPhase(n,.25,siPhase);sp.writeWeightStatistics(n,.25,"weighted",sample);
        sp.writeWeightStatistics(n,.25,"empty",empty);
    };
    {Step step(Step::Mode::Explicit);Species sp(Species::Mode::Explicit,{});
     check(!step.secondaryEmissionReady() && !step.impactReady() && !step.rfPhaseReady() && !sp.weightStatisticsReady(),"report streams initially unavailable");
     open(step,sp,0);write(step,sp,1);
     check(bytes(dir+"see.txt").empty() && bytes(dir+"phase.csv").empty() && bytes(dir+"weights.csv").empty(),"SEE RF and weights retain buffered writes");
     check(csvRows(dir+"impact.csv").size()==2,"impact row flushes immediately");
     write(step,sp,2);step.flushRfPhase();sp.flushWeightStatistics();
     step.writeSecondaryEmission(3,.25,emission,{7.,8.},true);
     check(bytes(dir+"see.txt").find("3,0.25,1,2,3,4,5,6,7,8")!=std::string::npos,"SEE explicit flush publishes the buffered report");}
    const auto impact=csvRows(dir+"impact.csv"),stats=csvRows(dir+"weights.csv"),rf=csvRows(dir+"phase.csv");
    check(impact.size()==3 && impact[1].size()==12 && impact[1][1]=="0.25"
        && impact[1].back()=="9007199254740993","impact schema precision and integer counts");
    for(int i=2;i<11;++i)check(std::stoi(impact[1][i])==i-1,"impact side and unknown ordering");
    check(stats.size()==5 && stats[0].size()==15 && stats[1][3]=="6" && stats[1][4]=="22"
        && stats[1][7]=="3" && stats[1][8]=="10" && stats[1][14]=="2","weight CSV mapping and schema");
    check(stats[2]==std::vector<std::string>{"1","0.25","empty","0","0","0","0","0","0","0","0","0","0","-1","-1"},"empty weight row exact layout");
    check(rf.size()==7 && rf[0][1]=="time_s" && rf[0].back()=="phi_SI_V" && rf[1].back()=="50","RF SI columns and values");
    check(bytes(dir+"see.txt").find("1,0.25,1,2,3,4,5,6,7,8")!=std::string::npos,"SEE columns and destructor flush");
    {Step step(Step::Mode::Explicit);Species sp(Species::Mode::Explicit,{});open(step,sp,2);write(step,sp,3);
     step.flushRfPhase();sp.flushWeightStatistics();}
    check(csvRows(dir+"impact.csv").back()[0]=="3" && csvRows(dir+"impact.csv").size()==3
        && csvRows(dir+"weights.csv").size()==5 && csvRows(dir+"phase.csv").size()==7,"all report CSV histories prune checkpoint rows");
    check(bytes(dir+"see.txt").find("1,0.25")==std::string::npos
        && bytes(dir+"see.txt").find("3,0.25,1,2,3,4,5,6,7,8")!=std::string::npos,"SEE deliberately truncates on restart");
    {Step step(Step::Mode::Explicit);step.openRfPhase(dir+"code.csv",0,false);step.writeRfPhase(1,2.,phase);step.flushRfPhase();}
    check(csvRows(dir+"code.csv")[0][1]=="time" && csvRows(dir+"code.csv")[0].back()=="phi"
        && csvRows(dir+"code.csv")[1][1]=="2" && csvRows(dir+"code.csv")[1].back()=="5","RF code time and potential remain caller-selected");
    {Step step(Step::Mode::Explicit);Species sp(Species::Mode::Explicit,{});open(step,sp,0);}
    check(csvRows(dir+"impact.csv").size()==1 && csvRows(dir+"weights.csv").size()==1
        && csvRows(dir+"phase.csv").size()==1,"fresh report session replaces history");
    {std::ofstream blocker(dir+"blocked");blocker<<"keep";}
    std::ostringstream warning;auto* savedError=std::cerr.rdbuf(warning.rdbuf());
    Step failed(Step::Mode::Explicit);Species failedSpecies(Species::Mode::Explicit,{});
    failed.openSecondaryEmission(dir+"blocked/see.txt");failed.openImpact(dir+"blocked/impact.csv",0);
    failed.openRfPhase(dir+"blocked/phase.csv",0,true);failedSpecies.openWeightStatistics(dir+"blocked/weights.csv",0);
    std::cerr.rdbuf(savedError);
    check(!failed.secondaryEmissionReady() && !failed.impactReady() && !failed.rfPhaseReady()
        && !failedSpecies.weightStatisticsReady() && bytes(dir+"blocked")=="keep","unavailable independent reports preserve blocker");
    check(warning.str()=="[ExplicitPIC] WARNING: failed to create see report: explicit_reports/blocked/see.txt\n"
        "[ExplicitPIC] WARNING: failed to create impact report: explicit_reports/blocked/impact.csv\n"
        "[ExplicitPIC] WARNING: failed to create RF phase report: explicit_reports/blocked/phase.csv\n"
        "[ExplicitPIC] WARNING: failed to create particle weight stats: explicit_reports/blocked/weights.csv\n","legacy report warnings preserved");
}
void checkExplicitEmbeddedBand() {
    using Metrics=Solvers::Poisson::CutCellMetrics;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=64;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph faces;faces.build(mesh);
    CellField<Vector2D> field(mesh,Vector2D(3.,4.));Metrics metrics;
    auto& cells=const_cast<std::vector<Metrics::CellMetric>&>(metrics.cellMetricsByIndex());
    cells.resize(4);cells[0].isCut=true;
    // Deliberately contradictory alpha: the explicit seed is isCut alone.
    for(auto& c:cells)c.alpha=0.;cells[0].alpha=1.;
    auto band=[&](int layers,double scale=2.) {
        return summarizeExplicitPicEmbeddedBandField(mesh,faces,&metrics,layers,field,scale);
    };
    auto seed=band(-3);
    check(seed.layers==0 && seed.cutCells==1 && seed.cells==1 && seed.meanVpm==10. && seed.rmsVpm==10.,
        "explicit band seeds use isCut and clamp negative layer counts");
    check(band(1).cells==3 && band(2).cells==4 && band(100).cells==4,"explicit band expands one face layer at a time");
    field.value[0]={0.,0.};auto full=band(2);
    check(full.meanVpm==7.5 && full.rmsVpm==std::sqrt(75.),"explicit unweighted cell mean and RMS");
    check(band(2,-2.).meanVpm==-7.5 && band(2,-2.).rmsVpm==-std::sqrt(75.),
        "explicit scale applies after reduction, including RMS sign");
    const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
    field.value[0]={nan,0.};
    check(band(2).cells==4 && std::isnan(band(2).rmsVpm) && std::isnan(band(2).meanVpm),
        "explicit nonfinite selected field propagates without implicit filtering");
    field.value[0]={inf,0.};check(std::isinf(band(0).rmsVpm),"explicit infinite field retained");
    field.value[0]={3.,4.};field.value[3]={nan,0.};check(band(0).rmsVpm==10.,"unselected nonfinite field ignored");
    field.value[3]={3.,4.};
    for(auto& segment:const_cast<std::vector<FaceSegment>&>(faces.segments()))
        if(segment.L>=0)segment.boundary=Boundary::west;
    check(band(100).cells==1,"explicit band does not cross boundary-marked neighbor faces");
    std::vector<double> alpha(4,1.);alpha[0]=.5;
    check(summarizePicEmbeddedBandField(mesh,faces,true,alpha,100,field,2.).cells==4,
        "implicit connectivity policy remains distinct");
    faces.build(mesh);cells[0].isCut=false;
    check(band(2).cutCells==0 && band(2).cells==0 && band(2).rmsVpm==0.,"valid empty explicit band");
    check(std::isnan(band(2,inf).rmsVpm),"valid empty band still receives legacy SI scale");
    auto absent=summarizeExplicitPicEmbeddedBandField(mesh,faces,nullptr,2,field,inf);
    check(absent.cutCells==0 && absent.cells==0 && absent.rmsVpm==0. && absent.meanVpm==0.,
        "unavailable metrics return zeros without applying nonfinite scale");
    cells[0].isCut=true;CellField<Vector2D> missing;
    check(summarizeExplicitPicEmbeddedBandField(mesh,faces,&metrics,2,missing,2.).cells==0,
        "absent field clears counts rather than retaining implicit denominator");
    auto saved=std::move(field.value);field.value=nullptr;
    auto nullField=band(2);field.value=std::move(saved);
    check(nullField.cells==0 && nullField.cutCells==0,"null storage with current field size rejected");
    cells.pop_back();check(band(2).cells==0,"stale metric size rejected before seed lookup");
    DictionaryMesh::RefineBox box;box.min={0.,.5};box.max={1.,1.};box.level=1;
    mesh.applyRefineBoxes({box});faces.build(mesh);
    cells.resize(mesh.cellIndexNum());
    check(band(100).cells==0,"stale field after AMR rejected");
    field.initialize(mesh,Vector2D(3.,4.));
    auto adapted=band(100);
    check(adapted.cutCells==1 && adapted.cells==mesh.cellIndexNum() && adapted.cells>4 && adapted.rmsVpm==10.,
        "explicit band uses current nonuniform AMR graph without retaining prior indices");

    using Series=PicStepSeries;const std::string path="explicit_band.csv";
    const std::string header="step,time_s,band_layers,cut_cells,band_cells,E_rms_Vpm,E_mean_Vpm\n";
    auto text=[&]() {std::ifstream file(path);return std::string(std::istreambuf_iterator<char>(file),{});};
    {Series s(Series::Mode::Explicit);s.openEmbeddedBand(path,0);
     check(s.embeddedBandReady() && bytes(path).empty(),"explicit band header initially buffered");
     s.writeEmbeddedBand(1,1e-11,{2,4,10.,7.5,1});s.flushEmbeddedBand();
     check(text()==header+"1,9.999999999999999e-12,2,1,4,10,7.5\n","explicit schema and defaultfloat precision");
     s.writeEmbeddedBand(2,.25,{0,0,0.,0.,0});s.flushEmbeddedBand();}
    {Series s(Series::Mode::Explicit);s.openEmbeddedBand(path,2);
     s.writeEmbeddedBand(3,.5,{1,3,6.,5.,1});s.flushEmbeddedBand();}
    const auto resumed=csvRows(path);
    check(resumed.size()==3 && resumed[1][0]=="1" && resumed[2][0]=="3","explicit band restart prunes checkpoint row inclusively");
    {Series s(Series::Mode::Explicit);s.openEmbeddedBand(path,0);s.flushEmbeddedBand();}
    check(text()==header,"fresh explicit band truncates old rows");
    {Series s(Series::Mode::Explicit);s.openEmbeddedBand("band_new_parent/report.csv",0);
     check(s.embeddedBandReady(),"shared CSV preparation creates missing parents");}
    {std::ofstream blocker("band_blocked_parent");blocker<<"keep";}
    std::ostringstream warning;auto* savedError=std::cerr.rdbuf(warning.rdbuf());
    {Series s(Series::Mode::Explicit);s.openEmbeddedBand("band_blocked_parent/report.csv",0);
     std::cerr.rdbuf(savedError);check(!s.embeddedBandReady(),"failed explicit band stream unavailable");}
    check(warning.str()=="[ExplicitPIC] WARNING: failed to create EB-band RMS csv: band_blocked_parent/report.csv\n"
        && bytes("band_blocked_parent")=="keep","explicit directory error is nonfatal and warns");
}
void checkImplicitDiagnosticSeries() {
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=64;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph faces;faces.build(mesh);
    CellField<double> amount(mesh,2.);const int n=mesh.cellIndexNum();check(n==4,"amount fixture");
    auto count=summarizePicSpeciesAmount(&amount,n,4.);
    check(count.stage==8. && count.physical==4. && count.relative==1.,"stage amounts are counts, without volume factors");
    check(summarizePicSpeciesAmount(nullptr,n,3.).relative==-1.
        && summarizePicSpeciesAmount(&amount,n+1,3.).stage==0.,"absent or stale-size amount field");
    check(summarizePicSpeciesAmount(nullptr,n,0.).relative==0.
        && std::isinf(summarizePicSpeciesAmount(&amount,n,0.).relative),"zero physical count branches");
    check(summarizePicSpeciesAmount(&amount,n,-4.).relative==-3.,"signed physical count retained");
    const double nan=std::numeric_limits<double>::quiet_NaN();
    amount.value[0]=nan;check(std::isnan(summarizePicSpeciesAmount(&amount,n,4.).relative),"nonfinite amount propagates");
    amount.value[0]=2.;
    CellField<Vector2D> field(mesh,Vector2D(3.,4.));std::vector<double> alpha(n,1.);alpha[0]=.5;
    auto band=[&](int layers){return summarizePicEmbeddedBandField(mesh,faces,true,alpha,layers,field,2.);};
    check(band(-3).layers==0 && band(0).cells==1 && band(0).rmsVpm==10. && band(0).meanVpm==10.,"band seed and SI conversion");
    check(band(1).cells==3 && band(2).cells==4 && band(9).cells==4,"face-neighbor layers saturate without duplicates");
    field.value[0]={nan,0.};auto invalid=band(2);
    check(invalid.cells==4 && invalid.meanVpm==7.5 && std::abs(invalid.rmsVpm-std::sqrt(75.))<1e-12,
        "nonfinite field is skipped in numerator but retained in legacy denominator");
    CellField<Vector2D> missing;
    auto empty=summarizePicEmbeddedBandField(mesh,faces,true,alpha,2,missing,2.);
    check(empty.cells==4 && empty.rmsVpm==0. && empty.meanVpm==0.,"missing field retains selected band size");
    check(summarizePicEmbeddedBandField(mesh,faces,false,alpha,2,field,2.).cells==0
        && summarizePicEmbeddedBandField(mesh,faces,true,{},2,field,2.).cells==0,"disabled metrics or wrong alpha size");
    alpha[0]=nan;check(band(2).cells==0,"nonfinite alpha is not a seed");
    alpha[0]=1.-5e-13;check(band(2).cells==0,"cut threshold tolerance");
    alpha[0]=0.;check(band(0).cells==1,"fully blocked seed retained");
    DictionaryMesh::RefineBox box;box.min={0.,.5};box.max={1.,1.};box.level=1;
    mesh.applyRefineBoxes({box});faces.build(mesh);field.initialize(mesh,Vector2D(3.,4.));
    alpha.assign(mesh.cellIndexNum(),1.);alpha[0]=.5;
    auto adapted=band(100);check(adapted.cells==mesh.cellIndexNum() && adapted.cells>4 && adapted.rmsVpm==10.,
        "report uses current AMR graph and fields without retaining prior indices");

    using Series=PicStepSeries;const std::string dir="implicit_diagnostics/";
    const std::vector<std::string> files{"physics.csv","amount.csv","species.csv","nested/band.csv"};
    auto open=[&](Series& s,int resume){s.openPhysics(dir+files[0],resume);s.openAmount(dir+files[1],resume);
        s.openAmountSpecies(dir+files[2],resume);s.openEmbeddedBand(dir+files[3],resume);};
    PicStepParticleBalance::Result charge{};
    charge.stageRelative=3.;charge.relative=4.;charge.sourceC=5.;charge.domainC=6.;charge.embeddedC=7.;
    charge.macroBeginC=8.;charge.macroEndC=9.;charge.stageGridC=11.;charge.stageMacroC=12.;charge.stageResidualC=-1.;
    const long long large=9007199254740993LL;
    auto write=[&](Series& s,int step){
        s.writePhysics(step,.25,{1.,2.,{10.,11.,12.,13.,14.,15.,16.,nan,18.,large,20.,999.}},charge);
        s.writeAmount(step,.25,charge);std::string name="rawSpecies";
        s.writeAmountSpecies(step,.25,name,{8.,4.,1.});name="changed";
        s.writeAmountSpecies(step,.25,"empty",{0.,0.,std::numeric_limits<double>::infinity()});
        s.writeEmbeddedBand(step,.25,{2,4,10.,10.});};
    {
        Series s(Series::Mode::Implicit);
        check(!s.physicsReady() && !s.amountReady() && !s.amountSpeciesReady() && !s.embeddedBandReady(),"unopened diagnostics are not ready");
        open(s,0);for(const auto& f:files)check(csvRows(dir+f).size()==1,"diagnostic headers immediately visible");
        write(s,1);write(s,2);
        for(const auto& f:files)check(csvRows(dir+f).size()==1,"small diagnostic rows stay buffered");
        s.flushPhysics();s.flushAmounts();s.flushEmbeddedBand();
    }
    const auto physics=csvRows(dir+files[0]),amountRows=csvRows(dir+files[1]),speciesRows=csvRows(dir+files[2]);
    check(physics.size()==3 && physics[0].size()==22 && physics[1].size()==22
        && physics[1][1]=="0.2500000000000000" && physics[1][20]==std::to_string(large),"physics schema, precision and integer macro counts");
    for(int i=2;i<20;++i)check(i==18?std::isnan(std::stod(physics[1][i])):std::stod(physics[1][i])==i-1.,"physics columns and charge ordering");
    check(std::stod(physics[1][21])==20. && physics[0][4]=="dQrel_stage" && physics[0][5]=="dQphys_rel","distinct stage/physical balance");
    check(amountRows[0].size()==6 && std::stod(amountRows[1][2])==11. && std::stod(amountRows[1][4])==-1.,"amount schema and values");
    check(speciesRows.size()==5 && speciesRows[1][2]=="rawSpecies" && std::isinf(std::stod(speciesRows[2][5])),"species views are synchronous, infinity preserved");
    check(csvRows(dir+files[3])[1]==std::vector<std::string>{"1","0.2500000000000000","2","4","10.0000000000000000","10.0000000000000000"},"band CSV exact row");
    {Series s(Series::Mode::Implicit);open(s,2);write(s,3);}
    for(const auto& f:files){const auto r=csvRows(dir+f);check(r[1][0]=="1" && r.back()[0]=="3"
        && r.size()==(f==files[2]?5u:3u),"restart prunes all species rows at/after resume and destruction flushes");}
    {Series s(Series::Mode::Implicit);open(s,0);}
    for(const auto& f:files)check(csvRows(dir+f).size()==1,"fresh session replaces prior diagnostic history");
    {std::ofstream blocker(dir+"blocked");blocker<<"keep";}
    Series failed(Series::Mode::Implicit);failed.openPhysics(dir+"blocked/physics.csv",0);
    failed.openAmount(dir+"blocked/amount.csv",0);failed.openAmountSpecies(dir+"blocked/species.csv",0);
    check(!failed.physicsReady() && !failed.amountReady() && !failed.amountSpeciesReady(),"failed streams are independent and unavailable");
    std::filesystem::create_directories(dir+"band_blocked/occupied");
    failed.openEmbeddedBand(dir+"band_blocked",0);check(!failed.embeddedBandReady(),"band open failure retains silent unready behavior");
    check(bytes(dir+"blocked")=="keep","failed open preserves unrelated file");
}
void checkExplicitConsole() {
    struct Console {ConsoleControl::Config saved=ConsoleControl::g_config;~Console(){ConsoleControl::g_config=saved;}} restore;
    ConsoleControl::setConfig("info",99);
    DictionaryMesh settings;settings.Nx=settings.Ny=2;settings.maxCellNum=32;settings.maxInitialIterations=0;
    QuadtreeMesh mesh(settings);CellField<double> phi(mesh,0.);phi.value[0]=-2.;phi.value[1]=4.;
    UnitsScales scales;scales.L0_m=2.;scales.T0_s=3.;scales.Phi0_V=10.;scales.v0_mps=2.;
    const auto units=UnitsContext::MakeSI(scales);
    std::vector<SpeciesConfig> configs(3);
    configs[0].name="electron";configs[0].m=2.*PhysicalConstants::eV_J;configs[0].q=-1.;configs[0].weight=3.;
    configs[1].name="ion";configs[1].q=1.;configs[2].name="neutral";
    PicSpeciesState particles(3,123);particles.initialize(configs,NormalizationScales{},mesh,true);
    Particle a;a.velocity={1.,2.};a.vz=3.;a.w=2.;a.rngState=123;
    Particle b;b.vz=1.;b.w=0.;b.rngState=456;
    particles.species()[0].particle={a,b};particles.species()[1].particle.resize(1);
    MCCCollisionOperator::StepStats collisions;collisions.n_attempts=50;collisions.n_collisions=41;
    collisions.n_null=9;collisions.n_spawned=4;collisions.max_nu_dt_majorant=.25;
    collisions.per_channel={{"e::ELASTIC",2},{"e::Excitation_1",3},{"e::IONIZATION",5},{"e::Attachment",7},{"i::charge-exchange",11},{"custom",13}};
    auto& e=collisions.per_projectile["electron"];e.n_attempts=17;e.n_collisions=8;e.n_null=9;e.n_spawned=2;
    e.max_nu_dt_majorant=.1;e.per_channel={{"Elastic",2},{"excitation_1",3}};
    auto& i=collisions.per_projectile["ion"];i.n_attempts=11;i.n_collisions=11;i.per_channel={{"CX",11}};
    auto& n=collisions.per_projectile["unknown"];n.n_attempts=13;n.n_collisions=13;n.per_channel={{"custom",13}};
    collisions.spawn_per_species={{"z",2},{"a",4}};
    PicBoundaryLosses::Totals losses;losses.domainMacroE={1,2,3,4};losses.domainMacroI={4,3,2,1};
    losses.domainChargeE={24.,48.,36.,48.};losses.domainChargeI={2.4,4.8,3.6,4.8};
    losses.unknownE=999;losses.periodicLeakE={999,999,999,999};
    losses.embeddedCounts.inside=2;losses.embeddedCounts.crossing=3;losses.embeddedCounts.unknown=4;losses.embeddedCounts.byId={{7,5}};
    PicExplicitConsoleInput physical{mesh,phi,&units,particles,collisions,losses,-8.,-10.,2.,4.,.5,true};
    auto values=summarizeExplicitPicConsole(physical,6.);
    const auto close=[](double x,double y){return std::abs(x-y)<1e-11;};
    check(values.state.phiMinV==-20. && values.state.phiMaxV==40. && close(values.state.meanElectronEV,24.8)
        && close(values.state.temperatureEV,24.8*2./3.) && close(values.state.electronTailFraction,.4),"explicit console SI potential and three-dimensional electron energy");
    check(values.macroCount==3 && values.species[0].macro==2 && values.species[0].physical==5.
        && values.species[1].physical==1. && values.species[2].macro==0 && values.wallLossMacro==20 && values.embeddedAbsorbed==5,
        "explicit populations and wall counts exclude unknown/periodic leaks and use particle weight fallback");
    check(values.chargeDifference==2. && values.chargeRelative==.2 && close(values.totalCurrent,11.)
        && values.electronCurrent==std::array<double,4>{1.,2.,3.,4.},"explicit deposited charge and physical side-current areas");
    check(values.state.elastic==2 && values.state.excitation==3 && values.state.ionization==5 && values.state.attachment==7
        && values.state.chargeExchange==11 && values.state.other==13 && values.state.electrons.att==17
        && values.state.ions.cx==11 && values.state.neutrals.other==13,"explicit collision names and configured projectile split");
    auto saved=particles.species()[0];particles.species()[0]=ParticleType("electron",configs[0].m,1.);
    check(summarizeExplicitPicConsole(physical,6.).state.electrons.att==17
        && summarizePicConsole(mesh,phi,&units,particles.species(),particles.weights(),collisions).ions.att==28,
        "explicit keeps configured projectile classes while implicit uses live charge");
    particles.species()[0]=std::move(saved);
    const std::unordered_set<std::string> overlap{"electron"};const PicConsoleValues::ProjectileGroups groups{overlap,overlap};
    check(summarizePicConsole(mesh,phi,&units,particles.species(),particles.weights(),collisions,&groups).electrons.att==17,
        "explicit configured electron names take precedence over overlapping ion names");
    auto noUnits=physical;noUnits.units=nullptr;noUnits.particleCharge=0.;noUnits.embeddedStatistics=false;
    values=summarizeExplicitPicConsole(noUnits,0.);
    check(values.state.phiMinV==-2. && values.state.phiMaxV==4. && close(values.state.meanElectronEV,6.2)
        && values.totalCurrent==0. && values.chargeRelative==0. && values.embeddedAbsorbed==0,"no-unit and invalid-time fallbacks");
    values=summarizeExplicitPicConsole(noUnits,2.);check(close(values.electronCurrent[0],3.) && close(values.electronCurrent[2],9.),"without units side areas ignore SI slab");
    auto zeroArea=physical;zeroArea.slabMeters=0.;check(summarizeExplicitPicConsole(zeroArea,6.).totalCurrent==0.,"nonpositive explicit area gate");

    PhysicsDiagnostics::ElectrodeCurrentGroups electrodes({{"driven",{1}}});electrodes.record(1,-12.,0.);
    ExplicitPIC::Diagnostics::PicDiagnostics::Result diagnostic;diagnostic.macro_total=999;diagnostic.rG_Linf_rel_rhs=0.;
    PicExplicitConsole::Sample sample{physical,electrodes,&diagnostic,{7,8,9.,10.,11.,12.},{6,1.,2.,3.,true,true,true,true,true,13}};
    PicStepSeries::ExplicitPerformance timing{};timing.totalSeconds=.25;timing.field.buildWeights_s=.1;
    timing.poissonSeconds[3]=.05;timing.amrAdaptSeconds=.025;
    const auto render=[&](const PicExplicitConsole& report,int step,double eta=3661.){
        std::ostringstream stream;stream<<std::fixed<<std::setprecision(9);
        report.writeStep(step,2.*step,2.,eta,timing,sample,stream);
        check(stream.precision()==9 && (stream.flags()&std::ios::floatfield)==std::ios::fixed,"explicit report preserves caller formatting");
        return stream.str();
    };
    PicConfig cfg;cfg.output.console.every=2;cfg.output.console.physicsDetail="HIGH";
    PicExplicitConsole high(cfg,0,9);cfg.output.console.every=100;cfg.output.console.physicsDetail="minimal";
    check(render(high,1).empty() && render(high,9).empty(),"explicit cadence snapshot does not force a non-cadence final report");
    auto line=render(high,2);
    check(line.find("Progress step=2/9 t_s=1.200e+01 dt_s=6.000e+00 np=3 stage=startup evt=amr-adapt,output ETA=01:01:01")!=std::string::npos
        && line.find("Particles electron=2(5.000e+00) ion=1(1.000e+00) neutral=0(0.000e+00)")!=std::string::npos,
        "explicit current population is independent of old sampled diagnostic counts");
    check(line.find("phi[V]={[-20,40]}")!=std::string::npos && line.find("eMeanE_eV=24.800 Te_eV=16.533 tail_e_gt_11.590eV=0.400")!=std::string::npos
        && line.find("Igrp_A{driven=-2.000e+00(e=-2.000e+00 i=0.000e+00)} Jtot_Apm2=1.100e+01")!=std::string::npos,"explicit state/circuit formats and units");
    check(line.find("N[s/m/d/e]=6/4/20/5")!=std::string::npos && line.find("MCC att=50 coll=41 null=9 spawn=4 byChan(el/ex/ion/att/cx/other)=2/3/5/7/11/13")!=std::string::npos
        && line.find("SEE ionImp_cath=7(9.000e+00) emit=8(1.000e+01) eEmitMean_eV=11.000 eEmitStd_eV=12.000")!=std::string::npos,
        "explicit source/sink, collision and SEE sections");
    check(line.find("Embedded absorb=5 inside=2 crossing=3 unk=4 byId(7:5) residual=13")!=std::string::npos
        && line.find("MCCSplit e(att/coll/null/spawn)=17/8/9/2")!=std::string::npos
        && line.find("spawnBySpecies=a:4,z:2")!=std::string::npos
        && line.find("Perf total_s=2.500e-01 top=buildWeights:1.000e-01|poisson:5.000e-02|amr:2.500e-02")!=std::string::npos,
        "explicit embedded order, projectile split, spawn ranking and performance labels");
    sample.diagnostic=nullptr;line=render(high,2,-1.);check(line.find("diag=pending")!=std::string::npos && line.find("ETA=N/A")!=std::string::npos,"missing shared diagnostics and invalid ETA");
    cfg.output.console.every=1;cfg.output.v2.console.format="MULTILINE";PicExplicitConsole minimal(cfg,10,30);
    sample.events.adapted=sample.events.outputWritten=false;line=render(minimal,11);
    check(line.find("[ExplicitPIC][Progress] step=11/30")!=std::string::npos && line.find("stage=field-init evt=advance")!=std::string::npos
        && line.find("[ExplicitPIC][Particles]")==std::string::npos && line.find("[ExplicitPIC][SEE]")==std::string::npos,"minimal multiline and resumed phase");
    sample.diagnostic=&diagnostic;sample.physical.gridCharge=-10.;line=render(minimal,21);
    check(line.find("stage=quasi-steady")!=std::string::npos,"late explicit phase reads shared sampled diagnostics");
    diagnostic.rG_Linf_rel_rhs=.1;check(render(minimal,21).find("stage=transient")!=std::string::npos,"changed shared diagnostic affects next report only");
    cfg.output.console.physicsDetail="unknown";cfg.output.v2.console.format="other";PicExplicitConsole medium(cfg,0,30);
    sample.events.collisionsEnabled=false;sample.events.seeEnabled=false;sample.emission={};
    line=render(medium,2);check(line.find("Particles ")!=std::string::npos && line.find("mcc=off")!=std::string::npos && line.find("Perf ")==std::string::npos,"unknown detail/layout fallback and collisions off");
    cfg.output.console.physicsShowCharge=cfg.output.console.physicsShowCurrent=cfg.output.console.physicsShowWall=cfg.output.console.physicsShowMCC=cfg.output.console.physicsShowPerf=false;
    PicExplicitConsole flags(cfg,0,30);line=render(flags,2);
    check(line.find("phi[V]")==std::string::npos && line.find("Circuit ")==std::string::npos && line.find("Wall ")==std::string::npos
        && line.find("mcc=")==std::string::npos && line.find("State eMeanE_eV=")!=std::string::npos,"explicit display switches retain mandatory sections");
    PhysicsDiagnostics::ElectrodeCurrentGroups ungrouped;
    PicExplicitConsole::Sample sides{sample.physical,ungrouped,nullptr,{},sample.events};
    std::ostringstream sideOutput;high.writeStep(2,4.,2.,0.,timing,sides,sideOutput);
    check(sideOutput.str().find("Je_Apm2{L=1.000e+00 R=2.000e+00 B=3.000e+00 T=4.000e+00}")!=std::string::npos,
        "ungrouped explicit circuit retains physical side currents");
    cfg.output.console.physicsEnabled=false;PicExplicitConsole disabled(cfg,0,30);
    check(render(disabled,2).empty(),"disabled summaries skip report");ConsoleControl::setConfig("warn",1);
    check(render(high,2).empty(),"explicit info gate remains live");
    check(particles.species()[0].particle[0].rngState==123 && particles.species()[0].particle[1].vz==1. && phi.value[0]==-2.,"console leaves physical state untouched");
}
void checkImplicitConsole() {
    static_assert(!std::is_copy_constructible_v<PicImplicitConsole> && !std::is_move_constructible_v<PicImplicitConsole>);
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);CellField<double> phi(mesh,0.);phi.value[0]=-2.;phi.value[1]=4.;
    UnitsScales scales;scales.Phi0_V=10.;scales.v0_mps=2.;const auto units=UnitsContext::MakeSI(scales);
    std::vector<ParticleType> species;species.emplace_back("electron",2.*PhysicalConstants::eV_J,-1.);
    species.emplace_back("ion",1.,1.);species.emplace_back("neutral",1.,0.);
    Particle a;a.velocity={1.,2.};a.vz=3.;a.w=2.;a.rngState=123;
    Particle b;b.velocity={0.,0.};b.vz=1.;b.w=0.;b.rngState=456;
    species[0].particle={a,b};species[1].particle.resize(1);
    std::vector<double> weights{3.,1.,1.};MCCCollisionOperator::StepStats collisions;
    collisions.per_channel={{"e::ELASTIC",2},{"e::Excitation_1",3},{"e::IONIZATION",5},{"e::Attachment",7},{"i::charge-exchange",11},{"custom",13}};
    collisions.n_attempts=50;collisions.n_collisions=41;collisions.n_null=9;collisions.n_spawned=4;collisions.max_nu_dt_majorant=.25;
    auto& electron=collisions.per_projectile["electron"];electron.n_attempts=17;electron.n_collisions=8;
    electron.n_null=9;electron.n_spawned=2;electron.max_nu_dt_majorant=.1;electron.per_channel={{"Elastic",2},{"excitation_1",3}};
    auto& ion=collisions.per_projectile["ion"];ion.n_attempts=11;ion.n_collisions=11;ion.per_channel={{"CX",11}};
    auto& unknown=collisions.per_projectile["not-a-species"];unknown.n_attempts=13;unknown.n_collisions=13;unknown.per_channel={{"custom",13}};
    collisions.spawn_per_species={{"z",2},{"a",4}};
    const auto values=summarizePicConsole(mesh,phi,&units,species,weights,collisions);
    const auto close=[](double x,double y){return std::abs(x-y)<1e-11;};
    check(values.phiMinV==-20. && values.phiMaxV==40. && close(values.meanElectronEV,24.8)
        && close(values.temperatureEV,24.8*2./3.) && close(values.electronTailFraction,.4),
        "implicit console uses SI potential, full 3D electron energy, weighted tail and default weight");
    check(values.elastic==2 && values.excitation==3 && values.ionization==5 && values.attachment==7
        && values.chargeExchange==11 && values.other==13 && values.electrons.att==17
        && values.ions.cx==11 && values.neutrals.other==13,"console preserves channel aliases and unknown-projectile grouping");
    check(species[0].particle[0].rngState==123 && species[0].particle[1].rngState==456
        && species[0].particle[0].velocity.y==2. && phi.value[0]==-2.,"console reductions leave physical inputs unchanged");

    PicStepParticleBalance::Result charge;charge.stageRelative=charge.relative=0.;
    PicBoundaryLosses::Totals losses;losses.domainMacroE={1,2,3,4};losses.domainMacroI={4,3,2,1};
    losses.unknownE=5;losses.unknownI=6;losses.embeddedMacro={5,6};losses.embeddedEnergyJ=7.;
    losses.embeddedCounts.byId={{9,3},{2,3},{5,7}};
    std::array<double,4> je{1.,2.,3.,4.},ji{.1,.2,.3,.4},qe{100.,101.,102.,103.},qi{10.,11.,12.,13.};
    PhysicsDiagnostics::ElectrodeCurrentGroups electrodes({{"driven",{1}}});electrodes.record(1,-4.,0.);
    PicElectrodeCurrent::Sample current;current.ready=true;current.conduction=1.23;current.displacement=9.87;
    PicImplicitConsole::Sample sample{mesh,phi,&units,species,weights,charge,losses,{je,ji,qe,qi},electrodes,current,collisions,
        {false,false,false,true,true,8,9}};
    PicImplicitStepTiming timing;timing.total_s=.25;timing.pm_s=.125;
    auto timingSummary=summarizePicImplicitTiming(timing);
    ExplicitPIC::Diagnostics::PicDiagnostics::Result diagnostic;
    diagnostic.speciesNames={"electron","ion","truncated"};diagnostic.speciesMacroCount={7,2};diagnostic.speciesPhysCount={5.,1.};
    const auto render=[&](PicImplicitConsole& report,int step,bool defaultSummary=false) {
        std::ostringstream stream;
        const auto count=report.writeStep(step,(step+1)*2.,2.,defaultSummary,timing,timingSummary,sample,stream);
        check(count==3,"console returns actual current population independently of sampled diagnostic counts");return stream.str();
    };
    PicConfig cfg;cfg.output.console.every=2;cfg.output.console.physicsDetail="HIGH";
    ImplicitPicOverrides controls;PicImplicitConsole high(cfg,controls,0,9,3,true);high.observeDiagnostics(diagnostic);
    diagnostic.speciesNames[0]="mutated";cfg.output.console.physicsDetail="minimal";cfg.output.console.every=100;
    check(render(high,0).empty(),"physics summary respects captured cadence");
    auto line=render(high,1);
    check(line.find("Progress step=2/9")!=std::string::npos && line.find("np=3")!=std::string::npos
        && line.find("Particles electron=7(5.000e+00) ion=2(1.000e+00)")!=std::string::npos
        && line.find("mutated")==std::string::npos && line.find("truncated")==std::string::npos,
        "console owns diagnostic snapshot and clamps mismatched species vectors");
    check(line.find("eMeanE_eV=2.480e+01")!=std::string::npos && line.find("Jtot_Apm2=1.100e+01")!=std::string::npos
        && line.find("driven=-2.000e+00(e=-2.000e+00 i=0.000e+00)")!=std::string::npos
        && line.find("Jcond_Apm2=1.230e+00")!=std::string::npos,"console circuit preserves grouped charge and nonfinite total-current fallback");
    check(line.find("N[s/m/d/e]=8/9/20/11")!=std::string::npos && line.find("qout_e[C]=1.00e+02")!=std::string::npos
        && line.find("byId(5:7,2:3,9:3)")!=std::string::npos && line.find("MCCSplit e(att/coll/null/spawn)=17/8/9/2")!=std::string::npos
        && line.find("Perf total_s=2.500e-01")!=std::string::npos,"console keeps separate implicit side charges, sorted boundary counts, collision splits and performance");
    check(render(high,8).find("step=9/9")!=std::string::npos,"last step always reports when physics summaries enabled");
    cfg.output.v2.console.format="MULTILINE";cfg.output.console.every=1;
    PicImplicitConsole minimal(cfg,controls,10,30,3,true);minimal.observeDiagnostics(diagnostic);
    line=render(minimal,10);
    check(line.find("[ImplicitPIC][Progress] step=11/30")!=std::string::npos && line.find("stage=field-init")!=std::string::npos
        && line.find("[ImplicitPIC][Particles]")==std::string::npos && line.find("[ImplicitPIC][Wall]")==std::string::npos
        && line.find("[ImplicitPIC][MCC]")==std::string::npos,"minimal multiline output and restored first-step phase");
    check(render(minimal,20).find("stage=quasi-steady")!=std::string::npos,"console classifies quiet late steps");
    sample.events.amrAdapted=true;
    check(render(minimal,21).find("stage=transient evt=amr-adapt")!=std::string::npos,"AMR prevents quasi-steady report");
    sample.events.amrAdapted=false;
    cfg.output.console.physicsEnabled=false;
    PicImplicitConsole fallbackCadence(cfg,controls,0,9,3,true);fallbackCadence.observeDiagnostics(diagnostic);
    check(render(fallbackCadence,0,false).empty() && !render(fallbackCadence,0,true).empty(),
        "disabled physics cadence retains caller's default summary gate");
    PicImplicitConsole silent(cfg,controls,0,9,3,false);silent.observeDiagnostics(diagnostic);
    check(render(silent,0,true).empty(),"info-disabled console remains silent while counting population");

    controls.consoleDqrelWarn=.1;controls.consoleDqrelRecover=.05;controls.consoleDqrelDebounceSteps=2;
    controls.consoleDqrelEventCooldownSteps=5;controls.consoleFallbackFracWarn=0.;
    controls.consoleRgL2Warn=1.;controls.consoleAmrSkipStreakWarn=2;
    PicImplicitConsole warnings(cfg,controls,0,50,3,true);warnings.observeDiagnostics(diagnostic);
    charge.stageRelative=.2;check(render(warnings,0).empty(),"charge warning debounces first exceedance");
    check(render(warnings,1).find("[WARN] |dQrel_stage|=2.000e-01")!=std::string::npos,"charge warning after consecutive exceedances");
    charge.stageRelative=.01;check(render(warnings,2).empty() && render(warnings,3).empty(),"cooldown suppresses recovery output while transitioning state");
    charge.stageRelative=.2;check(render(warnings,4).empty() && render(warnings,5).empty(),"cooldown suppresses repeated warning output");
    charge.stageRelative=.01;check(render(warnings,6).empty() && render(warnings,7).find("[RECOVER] dQrel recovered at step=8")!=std::string::npos,
        "suppressed transitions retain original cooldown timestamp semantics");
    diagnostic.rG_L2_per_active=2.;warnings.observeDiagnostics(diagnostic);
    check(render(warnings,8).find("[WARN] rG_L2_per_active=2.000e+00")!=std::string::npos,"Gauss residual warns without charge debounce");
    diagnostic.rG_L2_per_active=std::numeric_limits<double>::quiet_NaN();warnings.observeDiagnostics(diagnostic);
    check(render(warnings,9).find("[RECOVER] rG_L2_per_active recovered")!=std::string::npos,"legacy nonfinite Gauss value ends warning");
    sample.events.amrIndicator=sample.events.amrSkippedByGate=true;
    check(render(warnings,10).empty(),"first skipped indicator does not warn");
    sample.events.amrIndicator=false;check(render(warnings,11).empty(),"non-indicator step leaves skip streak intact");
    sample.events.amrIndicator=true;
    check(render(warnings,12).find("AMR gate skipped for 2 consecutive indicator steps at step=13")!=std::string::npos,
        "AMR streak counts only indicator steps");
    check(render(warnings,13).empty(),"AMR streak warns once at exact threshold");
    sample.events.amrSkippedByGate=false;render(warnings,14);
    sample.events.amrIndicator=false;

    controls.consoleFallbackFracWarn=.1;controls.consoleFallbackDebounceSteps=2;
    PicImplicitConsole noDiagnostic(cfg,controls,0,50,3,true);
    timingSummary.totalFallbackFraction=.2;
    check(render(noDiagnostic,0).empty() && render(noDiagnostic,1).find("[WARN] coupling_fallback_frac=2.000e-01")!=std::string::npos,
        "coupling warning updates even without a sampled physical diagnostic");
    timingSummary.totalFallbackFraction=std::numeric_limits<double>::quiet_NaN();
    check(render(noDiagnostic,2).empty() && render(noDiagnostic,3).find("[RECOVER] coupling_fallback_frac recovered")!=std::string::npos,
        "legacy nonfinite fallback fraction counts toward debounced recovery");
    controls.consoleDqrelEventCooldownSteps=0;controls.consoleDqrelDebounceSteps=1;
    PicImplicitConsole legacyCharge(cfg,controls,0,50,3,true);
    diagnostic.dQrel=.3;diagnostic.rG_L2_per_active=0.;legacyCharge.observeDiagnostics(diagnostic);
    charge.stageRelative=std::numeric_limits<double>::quiet_NaN();timingSummary.totalFallbackFraction=0.;
    check(render(legacyCharge,0).find("[WARN] |dQrel|=3.000e-01")!=std::string::npos,"missing stage closure falls back to stored diagnostic charge");
    charge.stageRelative=0.;cfg.output.console.physicsEnabled=true;cfg.output.console.every=1;
    PicImplicitConsole eta(cfg,{},0,30,1,true);diagnostic.dQrel=0.;eta.observeDiagnostics(diagnostic);
    for(int step=0;step<21;++step){timing.total_s=step+1;line=render(eta,step);}
    check(line.find("ETA=00:01:44")!=std::string::npos,"ETA retains only the last twenty completed step costs and original rounding");
    check(species[0].particle[0].rngState==123 && losses.embeddedCounts.byId.size()==3
        && collisions.n_collisions==41,"reporting does not mutate borrowed physical inputs");
}
void checkPerformanceSeries() {
    using Mode=PicStepSeries::Mode;
    PicStepSeries::ExplicitPerformance initial{{1.,2.,3.,4.,0.},{5.,6.,7.,8.},
        0.,0.,9.,10.,11.,23,24,25,true};
    {
        PicStepSeries series(Mode::Explicit);
        series.openPerformance("perf_explicit.csv",0);
        check(series.performanceReady(),"explicit performance stream missing");
        series.writePerformance(0,0.,initial);
        series.writePerformance(1,.125,initial);
        const auto rows=csvRows("perf_explicit.csv"); // explicit writes flush
        check(rows.size()==3 && rows[0].size()==20 && rows[1].size()==20,
            "explicit initialization must have the same 20-column layout as steps");
        check(rows[0][1]=="time_code" && std::stod(rows[2][1])==.125,"explicit code time metadata");
        check(rows[0][10]=="amr_indicator_s" && std::stod(rows[1][10])==0. &&
            std::stod(rows[1][11])==0. && std::stod(rows[1][12])==0. &&
            std::stod(rows[1][13])==9. && std::stod(rows[1][15])==11. &&
            rows[1][16]=="23" && rows[1][19]=="1","initial performance columns shifted");
    }
    const auto join=[](const std::vector<std::string>& row) {
        std::string line;for(const auto& field:row){if(!line.empty())line+=',';line+=field;}return line+'\n';
    };
    const auto current=csvRows("perf_explicit.csv");
    auto legacy=current;legacy[0][1]="time_s";
    legacy[1].erase(legacy[1].begin()+10,legacy[1].begin()+13);
    {std::ofstream file("perf_legacy.csv");for(const auto& row:legacy)file<<join(row);}
    {
        PicStepSeries series(Mode::Explicit);
        series.openPerformance("perf_legacy.csv",1);
        series.writePerformance(1,.25,initial);
    }
    const auto migrated=csvRows("perf_legacy.csv");
    check(migrated.size()==3 && migrated[0]==current[0] && migrated[1].size()==20 &&
        migrated[1][10]=="0" && migrated[1][13]==current[1][13] &&
        migrated[1][16]=="23" && std::stod(migrated[2][1])==.25,
        "legacy migration must preserve tokens and trim future rows before appending");
    const auto once=bytes("perf_legacy.csv");
    {PicStepSeries series(Mode::Explicit);series.openPerformance("perf_legacy.csv",2);}
    check(bytes("perf_legacy.csv")==once,"reopening migrated history must be idempotent");
    // Even malformed rows beyond the resume step are rejected before truncation.
    for(const auto& bad:std::vector<std::string>{
        "step,time_bogus\n0,1\n",join(legacy[0])+join(legacy[1])+"9,0\n",
        join(current[0])+join(legacy[1]),join(legacy[0])+"99999999999999999999,0\n"}) {
        {std::ofstream file("perf_bad.csv");file<<bad;}
        const auto before=bytes("perf_bad.csv");bool rejected=false;
        try {PicStepSeries series(Mode::Explicit);series.openPerformance("perf_bad.csv",1);}
        catch(const std::runtime_error&){rejected=true;}
        check(rejected && bytes("perf_bad.csv")==before,"incompatible performance history changed");
    }
    PicImplicitStepTiming timing;
    timing.mcc_pre_s=3.;timing.amr_s=5.;timing.packedFieldOutput_s=4.;
    timing.buildWeightsN_cwFallbackCount=2;timing.buildWeightsN_particleCount=8;
    timing.orbitAvg_sampleCwFallbackCount=3;timing.orbitAvg_sampleCount=12;
    timing.predictor_segments_total=24;timing.predictor_coarsen_hits=4;
    timing.predictor_coarsen_parent_climb=12;timing.orbitAvg_sampleCache_nnz=48;
    timing.orbitAvg_adaptAcceptCount=8;timing.orbitAvg_adaptRefineCount=4;
    timing.orbitAvg_adaptHitMaxDepthCount=2;timing.orbitAvg_adaptAcceptedDepthSum=24;
    const auto summary=summarizePicImplicitTiming(timing);
    check(summary.buildFallbackFraction==.25 && summary.orbitFallbackFraction==.25 &&
        summary.totalFallbackFraction==.25 && summary.totalFallbackCount==5,
        "fallback summary counts/fractions");
    check(std::string(summary.top[0].name)=="amr" && summary.top[0].t==5. &&
        std::string(summary.top[1].name)=="output" && std::string(summary.top[2].name)=="mcc_pre",
        "performance ranking must cover the original eight stages");
    const PicImplicitStepTiming zero;
    const auto zeroSummary=summarizePicImplicitTiming(zero);
    check(zeroSummary.buildFallbackFraction==0. && zeroSummary.orbitFallbackFraction==0. &&
        zeroSummary.totalFallbackFraction==0.,"empty timing fallback fractions must be zero");
    {
        PicStepSeries series(Mode::Implicit);
        series.openPerformance("perf_implicit.csv",0);
        series.openCouplingFallback("fallback.csv",0);
        check(csvRows("perf_implicit.csv").size()==1 && csvRows("fallback.csv").size()==1,
            "implicit performance headers must flush at open");
        series.writePerformance(1,.25,timing,summary,{4,true,false,true,8,32,"cg_fallback"});
        series.writeCouplingFallback(1,.25,timing,summary);
        series.writePerformance(2,.5,zero,zeroSummary,{0,false,true,false,0,0,"direct"});
        series.writeCouplingFallback(2,.5,zero,zeroSummary);
        series.flushPerformance();
        const auto rows=csvRows("perf_implicit.csv"),fallback=csvRows("fallback.csv");
        check(rows.size()==3 && rows[0].size()==45 && rows[1].size()==45 &&
            rows[0][1]=="time_s" && rows[1][2]=="4" && rows[1][24]=="1" &&
            rows[1][25]=="0" && rows[1][26]=="1" && rows[1][29]=="cg_fallback",
            "implicit timing schema/flags/solver mode");
        const std::array<double,9> ratios{{6.,1.,3.,3.,4.,2.,1.,.5,3.}};
        for(std::size_t j=0;j<ratios.size();++j) {
            check(std::stod(rows[1][30+j])==ratios[j],"implicit performance normalized counter");
            check(std::stod(rows[2][30+j])==0.,"zero-particle timing counter must remain finite");
        }
        check(rows[1][39]=="amr" && rows[1][41]=="output" && rows[1][43]=="mcc_pre",
            "timing CSV ranking differs from console summary");
        check(fallback.size()==3 && fallback[0].size()==9 && fallback[1][2]=="2" &&
            fallback[1][3]=="8" && std::stod(fallback[1][4])==.25 && fallback[1][8]=="5" &&
            fallback[2][8]=="0","fallback CSV counts/fractions");
    }
    {
        PicStepSeries series(Mode::Implicit);
        series.openPerformance("perf_implicit.csv",2);
        series.openCouplingFallback("fallback.csv",2);
        series.writePerformance(2,.75,timing,summary,{4,false,true,false,8,32,"cg"});
        series.writeCouplingFallback(2,.75,timing,summary);
        series.flushPerformance();
        const auto rows=csvRows("perf_implicit.csv"),fallback=csvRows("fallback.csv");
        check(rows.size()==3 && std::stod(rows[2][1])==.75 && rows[2][29]=="cg" &&
            fallback.size()==3 && std::stod(fallback[2][1])==.75,"implicit timing resume cadence");
    }
}
void checkPhysicsDiagnostics() {
    using namespace PhysicsDiagnostics;
    ElectrodeCurrentWriter unopenedElectrodes;ReactionDiagnosticsWriter unopenedReactions;
    EnergyHistogramWriter unopenedEnergy;
    check(!unopenedElectrodes.isOpen() && !unopenedReactions.isOpen() && !unopenedEnergy.isOpen(),
        "default streams must not activate diagnostic sampling or skip lazy open");
    unopenedReactions.writeStep(1,0.,1.,{});unopenedEnergy.writeStep(1,0.);
    check(!unopenedReactions.isOpen() && !unopenedEnergy.isOpen(),"unopened diagnostic write activated a stream");
    auto groups=makeElectrodeCurrentGroups({{"first",{-1,0,1,1,7}},{"second",{1,2}}});
    groups.record(0,99.,99.);groups.record(1,-4.,0.);groups.record(2,6.,0.);
    groups.record(7,0.,3.);groups.record(8,99.,99.);
    check(groups.buckets()[0].chargeC==-4. && groups.buckets()[0].electronChargeC==-4. &&
        groups.buckets()[0].neutralCountPhys==3. && groups.buckets()[1].ionChargeC==6.,
        "group ownership, invalid IDs, or signed charge classification changed");
    check(groups.indexOf(1)==0 && groups.indexOf(8)==groups.groups().size(),"first group owns repeated boundary");
    {
        ElectrodeCurrentWriter writer;check(writer.open("physics_groups",groups.groups(),true),"electrode open");
        writer.beginStep();writer.record(1,-4.,0.);writer.record(2,6.,0.);writer.record(7,0.,3.);
        writer.record(8,2.,1.);writer.writeStep(1,.5,.5);
        const auto rows=csvRows("physics_groups/electrodes/current_groups_step.csv");
        check(rows.size()==3 && rows[1][2]=="first" && rows[1][3]=="5" &&
            std::stod(rows[1][5])==-8. && std::stod(rows[2][9])==12.,"group CSV charge/current/count");
        const auto boundaries=csvRows("physics_groups/electrodes/current_boundaries_step.csv");
        bool ungrouped=false;for(const auto& row:boundaries)if(row[2]=="8")ungrouped=row[3]=="ungrouped";
        check(ungrouped,"per-boundary CSV must preserve unmatched positive boundaries");
    }
    groups.beginStep();check(groups.buckets()[0].chargeC==0. && groups.buckets()[1].ionChargeC==0.,"group reset");
    ReactionDatabase reactions;ReactionSet set;set.name="electron-Ar";set.projectile="electron";set.target="Ar";
    ReactionChannel channel;channel.name="b";channel.cross_section_key="dummy";channel.products={"ion","ion",""};
    set.channels.push_back(channel);channel.name="a";set.channels.push_back(channel);
    check(reactions.Register(set),"diagnostic reaction fixture");
    MCCCollisionOperator::StepStats stats;
    stats.n_particles=2;stats.n_particles_weighted=10.;stats.n_attempts=4;stats.n_attempts_weighted=20.;
    stats.n_collisions=3;stats.n_collisions_weighted=15.;stats.n_spawned=1;stats.n_spawned_weighted=5.;
    stats.energy_loss_J_total=4.;stats.per_channel["electron::a"]=3;
    stats.per_channel_weighted["electron::a"]=15.;stats.per_channel_energy_loss_J["electron::a"]=3.;
    stats.spawn_per_species["ion"]=1;stats.spawn_per_species_weighted["ion"]=5.;
    OutputConfig::DiagnosticsConfig::ReactionStatsConfig config;
    {
        ReactionDiagnosticsWriter writer;
        check(writer.open("physics_reactions",config,&reactions,nullptr) && writer.isOpen(),"reaction first open");
        writer.writeStep(1,.5,.5,stats);writer.writeStep(3,1.5,.5,stats);
        const auto rows=csvRows("physics_reactions/reactions/reaction_step.csv");
        check(rows.size()==3 && rows[1][0]=="1" && rows[1][8]=="3" &&
            std::stod(rows[1][10])==30. && std::stod(rows[1][13])==10.,"reaction first-step counters/frequencies");
        const auto channels=csvRows("physics_reactions/reactions/reaction_channel_step.csv");
        check(channels[0][2]=="electron_a_macro_collisions" && channels[0][5]=="electron_b_macro_collisions" &&
            channels[1][2]=="3" && channels[1][5]=="0","reaction column selection/order/zero fill");
        const auto energy=csvRows("physics_reactions/reactions/reaction_energy_step.csv");
        check(std::stod(energy[1][4])==3. && std::stod(energy[1][5])==6. && std::stod(energy[1][6])==.75,
            "reaction channel energy, power and fraction");
        const auto spawn=csvRows("physics_reactions/reactions/reaction_spawn_step.csv");
        check(spawn[0].size()==5 && spawn[1][2]=="1" && std::stod(spawn[1][4])==10.,"deduplicated spawn columns");
    }
    {
        ReactionDiagnosticsWriter writer;check(writer.open("physics_reactions",config,&reactions,nullptr,3),"reaction resume");
        writer.writeStep(3,2.,1.,stats);
        const auto rows=csvRows("physics_reactions/reactions/reaction_step.csv");
        check(rows.size()==3 && rows[1][0]=="1" && std::stod(rows[2][1])==2.,"reaction resume preserves first and trims future");
    }
    {
        ReactionDiagnosticsWriter writer;config.summary=config.perChannel=config.perProjectile=config.spawnedSpecies=false;
        check(writer.open("physics_disabled",config,nullptr,nullptr) && !writer.isOpen(),"all-disabled reaction streams falsely ready");
        writer.writeStep(1,0.,1.,stats);
        check(std::filesystem::is_empty("physics_disabled/reactions"),"disabled reactions created output files");
    }
    {
        EnergyHistogramWriter writer;OutputConfig::DiagnosticsConfig::EnergyHistogramConfig config;
        config.bins=2;config.maxEnergyEV=10.;std::vector<SpeciesConfig> species(1);species[0].name="electron";
        check(writer.open("physics_energy",config,species) && writer.isOpen(),"energy configured open");
        writer.beginStep();writer.observe("electron",-1.,"bulk",4.9,3.);
        writer.observe("electron",-1.,"bulk",5.,7.);writer.observe("electron",-1.,"bulk",20.,2.);
        writer.observe("electron",-1.,"domain",5.,99.);writer.writeStep(1,.5);
        const auto rows=csvRows("physics_energy/energy/histogram_step.csv");
        check(rows.size()==3 && rows[1][2]=="electron_bulk" && std::stod(rows[1][8])==3. &&
            std::stod(rows[2][8])==9.,"default energy samples, bins and bulk location filter");
        writer.beginStep();writer.writeStep(2,1.);
        check(std::stod(csvRows("physics_energy/energy/histogram_step.csv")[3][8])==0.,"histogram reset");
    }
    // Nonempty directories survive legacy reset and make the target open fail.
    for(const std::string leaf:{"electrodes/current_groups_step.csv","reactions/reaction_step.csv","energy/histogram_step.csv"}) {
        std::filesystem::create_directories(std::filesystem::path("physics_blocked")/leaf);
        std::ofstream(std::filesystem::path("physics_blocked")/leaf/"keep")<<"keep";
    }
    ElectrodeCurrentWriter blockedElectrodes;ReactionDiagnosticsWriter blockedReactions;EnergyHistogramWriter blockedEnergy;
    check(!blockedElectrodes.open("physics_blocked",{},false) && !blockedElectrodes.isOpen(),"failed electrode falsely ready");
    check(!blockedReactions.open("physics_blocked",true,false,false,false,{},{},{}) && !blockedReactions.isOpen(),
        "failed reaction falsely ready through an unopened sibling stream");
    check(!blockedEnergy.open("physics_blocked",2,10.,{}) && !blockedEnergy.isOpen(),"failed histogram falsely ready");
}
void checkPhysicsDiagnosticStage() {
    using namespace PhysicsDiagnostics;using Output=PicDiagnosticsOutput;using Mode=Output::Mode;
    static_assert(!std::is_copy_constructible_v<Output> && !std::is_move_constructible_v<Output>);
    SpeciesDatabase speciesDb;GasState gas;CrossSectionDB crossSections;ReactionDatabase reactions;
    ReactionSet set;set.name="electron-Ar";set.projectile="electron";set.target="Ar";
    ReactionChannel channel;channel.name="ionization";channel.products={"electron","ion"};set.channels={channel};
    check(reactions.Register(set),"diagnostic stage reaction fixture");
    MCCCollisionOperator operation(speciesDb,gas,crossSections,reactions);
    std::vector<SpeciesConfig> species(2);species[0].name="electron";species[1].name="ion";
    MCCCollisionOperator::StepStats stats;stats.n_collisions=2;stats.n_collisions_weighted=6.;
    stats.per_channel["electron::ionization"]=2;stats.per_channel_weighted["electron::ionization"]=6.;
    stats.spawn_per_species["ion"]=2;stats.spawn_per_species_weighted["ion"]=6.;
    const auto compareTree=[](const std::filesystem::path& a,const std::filesystem::path& b) {
        std::map<std::string,std::string> left,right;
        if(std::filesystem::exists(a))for(const auto& p:std::filesystem::recursive_directory_iterator(a))
            if(p.is_regular_file())left[std::filesystem::relative(p.path(),a).generic_string()]=bytes(p.path());
        if(std::filesystem::exists(b))for(const auto& p:std::filesystem::recursive_directory_iterator(b))
            if(p.is_regular_file())right[std::filesystem::relative(p.path(),b).generic_string()]=bytes(p.path());
        check(left==right,"diagnostic lifecycle CSVs match original independent writers exactly");
    };
    int combination=0;
    for(Mode mode:{Mode::Explicit,Mode::Implicit})for(int flags=0;flags<8;++flags)for(int cadence=0;cadence<3;++cadence) {
        const auto base=std::filesystem::path("physics_stage")/std::to_string(combination++);
        OutputConfig::DiagnosticsConfig config;config.enabled=false; // Individual parsed flags govern legacy IO.
        config.electrodeCurrent.enabled=(flags&1)!=0;config.reactions.enabled=(flags&2)!=0;config.energyHistogram.enabled=(flags&4)!=0;
        config.electrodeCurrent.includeDomain=true;config.electrodeCurrent.perBoundary=true;
        config.electrodeCurrent.groups={{"shared",{1,2,2}},{"other",{2,3}}};
        config.energyHistogram.bins=3;config.energyHistogram.maxEnergyEV=12.;
        config.energyHistogram.samples={{"all_bulk","all","bulk"},{"wall_e","electron","domain"},{"embedded_i","ion","embedded"}};
        config.electrodeCurrent.every=cadence==2?3:0;config.reactions.every=cadence==2?2:-1;config.energyHistogram.every=cadence==2?4:0;
        const int baseEvery=cadence==1?0:2;
        const auto execute=[&](int resumeStep,int first,int last) {
            auto settings=config;
            Output output(mode,settings,base/"actual",baseEvery,resumeStep,species,&reactions,nullptr);
            settings.electrodeCurrent.groups.clear();settings.electrodeCurrent.includeDomain=false;settings.energyHistogram.bins=999;
            check(output.includesDomainElectrodes() && output.consoleElectrodes().groups().size()==2,"diagnostic owner snapshots configuration");
            ElectrodeCurrentWriter electrodes;ReactionDiagnosticsWriter reactionWriter;EnergyHistogramWriter energy;
            auto console=makeElectrodeCurrentGroups(config.electrodeCurrent.groups);
            if(config.electrodeCurrent.enabled)check(electrodes.open(base/"reference",console.groups(),true,resumeStep),"reference electrode open");
            if(mode==Mode::Implicit && config.reactions.enabled)check(reactionWriter.open(base/"reference",config.reactions,&reactions,nullptr,resumeStep),"reference implicit reactions without MCC");
            if(config.energyHistogram.enabled)check(energy.open(base/"reference",config.energyHistogram,species,resumeStep),"reference histogram open");
            check(!output.sampleElectrodes() && !output.sampleEnergy(),"sampling remains inactive until beginStep");
            for(int step=first;step<=last;++step) {
                const MCCCollisionOperator* op=step>=2?&operation:nullptr;
                if(mode==Mode::Explicit && config.reactions.enabled && !reactionWriter.isOpen() && op)
                    check(reactionWriter.open(base/"reference",config.reactions,&reactions,op,resumeStep),"reference delayed reaction open");
                const auto due=[&](bool ready,int every){return ready && every>0 && step%every==0;};
                const bool e=due(electrodes.isOpen(),resolveEvery(baseEvery,config.electrodeCurrent.every));
                const bool r=due(reactionWriter.isOpen(),resolveEvery(baseEvery,config.reactions.every));
                const bool h=due(energy.isOpen(),resolveEvery(baseEvery,config.energyHistogram.every));
                output.beginStep(step,&reactions,op);
                check(output.sampleElectrodes()==e && output.sampleEnergy()==h,"diagnostic cadence uses healthy streams and positive intervals");
                if(e)electrodes.beginStep();console.beginStep();if(h)energy.beginStep();
                for(int boundary:{0,1,2,3,9}) {
                    const double q=boundary%2?-step:step,w=boundary==3?2.:0.;
                    output.recordElectrode(boundary,q,w);if(e)electrodes.record(boundary,q,w);console.record(boundary,q,w);
                }
                for(std::size_t i=0;i<console.buckets().size();++i) {
                    const auto& a=console.buckets()[i];const auto& b=output.consoleElectrodes().buckets()[i];
                    check(a.chargeC==b.chargeC && a.electronChargeC==b.electronChargeC && a.ionChargeC==b.ionChargeC && a.neutralCountPhys==b.neutralCountPhys,"console buckets reset each step independently of file cadence");
                }
                for(const char* location:{"domain","embedded","bulk"})for(std::size_t si=0;si<species.size();++si) {
                    const double charge=si?1.:-1.;
                    for(double value:{0.,3.9,4.,12.,40.}) {
                        output.observeEnergy(species[si].name,charge,location,value,double(step));
                        if(h)energy.observe(species[si].name,charge,location,value,double(step));
                    }
                }
                output.writeReactions(step,.25*step,.25,stats);if(r)reactionWriter.writeStep(step,.25*step,.25,stats);
                output.writeEnergy(step,.25*step);if(h)energy.writeStep(step,.25*step);
                output.writeElectrodes(step,.25*step,.25);if(e)electrodes.writeStep(step,.25*step,.25);
                compareTree(base/"actual",base/"reference");
            }
        };
        execute(0,1,7);execute(3,3,5);compareTree(base/"actual",base/"reference");
    }
    // Original failed-open policy: explicit retries reaction output at each
    // original prepare/begin phase; implicit attempts it only at construction.
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        const auto root=std::filesystem::path("physics_stage_blocked")/(mode==Mode::Explicit?"explicit":"implicit");
        for(const char* name:{"electrodes/current_groups_step.csv","reactions/reaction_step.csv","energy/histogram_step.csv"}) {
            std::filesystem::create_directories(root/name);std::ofstream(root/name/"keep")<<"keep";
        }
        OutputConfig::DiagnosticsConfig config;config.electrodeCurrent.enabled=config.reactions.enabled=config.energyHistogram.enabled=true;
        config.electrodeCurrent.groups={{"console",{1}}};config.reactions.perChannel=config.reactions.perProjectile=config.reactions.spawnedSpecies=false;
        struct Capture {std::ostringstream text;std::streambuf* old=std::cerr.rdbuf(text.rdbuf());~Capture(){std::cerr.rdbuf(old);}} capture;
        Output output(mode,config,root,1,0,species,&reactions,&operation);
        output.prepareExplicitReactions(&reactions,&operation);
        output.beginStep(1,&reactions,&operation);output.recordElectrode(1,-2.,0.);
        check(!output.sampleElectrodes() && !output.sampleEnergy() && output.consoleElectrodes().buckets()[0].chargeC==-2.,"failed file IO keeps console accounting active");
        output.beginStep(2,&reactions,&operation);
        const std::string prefix=mode==Mode::Explicit?"[ExplicitPIC]":"[ImplicitPIC]";
        const auto warning=[&](const char* kind){return prefix+" WARNING: failed to create "+kind+" diagnostics under "+root.string()+"\n";};
        const std::string expected=mode==Mode::Explicit
            ?warning("electrode")+warning("energy histogram")+warning("reaction")+warning("reaction")+warning("reaction")
            :warning("electrode")+warning("reaction")+warning("energy histogram");
        check(capture.text.str()==expected,"failed diagnostic open warning order and explicit-only retry policy");
    }
    OutputConfig::DiagnosticsConfig onlyReactions;onlyReactions.reactions.enabled=true;
    Output noMcc(Mode::Explicit,onlyReactions,"physics_stage_no_mcc",1,0,species,&reactions,nullptr);
    noMcc.beginStep(1,&reactions,nullptr);noMcc.writeReactions(1,1.,1.,stats);
    check(!std::filesystem::exists("physics_stage_no_mcc"),"explicit reaction-only output stays absent without MCC");
}
void checkMaskDiagnostics(const QuadtreeMesh& mesh) {
    CellField<double> phi(mesh,1.23456789);
    std::vector<unsigned char> valid(mesh.cellIndexNum(),1);valid[0]=0;
    std::vector<double> alpha(mesh.cellIndexNum(),.3);
    alpha[0]=-1.;alpha[1]=2.;alpha[2]=std::numeric_limits<double>::quiet_NaN();
    alpha[3]=std::numeric_limits<double>::infinity();
    std::vector<LevelSet> levelSets;levelSets.emplace_back(mesh);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)
        levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.5;
    FaceSegmentGraph graph;graph.build(mesh);
    for(int variant=0;variant<3;++variant) {
        const std::vector<double> wrongAlpha{0.};
        PicFieldOutput::writeMaskDiagnostics(mesh,20+variant,valid,phi,graph,nullptr,
            variant==2?&wrongAlpha:&alpha,variant==2?nullptr:&levelSets,{-1,0,99},variant==0);
        const auto rows=csvRows("output/meta/packed_mask_diag_"+std::to_string(20+variant)+".csv");
        check(rows.size()==static_cast<std::size_t>(mesh.cellIndexNum()+1),"diagnostic must include invalid output cells");
        check(rows[0]==std::vector<std::string>({"cellIndex","x_m","y_m","isLeaf","outputValid","cutActive",
            "faceSegCount","poissonActive","alpha","allInsideAnyLS","phiCell","isCut","A_gamma","d_gamma",
            "gammaBoundaryId","gammaLineIndex"}),"mask diagnostic schema");
        int insideCount=0;
        for(int ci=0;ci<mesh.cellIndexNum();++ci) {
            const auto& row=rows[ci+1];const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);
            check(row.size()==16 && std::stoi(row[0])==ci,"mask diagnostic row order/width");
            check(std::stod(row[1])==cell.centerCoord().x && std::stod(row[2])==cell.centerCoord().y,
                "diagnostic retains raw mesh coordinates despite its SI column labels");
            check(row[3]=="1" && row[4]==(ci==0?"0":"1") && row[5]=="0","mask flags and absent cut metrics");
            check(row[7]=="1" && row[10]=="1.23457","face activity and default stream precision");
            const double expectedAlpha=variant==2?1.:(ci==0?0.:(ci<=3?1.:.3));
            check(std::stod(row[8])==expectedAlpha,"diagnostic alpha clamp/non-finite/size fallback");
            const bool inside=variant==0 ? cell.centerCoord().x+cell.size().x*.5<.5
                : variant==1 && cell.centerCoord().x-cell.size().x*.5>.5;
            check(row[9]==(inside?"1":"0"),"strict level-set sign at all four corners, invalid index/null geometry handling");
            insideCount+=inside?1:0;
            check(row[11]=="0" && row[12]=="0" && row[13]=="0" && row[14]=="0" && row[15]=="-1",
                "absent cut metric diagnostic defaults");
        }
        if(variant<2)check(insideCount>0,"both sign policies need an inside cell");
    }
    std::filesystem::create_directory("output/meta/packed_mask_diag_23.csv");
    PicFieldOutput::writeMaskDiagnostics(mesh,23,valid,phi,graph,nullptr,nullptr,nullptr,{},true);
    check(std::filesystem::is_directory("output/meta/packed_mask_diag_23.csv"),"diagnostic open failure remains silent");
    check(phi.value[0]==1.23456789 && valid[0]==0 && alpha[0]==-1.,"mask diagnostic modified physical inputs");
}
void checkProfiles(const QuadtreeMesh& mesh) {
    const auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(UnitsScales{}));
    const auto norm=NormalizationContext::MakeFromUnits(units,NormalizationConfig{});
    CellField<double> phi(mesh,1.),rho(mesh,0.);CellField<Vector2D> electric(mesh,Vector2D{2.,0.});
    std::vector<ParticleType> species;species.emplace_back("electron",2.*PhysicalConstants::eV_J,-1.);
    const std::vector<double> weights{1.};std::vector<CouplingWeightsCache> caches(1);
    caches[0].buildFromCSR(0,{0},{},{});
    const PicProfiles::Sample sample{mesh,units.get(),&norm,phi,electric,species,weights,nullptr};
    auto setPhi=[&](double value) {for(int ci=0;ci<mesh.cellIndexNum();++ci)phi.value[ci]=value;};
    int calls=0;
    auto build=[&](std::size_t si,CouplingWeightsCache& cache) {check(si==0,"profile species callback index");++calls;cache=caches[si];};
    for(auto mode:{PicFieldMode::Explicit,PicFieldMode::Implicit}) {
        const std::string dir=mode==PicFieldMode::Explicit?"profiles/explicit":"profiles/implicit";
        const auto average=dir+"/profiles_avg.csv";setPhi(1.);calls=0;
        std::string previousAverage;
        {
            PicProfiles output(mode,dir,0,2,{true,2,true,2});
            check(output.instantaneousDue(0) && !output.instantaneousDue(1) && output.instantaneousDue(2),"profile output cadence");
            output.writeInstantaneous(1,.1,.1,sample,rho,build,true);check(calls==0,"off-cadence profile must not invoke cache builder");
            output.accumulate(1,sample,caches);check(!std::filesystem::exists(average),"average requires a full window");
            setPhi(3.);output.accumulate(2,sample,caches);output.writeInstantaneous(2,.2,.2,sample,rho,build,true);
            check(calls==2,"instantaneous profile rebuilds density and electron energy, including empty species");
            const auto first=csvRows(average);check(first.size()==3,"profile average configured bin count");
            for(std::size_t row=1;row<first.size();++row)
                check(std::stod(first[row][2])==2. && std::stod(first[row][3])==1.,"profile rolling mean/population deviation");
            const auto initial=bytes(average);setPhi(5.);output.accumulate(3,sample,caches);
            check(bytes(average)==initial,"off-cadence average must not rewrite its snapshot");
            output.accumulate(4,sample,{});check(bytes(average)==initial,"cache-count mismatch must skip accumulation and writing");
            output.writeInstantaneous(4,.4,.4,sample,rho,build,true);
            setPhi(9.);output.accumulate(6,sample,caches);
            const auto last=csvRows(average);
            for(std::size_t row=1;row<last.size();++row)
                check(std::stod(last[row][2])==7. && std::stod(last[row][3])==2.,"window replacement or skipped-cache sample changed average");
            previousAverage=bytes(average);
        }
        auto instant=csvRows(dir+"/profile1d.csv");check(instant.size()==5 && instant[1][0]=="2" && instant[3][0]=="4","instantaneous rows and cadence");
        {
            PicProfiles resumed(mode,dir,2,2,{true,2,true,2});setPhi(11.);
            resumed.accumulate(4,sample,caches);
            check(bytes(average)==previousAverage,"restart must reset the average window, not rewrite it before ready");
            resumed.accumulate(6,sample,caches);resumed.writeInstantaneous(4,.4,.4,sample,rho,build,true);
        }
        instant=csvRows(dir+"/profile1d.csv");check(instant.size()==3 && instant[1][0]=="4" && instant[2][0]=="4",
            "resume must trim the boundary step and later rows, then preserve one header");
        const auto reset=csvRows(average);
        for(std::size_t row=1;row<reset.size();++row)
            check(std::stod(reset[row][2])==11. && std::stod(reset[row][3])==0.,"resumed averages retained old samples");
        {
            PicProfiles defaults(mode,dir+"/defaults",0,1,{true,0,true,1});
            defaults.accumulate(1,sample,caches);defaults.writeInstantaneous(1,.1,.1,sample,rho,build,false);
        }
        check(csvRows(dir+"/defaults/profile1d.csv").size()==2 && csvRows(dir+"/defaults/profiles_avg.csv").size()==129,
            "instantaneous nonpositive bins clamp to one; average retains 128-bin fallback");
    }
    calls=0;
    PicProfiles disabled(PicFieldMode::Explicit,"profiles/disabled",0,1,{false,2,true,0});
    disabled.accumulate(1,sample,caches);disabled.writeInstantaneous(1,.1,.1,sample,rho,build,true);
    check(!disabled.instantaneousDue(1) && calls==0 && !std::filesystem::exists("profiles/disabled"),"disabled profiles must not create files or build caches");
    std::filesystem::create_directories("profiles/blocked/profile1d.csv/keep");
    std::ostringstream warning;auto* old=std::cerr.rdbuf(warning.rdbuf());
    {PicProfiles blocked(PicFieldMode::Implicit,"profiles/blocked",0,1,{true,2,false,0});check(!blocked.instantaneousDue(1),"failed stream must disable instant writes");}
    std::cerr.rdbuf(old);check(warning.str().find("[ImplicitPIC] WARNING: failed to create profile1d csv:")!=std::string::npos,"profile-open warning identity");
}
void checkSpeciesSeries() {
    using Series=PicSpeciesSeries;
    const std::vector<Series::Species> species{{"A-B",-1.},{"a b",1.},{"a_b_2",0.},
        {"PhotoElectron",1.},{"!!!",0.},{"...",0.}};
    const std::vector<std::string> tokens{"a_b","a_b_2","a_b_2_2","photoelectron","species","species_2"};
    const std::vector<long long> macro{9007199254740993LL,7};
    const std::vector<double> physical{1.25,2.5},density{3.25,4.5};
    const Series::Moments moments{2,macro,physical,density,6.,4.};
    const Series::BudgetRow budget{{9007199254740993LL,.25},{2,2.25},{3,3.25},{4,4.25},{5,5.25},
        {6,6.25},{7,7.25},{8,8.25},{-9,-9.25},{-10,-10.25}};
    for(auto mode:{Series::Mode::Explicit,Series::Mode::Implicit}) {
        const std::string dir=(mode==Series::Mode::Explicit)?"species/explicit/":"species/implicit/";
        Series output(mode,species);
        check(!output.primaryReady() && !output.detailReady() && !output.budgetReady(),"unopened CSV readiness");
        output.openPrimary(dir+"primary.csv",0);output.openDetail(dir+"detail.csv",0);output.openBudget(dir+"budget.csv",0);
        check(output.primaryReady() && output.detailReady() && output.budgetReady(),"opened CSV readiness");
        if(mode==Series::Mode::Implicit)check(csvRows(dir+"primary.csv").size()==1,"implicit header must be flushed immediately");
        for(int step=1;step<=3;++step) {
            output.writePrimary(step,step*.25,moments);output.writeDetail(step,step*.25,moments);
            output.beginBudget(step,step*.25);for(std::size_t si=0;si<species.size();++si)output.addBudget(budget);output.endBudget();
        }
        check(bytes(dir+"primary.csv")==bytes(dir+"detail.csv"),"primary/detail schema or formatting diverged");
        auto rows=csvRows(dir+"primary.csv");check(rows.size()==4 && rows[0].size()==32,"species schema width");
        for(std::size_t si=0;si<tokens.size();++si) {
            check(rows[0][2+5*si]==tokens[si]+"_macro_count" && rows[0][6+5*si]==tokens[si]+"_Te_eV","species column order");
            check(std::stod(rows[1][5+5*si])==(si==3?6.:0.) && std::stod(rows[1][6+5*si])==(si==3?4.:0.),
                "electron name priority or electron summary on a missing species changed");
            if(si>=2)for(int j=0;j<3;++j)check(std::stod(rows[1][2+5*si+j])==0.,"missing species zero fallback");
        }
        check(rows[1][2]=="9007199254740993" && rows[1][1]=="0.2500000000000000","integer counts or fixed precision lost");
        check(rows[1][3]=="1.2500000000000000" && rows[1][4]=="3.2500000000000000","prepared physical values changed");
        auto budgets=csvRows(dir+"budget.csv");check(budgets.size()==4 && budgets[0].size()==122,"budget schema width");
        for(std::size_t si=0;si<tokens.size();++si) {
            check(budgets[0][2+20*si]==tokens[si]+"_macro_begin" && budgets[0][21+20*si]==tokens[si]+"_phys_residual","budget column order");
            check(budgets[1][2+20*si]=="9007199254740993" && budgets[1][3+20*si]=="0.2500000000000000","budget count types");
            for(int pair=1;pair<10;++pair) {
                const double expected=pair<8?pair+1.:-(pair+1.);
                check(std::stod(budgets[1][2+20*si+2*pair])==expected
                    && std::stod(budgets[1][3+20*si+2*pair])==expected+(expected>0?.25:-.25),"budget pairs reordered");
            }
        }
    }
    check(bytes("species/explicit/primary.csv")==bytes("species/implicit/primary.csv"),"integrator CSV formats differ");
    const auto original=bytes("species/implicit/primary.csv");
    {
        std::ofstream file("species/implicit/primary.csv",std::ios::app);file<<"# preserved comment\n999,stale\n";
    }
    {
        Series resumed(Series::Mode::Implicit,species);resumed.openPrimary("species/implicit/primary.csv",2);
        resumed.writePrimary(2,.5,moments);
        const auto rows=csvRows("species/implicit/primary.csv");
        check(rows.size()==4 && rows[0][0]=="step" && rows[1][0]=="1" && rows[2][0]=="# preserved comment" && rows[3][0]=="2",
            "resume must preserve earlier rows/comments, prune later rows and avoid a duplicate header");
    }
    {
        Series fresh(Series::Mode::Implicit,species);fresh.openPrimary("species/implicit/primary.csv",0);
        fresh.writePrimary(1,.25,moments);
        const auto rows=csvRows("species/implicit/primary.csv");check(rows.size()==2,"fresh CSV did not replace prior contents");
        check(original.rfind(bytes("species/implicit/primary.csv"),0)==0,"fresh header/first row changed");
    }
    for(int variant=0;variant<2;++variant) {
        const auto path="species/fallback"+std::to_string(variant)+".csv";
        Series output(Series::Mode::Implicit,{{"other",variant?1.:-1.},{"also",-double(!variant)}});
        output.openPrimary(path,0);output.writePrimary(1,.25,moments);
        const auto row=csvRows(path)[1];check(std::stod(row[5])==(variant?0.:6.) && std::stod(row[10])==0.,"negative-charge fallback or absent electron");
    }
    std::filesystem::create_directories("species/blocked.csv/child");
    Series failed(Series::Mode::Implicit,species);
    std::ostringstream warning;auto* previous=std::cerr.rdbuf(warning.rdbuf());
    failed.openDetail("species/blocked.csv",0);std::cerr.rdbuf(previous);
    check(!failed.detailReady() && warning.str().find("[ImplicitPIC] WARNING: failed to create species detail csv:")!=std::string::npos,
        "original CSV-open warning/readiness");
    failed.openPrimary("species/usable.csv",0);failed.writePrimary(1,.25,moments);
    check(csvRows("species/usable.csv").size()==2,"failed detail stream changed primary ownership");
    {
        Series unfinished(Series::Mode::Implicit,{{"e",-1.}});unfinished.openBudget("species/partial.csv",0);
        unfinished.beginBudget(1,.25);unfinished.addBudget(budget);
    }
    check(bytes("species/partial.csv").back()!='\n',"destruction must not synthesize a completed budget row");
}
void checkBoundarySeries() {
    using Series=PicSpeciesSeries;
    const std::uint64_t largeUnsigned=18446744073709551000ULL;
    const long long largeSigned=9007199254740993LL;
    const std::vector<Series::Species> species{{"e-",-1.},{"e ",1.}};
    for(auto mode:{Series::Mode::Explicit,Series::Mode::Implicit})for(bool si:{false,true}) {
        const bool implicit=mode==Series::Mode::Implicit;
        const auto dir=std::string("boundaries/")+(implicit?"implicit":"explicit")+(si?"_si/":"_code/");
        {
            Series output(mode,species);
            check(!output.boundaryReady() && !output.boundaryEnergyReady() && !output.fluxReady(),"unopened boundary readiness");
            output.openBoundary(dir+"primary.csv",0);output.openBoundaryEnergy(dir+"energy.csv",0);output.openFlux(dir+"flux.csv",0,si);
            check(output.boundaryReady() && output.boundaryEnergyReady() && output.fluxReady(),"opened boundary readiness");
            if(implicit)check(csvRows(dir+"flux.csv").size()==1,"implicit flux header must flush immediately");
            for(int step=1;step<=2;++step) {
                output.beginBoundary(step,step*.25);output.beginBoundaryEnergy(step,step*.25);
                for(int speciesIndex=0;speciesIndex<2;++speciesIndex)for(int side=0;side<5;++side) {
                    const int offset=100*speciesIndex+side;
                    if(side<4) {
                        if(implicit)output.addBoundary(Series::IntegerFlux{largeUnsigned+offset,1.25,2.5,3.75,-1e-12});
                        else output.addBoundary(Series::FloatingFlux{double(offset+3),1.25,2.5,3.75,-1e-12});
                    }
                    output.addBoundaryEnergy({{largeSigned+offset,1.25},2.5,3.75,4.5,
                        side==4?std::numeric_limits<double>::quiet_NaN():5.25,side==4?0.:6.5});
                }
                output.endBoundary();output.endBoundaryEnergy();
                if(implicit)output.writeFlux(step,step*.25,"left","e raw-name",Series::IntegerFlux{largeUnsigned,1.25,2.5,3.75,-1e-12});
                else output.writeFlux(step,step*.25,"L","e raw-name",Series::FloatingFlux{3.,1.25,2.5,3.75,-1e-12});
                output.flushFlux();
            }
            const auto boundary=csvRows(dir+"primary.csv"),energy=csvRows(dir+"energy.csv"),flux=csvRows(dir+"flux.csv");
            check(boundary.size()==3 && boundary[0].size()==42 && energy[0].size()==72 && flux[0].size()==9,"boundary schema widths");
            const std::array<std::string,5> sides{"l","r","b","t","embedded"};
            for(int sp=0;sp<2;++sp)for(int side=0;side<5;++side) {
                const auto prefix=std::string(sp?"e_2_":"e_")+sides[side];const int offset=100*sp+side;
                if(side<4) {
                    const auto index=2+20*sp+5*side;
                    check(boundary[0][index]==prefix+"_macro_out" && boundary[0][index+4]==prefix+"_current_Apm2","boundary side/species order");
                    const auto expected=implicit?std::to_string(largeUnsigned+offset):std::to_string(offset+3)+".0000000000000000";
                    check(boundary[1][index]==expected && boundary[1][index+4]=="-0.0000000000010000","boundary count type or current format");
                }
                const auto index=2+35*sp+7*side;
                check(energy[0][index]==prefix+"_macro_out" && energy[0][index+6]==prefix+"_surface_area_m2","energy side/species order");
                check(energy[1][index]==std::to_string(largeSigned+offset),"energy signed integer precision");
                for(int j=1;j<5;++j)check(std::stod(energy[1][index+j])==std::array<double,4>{1.25,2.5,3.75,4.5}[j-1],"energy values reordered");
                if(side==4)check(std::isnan(std::stod(energy[1][index+5])) && std::stod(energy[1][index+6])==0.,"embedded undefined flux/zero area");
                else check(std::stod(energy[1][index+5])==5.25 && std::stod(energy[1][index+6])==6.5,"energy flux/area order");
            }
            check(flux[0][1]==((implicit||si)?"time_s":"time"),"flux unit-dependent header");
            check(flux[0][6]==((implicit||si)?"flux_macro_1pm2s":"flux_macro"),"flux units label");
            check(flux[1][2]==(implicit?"left":"L") && flux[1][3]=="e raw-name","raw flux side/species names");
            check(flux[1][1]==(implicit?"0.2500000000000000":"0.25") && flux[1][4]==(implicit?std::to_string(largeUnsigned):"3"),"long-form numeric format");
        }
        {
            Series resumed(mode,species);resumed.openFlux(dir+"flux.csv",2,si);
            if(implicit)resumed.writeFlux(2,.5,"left","e raw-name",Series::IntegerFlux{largeUnsigned,1.25,2.5,3.75,-1e-12});
            else resumed.writeFlux(2,.5,"L","e raw-name",Series::FloatingFlux{3.,1.25,2.5,3.75,-1e-12});
            resumed.flushFlux();const auto rows=csvRows(dir+"flux.csv");
            check(rows.size()==3 && rows[1][0]=="1" && rows[2][0]=="2","flux resume boundary or duplicate header");
        }
        for(const auto& role:{"primary","energy","flux"})std::filesystem::create_directories(dir+role+"_blocked.csv/child");
        Series failed(mode,species);std::ostringstream warnings;auto* old=std::cerr.rdbuf(warnings.rdbuf());
        failed.openBoundary(dir+"primary_blocked.csv",0);failed.openBoundaryEnergy(dir+"energy_blocked.csv",0);failed.openFlux(dir+"flux_blocked.csv",0,si);
        std::cerr.rdbuf(old);
        check(!failed.boundaryReady() && !failed.boundaryEnergyReady() && !failed.fluxReady(),"failed boundary stream readiness");
        check(warnings.str().find("boundary primary csv:")!=std::string::npos && warnings.str().find("boundary energy csv:")!=std::string::npos
            && warnings.str().find(implicit?"boundary flux csv:":"flux report:")!=std::string::npos,"boundary warning format");
    }
}
void checkStepDiagnostics() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    std::vector<std::array<double,4>> wall{{2.,3.,0.,1.},{1.,0.,0.,1.}};
    std::vector<double> embedded{6.,0.};
    std::unordered_map<std::string,double> reactions{{"elastic",5.},{"ionization",3.}};
    PicProcessSummaryInput input{-3.,1.,6.,2.,8.,6.,8.,.5,2,wall,embedded,reactions};
    auto summary=summarizePicProcesses(input);
    check(summary.displacementFraction==.25 && summary.dominantSource=="drive" && summary.dominantSourceFraction==.75,"process current/source fractions");
    check(std::string(summary.dominantBoundary)=="embedded" && summary.dominantBoundaryPowerW==12.
        && summary.dominantBoundaryFraction==3./7.,"summed species boundary selection/window");
    check(summary.dominantReactionChannel=="elastic" && summary.dominantReactionPowerW==10.
        && summary.dominantReactionFraction==.625,"reaction selection/window");
    check(wall[0][0]==2. && embedded[0]==6. && reactions.at("elastic")==5.,"diagnostics must not mutate accumulators");

    input.injectionEnergyJ=6.;wall={{{3.,3.,0.,0.}},{{0.,0.,0.,0.}}};embedded={0.,0.};
    reactions={{"first candidate",1.},{"second candidate",1.}};
    const auto firstReaction=reactions.begin()->first;summary=summarizePicProcesses(input);
    check(summary.dominantSource=="drive" && summary.dominantSourceFraction==.5,"source tie favors drive");
    check(std::string(summary.dominantBoundary)=="L" && summary.dominantReactionChannel==firstReaction,"boundary/map iteration tie order");
    reactions.clear();
    check(summary.dominantReactionChannel==firstReaction,"result must own selected reaction name");

    input.conductionCurrent=nan;input.displacementCurrent=std::numeric_limits<double>::infinity();
    input.driveEnergyJ=-1.;input.injectionEnergyJ=nan;input.timeWindowS=0.;
    input.domainAbsorptionJ=0.;input.embeddedAbsorptionJ=0.;input.reactionEnergyLossJ=0.;
    wall={{{0.,0.,0.,0.}},{{0.,0.,0.,0.}}};summary=summarizePicProcesses(input);
    check(summary.displacementFraction==0. && summary.dominantSource=="none" && summary.dominantSourceFraction==0.,"nonfinite current/negative source fallback");
    check(std::string(summary.dominantBoundary)=="L" && summary.dominantBoundaryPowerW==0. && summary.dominantBoundaryFraction==0.
        && summary.dominantReactionChannel=="none" && summary.dominantReactionPowerW==0. && summary.dominantReactionFraction==0.,"empty/zero window fallback");
    wall={{{1e16,1e16,0.,0.}},{{1.,2.,0.,0.}},{{1.,0.,0.,0.}}};embedded={0.,0.,0.};
    input.speciesCount=3;input.timeWindowS=1.;summary=summarizePicProcesses(input);
    check(std::string(summary.dominantBoundary)=="R","species energy accumulation order changed");

    const double upper=1.-1e-12;
    std::vector<double> alpha{0.,1e-12,std::nextafter(1e-12,1.),.5,std::nextafter(upper,0.),upper,1.,nan};
    const auto geometry=summarizePicGeometry(true,alpha,static_cast<int>(alpha.size()),nullptr,1);
    check(geometry.activeCutCellCount==3,"strict cut-cell thresholds or nonfinite fallback");
    check(geometry.movingRegionCount==0 && geometry.movingLevelsetCount==0 && geometry.movingRegionTotalAreaM2==0.
        && geometry.movingRegionMeanCos==1. && geometry.movingRegionMeanEnergyEV==0.,"absent moving geometry defaults");
    check(summarizePicGeometry(false,alpha,8,nullptr,1).activeCutCellCount==0
        && summarizePicGeometry(true,alpha,7,nullptr,1).activeCutCellCount==0,"cut metrics disabled/mismatched size");
    DictionaryGeometry dictionary;MovingEBConfig config;config.enabled=true;
    MovingEBRegion region;region.name="first";region.motionMode="rigidShift";region.levelSetIndex=3;region.startStep=2;region.endStep=4;
    config.regions.push_back(region);region.name="second";config.regions.push_back(region);
    region.name="fixed";region.motionMode="fixed";region.levelSetIndex=4;config.regions.push_back(region);
    MovingEBController controller(dictionary,config);
    const auto active=summarizePicGeometry(true,alpha,8,&controller,2);
    check(active.movingRegionCount==3 && active.movingLevelsetCount==1,"active region count/deduplicated moving level sets");
    check(active.movingRegionTotalAreaM2==0. && active.movingRegionMeanEnergyEV==0. && active.movingRegionMeanCos==1.
        && active.movingRegionMeanAbsSpeedMps==0.,"zero-area active-region fallback");
    check(summarizePicGeometry(true,alpha,8,&controller,1).movingRegionCount==0
        && summarizePicGeometry(true,alpha,8,&controller,4).movingRegionCount==3
        && summarizePicGeometry(true,alpha,8,&controller,5).movingRegionCount==0,"inclusive geometry activity window");
}
void checkDetailSeries() {
    using Series=PicStepSeries;
    using Result=ExplicitPIC::Diagnostics::PicDiagnostics::Result;
    Result row;row.nActive=13;row.rG_L2=.125;row.energyW=7.;row.Qgrid=8.;row.Qmacro=9.;
    row.speciesNames={"Ar+","AR+","e-","Ar*","???"};
    row.speciesNmean_1pm3={11.,12.,13.,14.,15.};
    row.speciesMacroCount={21,22,9007199254740993LL,24,25};
    row.speciesPhysCount={31.,32.,33.,34.,35.};
    row.macro_total=9007199254741085LL;row.phys_total=165.;
    row.electron_mean_energy_eV=4.5;row.Te_e_mean_eV=3.;
    UnitsScales scales;scales.T0_s=2.;scales.Phi0_V=3.;
    const auto units=UnitsContext::MakeSI(scales);
    for(bool si:{false,true}) {
        const std::string path=si?"detail/si.csv":"detail/code.csv";
        auto write=[&](Series& series,int step,const Result& values,int resume) {
            series.writeDetail(path,step,.25,values,&units,si,resume);
        };
        {
            Series series(Series::Mode::Explicit);
            for(int step:{0,1,2,4})write(series,step,row,0);
        } // A session owns and closes even a partially filled buffer.
        auto original=csvRows(path);check(original.size()==5,"detail destruction flushes all rows");
        const auto header=original[0];const std::size_t prefix=si?24:22;
        check(header.size()==prefix+19 && header[prefix]=="n_arplus_mean_1pm3"
            && header[prefix+1]=="n_arplus_2_mean_1pm3" && header[prefix+2]=="n_eminus_mean_1pm3"
            && header[prefix+3]=="n_arstar_mean_1pm3" && header[prefix+4]=="n_unnamed_mean_1pm3",
            "detail keeps legacy charged/excited species column keys and order");
        for(std::size_t i=1;i<original.size();++i)check(original[i].size()==header.size(),"fresh detail row width");
        check(original[1][prefix+13]=="9007199254740993","detail integer counts must not round through double");
        if(si)check(original[1][2]=="0.5" && original[1][19]=="21","detail time and energy SI conversion");
        else check(original[1][1]=="0.25" && original[1][17]=="7","detail code units layout");
        {
            Series resumed(Series::Mode::Explicit);
            write(resumed,2,row,2);
            auto reordered=row;
            reordered.speciesNames={"e-","Ar+"};reordered.speciesNmean_1pm3={103.,101.};
            reordered.speciesMacroCount={203,201};reordered.speciesPhysCount={303.,301.};
            write(resumed,3,reordered,2);
        }
        const auto restored=csvRows(path);
        check(restored.size()==5 && restored[0]==header && restored[1]==original[1] && restored[2]==original[2]
            && restored[3]==original[3],"resume restores species schema and replaces rows at/after the checkpoint step");
        const auto& last=restored[4];check(last.size()==header.size(),"resumed detail row keeps every species column");
        check(last[prefix]=="101" && last[prefix+1]=="0" && last[prefix+2]=="103"
            && last[prefix+3]=="0" && last[prefix+4]=="0" && last[prefix+9]=="201"
            && last[prefix+11]=="0" && last[prefix+13]=="203", "later rows map raw names and zero missing species");
        const auto before=bytes(path);
        for(int mismatch:{0,1}) {
            auto changed=row;if(mismatch==0)changed.speciesNames[0]="Ne+";
            bool rejected=false;
            try {
                Series incompatible(Series::Mode::Explicit);
                incompatible.writeDetail(path,2,.25,changed,&units,mismatch==0?si:!si,2);
            } catch(const std::runtime_error&) {rejected=true;}
            check(rejected && bytes(path)==before,"incompatible species/units header rejected before trimming existing history");
        }
        {
            Series fresh(Series::Mode::Explicit);write(fresh,0,row,0);
        }
        check(csvRows(path).size()==2,"second fresh session at same path must discard old run");
        {
            Series buffered(Series::Mode::Explicit);
            for(int step=0;step<64;++step)write(buffered,step,row,0);
            check(csvRows(path).size()==65,"detail keeps its 64-row flush boundary");
            bool rejected=false;
            try {buffered.writeDetail(path+"_other",64,.25,row,&units,si,0);}
            catch(const std::runtime_error&) {rejected=true;}
            check(rejected && !std::filesystem::exists(path+"_other"),"a detail writer cannot silently switch output files");
        }
    }
    {
        Result empty;Series missing(Series::Mode::Explicit);
        missing.writeDetail("detail/new_resume.csv",5,1.,empty,nullptr,true,5);
    }
    auto empty=csvRows("detail/new_resume.csv");
    check(empty.size()==2 && empty[0].size()==26 && empty[1].size()==26 && empty[0][1]=="time",
        "missing history and zero species use a complete code-unit header when units are absent");
    {
        auto shortRow=row;shortRow.speciesPhysCount.resize(2);
        Series partial(Series::Mode::Explicit);partial.writeDetail("detail/partial.csv",0,0.,shortRow);
    }
    check(csvRows("detail/partial.csv")[0].size()==32,"initial species layout uses all available vector lengths");
    std::filesystem::create_directories("detail/blocked/child");
    Series failed(Series::Mode::Explicit);failed.writeDetail("detail/blocked",0,0.,row);
    check(std::filesystem::is_directory("detail/blocked/child"),"optional detail open failure preserves original no-write policy");
}
void checkStepSeries() {
    using Series=PicStepSeries;
    const long long large=9007199254740993LL;
    const auto count=static_cast<std::size_t>(large);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    for(auto mode:{Series::Mode::Explicit,Series::Mode::Implicit}) {
        const bool implicit=mode==Series::Mode::Implicit;
        const std::string dir=implicit?"step_series/implicit/":"step_series/explicit/";
        auto open=[&](Series& output,int resume) {
            output.openOverview(dir+"overview.csv",resume);
            output.openProcess(dir+"process.csv",resume);
            output.openGeometry(dir+"geometry.csv",resume);
            output.openPeriodic(dir+"periodic.csv",resume);
        };
        auto write=[&](Series& output,int step) {
            if(implicit)output.writeOverview(step,.25,Series::ImplicitOverview{1.,2.,
                {7.,8.,9.,10.,11.,12.,13.,nan,15.,large,17.,18.}});
            else output.writeOverview(step,.25,Series::ExplicitOverview{1.,2.,3.,4.,5.,6.,
                {7.,8.,9.,10.,11.,12.,13.,nan,15.,3.5,17.,18.}});
            std::string source="injection",reaction="raw-channel";
            output.writeProcess(step,.25,{1.,2.,3.,4.,5.,6.,7.,.125,source,.5,"embedded",12.,.75,reaction,15.,.875});
            source="changed";reaction="changed"; // The synchronous writer must not retain the views.
            output.writeGeometry(step,.25,{1,count,3.,4u,5u,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.});
            output.writePeriodic(step,.25,.125,2.,true,false,step==2);
        };
        {
            Series output(mode);
            check(!output.overviewReady() && !output.processReady() && !output.geometryReady(),"unopened step streams");
            check(!output.periodicReady(),"unopened periodic stream");
            open(output,0);
            check(output.overviewReady() && output.processReady() && output.geometryReady(),"opened step streams");
            check(output.periodicReady(),"opened periodic stream");
            if(implicit)for(const auto* file:{"overview","process","geometry","periodic"})
                check(csvRows(dir+file+".csv").size()==1,"implicit step headers must flush");
            write(output,1);write(output,2);
            if(!implicit)check(csvRows(dir+"overview.csv").size()==3,"explicit overview rows must flush");
            check(csvRows(dir+"process.csv").size()==3 && csvRows(dir+"geometry.csv").size()==3,"summary rows must flush");
            check(csvRows(dir+"periodic.csv").size()==3,"periodic rows flush on implicit report/each explicit step");
        }
        const auto overview=csvRows(dir+"overview.csv"),process=csvRows(dir+"process.csv"),geometry=csvRows(dir+"geometry.csv");
        const auto periodic=csvRows(dir+"periodic.csv");
        check(periodic[0]==std::vector<std::string>{"step","time_s","phiJumpX_V","phiJumpY_V","periodicX","periodicY"}
            && periodic[1]==std::vector<std::string>{"1","0.2500000000000000","0.1250000000000000","2.0000000000000000","1","0"},
            "periodic CSV layout and precision");
        check(overview.size()==3 && overview[0].size()==(implicit?16:20),"distinct overview schema widths");
        check(overview[1].size()==overview[0].size() && overview[1][1]=="0.2500000000000000","overview row width/precision");
        const int first=implicit?4:8;
        check(overview[0][first]=="field_energy_J" && overview[0][first+9]=="macro_total","overview layout");
        for(int i=2;i<first;++i)check(std::stod(overview[1][i])==i-1.,"residual column order");
        for(int i=0;i<12;++i) {
            const auto& value=overview[1][first+i];
            if(i==7)check(std::isnan(std::stod(value)),"undefined overview current");
            else if(i==9)check(value==(implicit?std::to_string(large):"3.5000000000000000"),"overview macro count type/precision");
            else check(std::stod(value)==i+7.,"overview value order");
        }
        check(process[0].size()==18 && process[1].size()==18 && process[0][10]=="dominant_source","process schema");
        const std::vector<std::string> expectedProcess{"1","0.2500000000000000","1.0000000000000000","2.0000000000000000",
            "3.0000000000000000","4.0000000000000000","5.0000000000000000","6.0000000000000000","7.0000000000000000",
            "0.1250000000000000","injection","0.5000000000000000","embedded","12.0000000000000000",
            "0.7500000000000000","raw-channel","15.0000000000000000","0.8750000000000000"};
        check(process[1]==expectedProcess,"process values, names or formatting changed");
        check(geometry[0].size()==20 && geometry[1].size()==20 && geometry[0][3]=="active_cut_cell_count","geometry schema");
        check(geometry[1][2]=="1" && geometry[1][3]==std::to_string(count) && geometry[1][5]=="4" && geometry[1][6]=="5","geometry integer counts");
        for(int i=4;i<20;++i)check(std::stod(geometry[1][i])==i-1.,"geometry value order");
        {
            Series resumed(mode);open(resumed,2);write(resumed,3);resumed.flushPeriodic();
            check(csvRows(dir+"periodic.csv").size()==3,"explicit end-of-step periodic flush");
        }
        for(const auto* file:{"overview","process","geometry","periodic"}) {
            const auto rows=csvRows(dir+file+".csv");
            check(rows.size()==3 && rows[1][0]=="1" && rows[2][0]=="3","step-series resume pruning/header");
            std::filesystem::create_directories(dir+file+"_blocked.csv/child");
        }
        Series failed(mode);std::ostringstream warnings;auto* old=std::cerr.rdbuf(warnings.rdbuf());
        failed.openOverview(dir+"overview_blocked.csv",0);failed.openProcess(dir+"process_blocked.csv",0);failed.openGeometry(dir+"geometry_blocked.csv",0);
        failed.openPeriodic(dir+"periodic_blocked.csv",0);
        std::cerr.rdbuf(old);
        check(!failed.overviewReady() && !failed.processReady() && !failed.geometryReady(),"failed step streams");
        check(!failed.periodicReady(),"failed periodic stream");
        for(const auto* role:{"overview primary","process summary","geometry coupling","periodic diagnostics"})
            check(warnings.str().find(std::string(implicit?"[ImplicitPIC]":"[ExplicitPIC]")+" WARNING: failed to create "+role+" csv:")!=std::string::npos,"step-series warning");
    }
}
void checkParticleManagementSeries() {
    using Series=PicStepSeries;
    PicSpawnManagement::Stats spawn;spawn.in=7;spawn.kept=3;spawn.drop=4;spawn.physIn=14.;spawn.physKept=15.;
    ParticleManagementReport report;report.macroBefore=5;report.physBefore=15.;report.macroAfter=3;report.physAfter=15.;
    report.counts={2,1,1,0};report.weights={2.,3.,4.,5.,.25};
    for(auto mode:{Series::Mode::Explicit,Series::Mode::Implicit}) {
        const bool implicit=mode==Series::Mode::Implicit;const std::string file=implicit?"pm/implicit.csv":"pm/explicit.csv";
        {
            Series output(mode);check(!output.particleManagementReady(),"unopened population CSV");output.openParticleManagement(file,0);check(output.particleManagementReady(),"population CSV opens");
            if(implicit)check(csvRows(file).size()==1,"implicit population header flushes");
            output.writeSpawnManagement(1,.25,"mcc",spawn,.2);output.writeParticleManagement(1,.25,"electron",report);output.flushParticleManagement();
            auto rows=csvRows(file);check(rows.size()==3 && rows[0].size()==21 && rows[1].size()==21 && rows[2].size()==21,"population CSV column contract");
            check(rows[1][1]==(implicit?"0.2500000000000000":"0.25") && rows[2][1]==rows[1][1],"original population float formats");
            check(rows[1][2]=="spawn" && rows[1][3]=="mcc" && rows[1][12]=="7" && rows[1][13]=="3" && rows[1][14]=="4","spawn column mapping");
            check(std::stod(rows[1][7])==15. && std::stod(rows[1][15])==.2,"spawn physical count and effective probability");
            check(rows[2][2]=="species" && rows[2][3]=="electron" && rows[2][8]=="2" && rows[2][9]=="1" && rows[2][10]=="1","resampling operation counts");
            check(std::stod(rows[2][16])==2. && std::stod(rows[2][17])==3. && std::stod(rows[2][18])==4. && std::stod(rows[2][19])==5. && std::stod(rows[2][20])==.25,"weight-tail column mapping");
            output.writeParticleManagement(3,.75,"electron",report);
        }
        check(csvRows(file).size()==4,"population destructor flush");
        {Series resumed(mode);resumed.openParticleManagement(file,3);resumed.writeParticleManagement(4,1.,"electron",report);}
        const auto rows=csvRows(file);check(rows.size()==4 && rows[1][0]=="1" && rows[2][0]=="1" && rows[3][0]=="4","population restart prunes all rows at and after resume step");
        const auto blocked=file+"_blocked";std::filesystem::create_directories(blocked+"/child");
        Series failed(mode);std::ostringstream warning;auto* previous=std::cerr.rdbuf(warning.rdbuf());failed.openParticleManagement(blocked,0);std::cerr.rdbuf(previous);
        check(!failed.particleManagementReady() && warning.str().find(implicit?"[ImplicitPIC] WARNING: failed to create particle management csv:":"[ExplicitPIC] WARNING: failed to create particle management report:")!=std::string::npos,"population failure warning and gate");
    }
}
void checkCurrentBalanceSeries() {
    using Series=PicStepSeries;
    PicExplicitCurrentBalance explicitRow;
    explicitRow.ionInjectionA=1.;explicitRow.electronInjectionA=2.;explicitRow.injectionMagnitudeA=3.;
    explicitRow.domainMagnitudeA=4.;explicitRow.embeddedMagnitudeA=5.;explicitRow.magnitudeResidualA=6.;
    explicitRow.magnitudeRelative=7.;explicitRow.injectionSignedA=8.;explicitRow.domainSignedA=9.;
    explicitRow.embeddedSignedA=10.;explicitRow.signedResidualA=11.;explicitRow.signedRelative=12.;
    PicStepParticleBalance balance(10.);balance.addSource(4.);balance.addDomain(2.);balance.addEmbedded(-1.);balance.finish(14.);
    for(auto mode:{Series::Mode::Explicit,Series::Mode::Implicit}) {
        const bool implicit=mode==Series::Mode::Implicit;
        const std::string file=implicit?"charge/implicit.csv":"charge/explicit.csv";
        const auto write=[&](Series& output,int step) {
            if(implicit) output.writeCurrentBalance(step,.5,balance.current());
            else output.writeCurrentBalance(step,.5,explicitRow);
        };
        {
            Series output(mode);check(!output.currentBalanceReady(),"unopened current balance gate");
            output.openCurrentBalance(file,0);check(output.currentBalanceReady(),"current balance stream opens");
            if(implicit) check(csvRows(file).size()==1,"implicit header flushes on open");
            write(output,3);if(implicit) output.flushCurrentBalance();
            auto rows=csvRows(file);check(rows.size()==2 && rows[0].size()==rows[1].size(),"charge header width and legacy flush points");
            check(rows[1][1]==(implicit?"0.5000000000000000":"0.5"),"implicit fixed versus explicit defaultfloat");
            if(implicit) {
                check(rows[0][2]=="Qsource_injected_C" && rows[0].back()=="dQphys_rel" && rows[0].size()==10,"implicit charge columns");
                const std::array<double,8> expected{4.,2.,-1.,10.,14.,4.,1.,1./14.};
                for(std::size_t i=0;i<expected.size();++i) check(std::abs(std::stod(rows[1][i+2])-expected[i])<1e-15,"implicit charge row order");
            } else {
                check(rows[0][2]=="Ii_inj_A" && rows[0].back()=="deltaI_signed_rel" && rows[0].size()==14,"explicit current columns");
                for(int i=1;i<=12;++i) check(std::stod(rows[1][i+1])==i,"explicit current row order");
            }
            write(output,6);
        }
        check(csvRows(file).size()==3,"charge destructor flush");
        {Series resumed(mode);resumed.openCurrentBalance(file,6);write(resumed,7);}
        auto rows=csvRows(file);check(rows.size()==3 && rows[1][0]=="3" && rows[2][0]=="7","inclusive charge CSV restart pruning");
        const auto blocked=file+"_blocked";std::filesystem::create_directories(blocked+"/child");
        Series failed(mode);std::ostringstream warnings;auto* previous=std::cerr.rdbuf(warnings.rdbuf());
        failed.openCurrentBalance(blocked,0);std::cerr.rdbuf(previous);
        check(!failed.currentBalanceReady() && warnings.str().find(implicit?"[ImplicitPIC] WARNING: failed to create current balance csv:":"[ExplicitPIC] WARNING: failed to create current-balance csv:")!=std::string::npos,"failed charge stream warning/gate");
    }
}
void checkPowerSeries() {
    using Series=PicStepSeries;
    const Series::Power row{1.,2.,3.,4.,5.,6.,7.,8.,{9.,10.,11.},{12.,13.,14.},15.,17.,18.,19.,20.,
        {22.,31.,32.,16.,21.,23.},{27.,28.,29.,30.,24.,25.,26.},33.,34.};
    const std::array<std::string,3> expected{
        "7,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,24,25,26,27,28,29,30,31,32,33,34",
        "7,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,24,25,26,0,28,29,30,31,32,0,0",
        "7,0.5,1,2,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,28,29,30,31,32,33,34"};
    for(int kind=0;kind<3;++kind) {
        const auto mode=kind==2?Series::Mode::Implicit:Series::Mode::Explicit;
        const std::string file="power/"+std::to_string(kind)+".csv";
        {
            Series output(mode);check(!output.powerReady(),"unopened power stream must not enable sampling");
            output.openPower(file,0,kind!=1);check(output.powerReady(),"power report open");
            output.writePower(7,.5,row);output.flushPower();
            std::ifstream stream(file);std::string header,line;std::getline(stream,header);std::getline(stream,line);
            check(line==expected[kind],"power column order, mode layout, defaultfloat precision or code-unit zeros changed");
            check(header.find(kind==1?"step,time,":"step,time_s,")==0,"power header time unit");
            check(header.find(kind==1?"field_energy,field_dWdt":"field_energy_J,field_dWdt_W")!=std::string::npos,"power field unit header");
            check((header.find("total_dWdt_W")!=std::string::npos)==(kind==2),"implicit total derivative column");
            output.writePower(9,.75,row);
        }
        auto rows=csvRows(file);check(rows.size()==3 && rows[0].size()==rows[1].size(),"power header width and destructor flush");
        {
            Series resumed(mode);resumed.openPower(file,9,kind!=1);resumed.writePower(10,1.,row);
        }
        rows=csvRows(file);check(rows.size()==3 && rows[1][0]=="7" && rows[2][0]=="10","power inclusive restart trimming without duplicate header");
        const auto blocked=file+"_blocked";std::filesystem::create_directories(blocked+"/child");
        Series failed(mode);std::ostringstream warnings;auto* previous=std::cerr.rdbuf(warnings.rdbuf());
        failed.openPower(blocked,0,kind!=1);std::cerr.rdbuf(previous);
        check(!failed.powerReady() && warnings.str().find(kind==2?"[ImplicitPIC] WARNING: failed to create power report:":"[ExplicitPIC] WARNING: failed to create power step report:")!=std::string::npos,
            "failed power stream readiness/warning");
    }
}
void checkWaveformSeries() {
    using Series=PicStepSeries;
    for(auto mode:{Series::Mode::Explicit,Series::Mode::Implicit}) {
        const bool implicit=mode==Series::Mode::Implicit;
        const std::string file=implicit?"waveform/implicit.csv":"waveform/explicit.csv";
        {
            Series output(mode);check(!output.waveformDue(2),"disabled waveform does not write");
            output.openWaveform(file,0,2);
            check(!output.waveformDue(1) && output.waveformDue(2) && !output.waveformDue(3),"waveform CSV cadence");
            output.writeWaveform(2,.25,2.5,3.,-1.);
            output.writeWaveform(4,.5,5.,6.,-2.);
        }
        auto rows=csvRows(file);
        check(rows.size()==3 && rows[0]==std::vector<std::string>{"step","t_si_s","V_drive_V","J_cond_Am2","J_disp_Am2","J_total_Am2"},
            "legacy waveform header/destructor flush");
        check(rows[1][1]=="2.50000000000000000e-01" && rows[1][2]=="2.50000000000000000e+00"
            && std::stod(rows[1][3])==3. && std::stod(rows[1][4])==-1. && std::stod(rows[1][5])==2.,"waveform precision and current sum");
        {
            Series resumed(mode);resumed.openWaveform(file,4,2);resumed.writeWaveform(6,.75,7.,8.,-3.);
        }
        rows=csvRows(file);check(rows.size()==3 && rows[1][0]=="2" && rows[2][0]=="6","waveform inclusive resume pruning");
        const auto blocked=file+"_blocked";std::filesystem::create_directories(blocked+"/child");
        Series failed(mode);std::ostringstream warnings;auto* previous=std::cerr.rdbuf(warnings.rdbuf());
        failed.openWaveform(blocked,0,2);std::cerr.rdbuf(previous);
        check(!failed.waveformDue(2) && warnings.str().find(std::string(implicit?"[ImplicitPIC]":"[ExplicitPIC]")+" WARNING: failed to create current waveform csv:")!=std::string::npos,
            "waveform failed open warning and disabled writes");
    }
}
void checkExplicitOutputStage() {
    struct Location {std::filesystem::path old=std::filesystem::current_path();~Location(){std::filesystem::current_path(old);}} location;
    struct Logs {std::ostringstream out,err;std::streambuf* o=std::cout.rdbuf(out.rdbuf());std::streambuf* e=std::cerr.rdbuf(err.rdbuf());
        ~Logs(){std::cout.rdbuf(o);std::cerr.rdbuf(e);}} logs;
    const auto root=location.old/"explicit_output_stage";std::filesystem::create_directories(root);
    struct Display {int elements;std::vector<std::string> names;std::vector<std::vector<double>> nodes;};
    const auto readDisplay=[](int step) {
        std::ostringstream tag;tag<<std::setw(6)<<std::setfill('0')<<step;
        std::ifstream f("output/fields/fields_"+tag.str()+".dat");std::string line;Display display{};
        std::getline(f,line);check(line=="Title= \"fields_"+tag.str()+"\"","explicit output title and unshifted step");
        std::getline(f,line);const std::regex quoted("\"([^\"]+)\"");
        for(auto it=std::sregex_iterator(line.begin(),line.end(),quoted);it!=std::sregex_iterator();++it)display.names.push_back((*it)[1]);
        std::getline(f,line);std::smatch m;check(std::regex_match(line,m,std::regex("ZONE N = ([0-9]+), E = ([0-9]+), F = FEPOINT, ET = QUADRILATERAL")),"explicit stage zone");
        display.elements=std::stoi(m[2]);display.nodes.resize(std::stoi(m[1]),std::vector<double>(display.names.size()));
        for(auto& row:display.nodes)for(auto& value:row)check(bool(f>>value),"explicit stage complete values");
        for(int i=0;i<display.elements*4;++i){int label;check(bool(f>>label) && label>0 && label<=static_cast<int>(display.nodes.size()),"explicit connectivity uses prepared vertices");}
        return display;
    };
    for(bool profileSI:{false,true}) {
        const auto dir=root/(profileSI?"profile_si":"profile_code");std::filesystem::create_directory(dir);std::filesystem::current_path(dir);
        DictionaryMesh settings;settings.Nx=settings.Ny=2;settings.maxCellNum=64;settings.maxInitialIterations=0;
        QuadtreeMesh mesh(settings);FaceSegmentGraph graph;graph.build(mesh);
        PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);
        UnitsScales scales;scales.L0_m=2.;scales.T0_s=3.;scales.Phi0_V=6.;
        const auto units=std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));const auto norm=NormalizationContext::MakeFromUnits(units,{});
        CellField<double> phi(mesh,2.),rho(mesh,5.);CellField<Vector2D> electric(mesh,{4.,-2.});
        std::vector<ParticleType> species;species.emplace_back("electron",2.*PhysicalConstants::eV_J,-1.);species.emplace_back("Ar+",1.,1.);
        const std::vector<double> weights{2.,3.};
        for(int ci=0;ci<4;++ci){Particle p;p.position=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();p.inCellL=mesh.findCellLabel(p.position);
            p.velocity={1.,2.};p.vz=3.;p.w=2.;p.trackId=100+ci;p.rngState=200+ci;species[0].particle.push_back(p);}species[0].number=4;
        PicCurrentWeights::Policy coupling(mesh,DepositionConfig{},PicCurrentWeights::Policy::Mode::Explicit,false,false);
        std::vector<std::string> fields{"E","phi","rho","psi","n","te","electron_mean_energy","w_mean","w_sum","w_max","n_macro"};
        PicFieldOutput::ExplicitInput input{mesh,graph,embedded,phi,rho,electric,species,weights,units.get(),&norm,coupling,fields};
        {
            PicProfiles profiles(PicFieldMode::Explicit,"output/fields",0,2,{true,2,true,2});
            for(int step:{0,1,2}) {
                phi.initialize(mesh,2.+step);PicFieldOutput::writeExplicit(step,.125,input,profiles,profileSI);
                const auto display=readDisplay(step);
                check(display.names==std::vector<std::string>{"X","Y","phi_SI_V","rho_SI_Cpm3","psi_SI_m","n_electron_SI_1pm3","n_arplus_SI_1pm3",
                    "electron_mean_energy_eV","Te_e_eV","w_mean","w_sum","w_max","n_macro","Ex_SI_Vpm","Ey_SI_Vpm"},"explicit stage retains fixed order and species naming");
                check(display.elements==4 && display.nodes.size()==9,"explicit vertex layout");
                for(const auto& row:display.nodes)check(row[0]>=0. && row[0]<=2. && row[1]>=0. && row[1]<=2.
                    && row[2]==6.*(2.+step) && row[3]==6.64064e-11 && row[4]==0. && row[5]==2. && row[6]==0.
                    && row[7]==6.22222 && row[8]==4.14815 && row[9]==2. && row[10]==2. && row[11]==2. && row[12]==1.
                    && row[13]==12. && row[14]==-6.,"analytic explicit SI fields, 3D energy and particle statistics at original stream precision");
            }
            check(!std::filesystem::exists("output/fields/profiles_avg.csv"),"field output must not advance solve-stage average window");
            DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
            graph.build(mesh);coupling.refreshAfterMeshChange();phi.initialize(mesh,6.);rho.initialize(mesh,5.);electric.initialize(mesh,{4.,-2.});
            for(auto& p:species[0].particle)p.inCellL=mesh.findCellLabel(p.position);
            PicFieldOutput::writeExplicit(3,.125,input,profiles,profileSI);
            const auto adapted=readDisplay(3);check(adapted.elements==16 && adapted.nodes.size()==25,"explicit stage regenerates coordinates and connectivity after real AMR");
            for(const auto& row:adapted.nodes)check(row[0]>=0. && row[0]<=2. && row[1]>=0. && row[1]<=2. && row[2]==36.,"post-AMR SI coordinates and current potential");
            for(bool missingUnits:{false,true}) {
                auto invalid=input;if(missingUnits)invalid.units=nullptr;else invalid.normalization=nullptr;
                bool failed=false;try{PicFieldOutput::writeExplicit(6,.125,invalid,profiles,profileSI);}catch(const std::runtime_error&){failed=true;}
                check(failed && !std::filesystem::exists("output/fields/fields_000006.dat"),"missing SI contexts fail before profile/field output");
            }
            std::filesystem::create_directory("output/fields/fields_000004.dat");bool failed=false;
            try{PicFieldOutput::writeExplicit(4,.125,input,profiles,profileSI);}catch(const std::runtime_error&){failed=true;}
            check(failed && std::filesystem::is_directory("output/fields/fields_000004.dat"),"field-open failure propagates after profile sampling");
        }
        const auto profile=csvRows("output/fields/profile1d.csv");check(profile.size()==7,"profiles retain step-zero, cadence and rows preceding field-open failure");
        const auto column=[&](const std::string& name){auto i=std::find(profile[0].begin(),profile[0].end(),name);check(i!=profile[0].end(),"profile column missing");return static_cast<std::size_t>(i-profile[0].begin());};
        check((std::find(profile[0].begin(),profile[0].end(),"time_SI_s")!=profile[0].end())==profileSI,"profile SI option does not change packed SI fields");
        for(std::size_t row=1;row<profile.size();++row){const int step=2*static_cast<int>((row-1)/2);check(profile[row][0]==std::to_string(step)
            && std::stod(profile[row][column("time_code")])==step*.125 && std::stod(profile[row][column("phi_code")])==2.+step,"profile samples original step and field before failed packed write");
            if(profileSI)check(std::stod(profile[row][column("time_SI_s")])==step*.375 && std::stod(profile[row][column("phi_V")])==6.*(2.+step),"profile SI time/potential conversion");}
        PicProfiles disabled(PicFieldMode::Explicit,"unused",0,1,{});
        fields.clear();PicFieldOutput::writeExplicit(7,.125,input,disabled,profileSI);const auto defaults=readDisplay(7);
        check(std::find(defaults.names.begin(),defaults.names.end(),"phi_SI_V")!=defaults.names.end()
            && std::find(defaults.names.begin(),defaults.names.end(),"n_macro")!=defaults.names.end(),"empty explicit selection retains defaults");
        fields={"unknown"};PicFieldOutput::writeExplicit(8,.125,input,disabled,profileSI);
        check(readDisplay(8).names==std::vector<std::string>{"X","Y","n_electron_SI_1pm3","n_arplus_SI_1pm3"},"unknown explicit selection still emits mandatory species densities");
        check(species[0].particle[0].trackId==100 && species[0].particle[0].rngState==200
            && species[0].particle[0].vz==3. && phi.value[0]==6. && rho.value[0]==5.,"explicit output preserves physical inputs");
    }
    check(logs.err.str().find("SI outputs require units.dict and normalization")!=std::string::npos,"original missing-context error retained");
}
void checkImplicitOutputStage() {
    using Queue=PicFieldOutputQueue;
    struct Location {std::filesystem::path old=std::filesystem::current_path();~Location(){std::filesystem::current_path(old);}} location;
    struct Logs {std::ostringstream out;std::streambuf* old=std::cout.rdbuf(out.rdbuf());~Logs(){std::cout.rdbuf(old);}} logs;
    const auto stageRoot=location.old/"implicit_output_stage";std::filesystem::create_directories(stageRoot);
    struct Display {int elements;std::vector<std::string> names;std::vector<std::vector<double>> nodes;};
    const auto readDisplay=[](const std::filesystem::path& path) {
        std::ifstream f(path);std::string line;Display display{};
        check(bool(std::getline(f,line)) && line.rfind("Title= \"fields/fields_",0)==0,"implicit output title");
        std::getline(f,line);const std::regex quoted("\"([^\"]+)\"");
        for(auto it=std::sregex_iterator(line.begin(),line.end(),quoted);it!=std::sregex_iterator();++it)display.names.push_back((*it)[1]);
        std::getline(f,line);std::smatch m;check(std::regex_match(line,m,std::regex("ZONE N = ([0-9]+), E = ([0-9]+), F = FEPOINT, ET = QUADRILATERAL")),"packed stage zone");
        display.elements=std::stoi(m[2]);display.nodes.resize(std::stoi(m[1]),std::vector<double>(display.names.size()));
        for(auto& row:display.nodes)for(auto& value:row)check(bool(f>>value),"packed stage complete values");
        for(int i=0;i<display.elements*4;++i){int label;check(bool(f>>label) && label>0 && label<=static_cast<int>(display.nodes.size()),"packed stage connectivity");}
        std::string extra;check(!(f>>extra),"packed stage trailing records");return display;
    };
    for(bool cells:{false,true})for(bool fast:{false,true})for(bool batch:{false,true}) {
        const auto dir=stageRoot/(std::string(cells?"cells":"vertices")+(fast?"_fast":"_csr")+(batch?"_batch":"_immediate"));
        std::filesystem::create_directory(dir);std::filesystem::current_path(dir);
        DictionaryMesh settings;settings.Nx=settings.Ny=2;settings.maxCellNum=64;settings.maxInitialIterations=0;
        QuadtreeMesh mesh(settings);FaceSegmentGraph graph;graph.build(mesh);
        PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Implicit);Solvers::Poisson::CutCellMetrics metrics;
        auto& mask=const_cast<std::vector<unsigned char>&>(metrics.activeMaskByCellIndex());mask.assign(4,1);
        auto& metricCells=const_cast<std::vector<Solvers::Poisson::CutCellMetrics::CellMetric>&>(metrics.cellMetricsByIndex());metricCells.resize(4);
        UnitsScales scales;scales.L0_m=2.;scales.Phi0_V=6.;
        const auto units=std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));
        const auto norm=NormalizationContext::MakeFromUnits(units,{});
        CellField<double> phi(mesh,2.),rho(mesh,5.);CellField<Vector2D> electric(mesh,{4.,-2.});
        std::vector<CellField<double>> amounts;amounts.emplace_back(mesh,7.);
        std::vector<double> mu(4,7.),kappa(4,3.),segments(graph.segments().size(),3.),weights{2.};
        std::vector<ParticleType> species;species.emplace_back("electron",2.*PhysicalConstants::eV_J,-1.);
        for(int i:{0,3}) {Particle p;p.position=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();p.inCellL=mesh.findCellLabel(p.position);
            p.velocity={1.,2.};p.vz=3.;p.w=2.;p.trackId=10+i;p.rngState=20+i;species[0].particle.push_back(p);}
        PicCurrentWeights::Policy coupling(mesh,DepositionConfig{},PicCurrentWeights::Policy::Mode::Implicit,false,false);
        std::vector<std::string> levelNames;
        Queue::Input input{graph,embedded,nullptr,levelNames,cells?&metrics:nullptr,false,phi,rho,electric,amounts,
            species,weights,&norm,coupling,mu,kappa,segments,false,false};
        ImplicitPicOverrides overrides;overrides.packedFlushEvery=3;overrides.outputFastDensity=fast;
        overrides.outputCellConstant=cells;overrides.outputMaskDiag=true;overrides.packedBatchWithAmr=batch;
        auto options=Queue::Options::fromImplicit(overrides);
        std::vector<std::string> fields{"muFace","E","phi","rho","n","te","psi_composite"};
        Queue queue(mesh,units.get(),fields,options);fields.clear();options.fastDensity=!fast;overrides.outputFastDensity=!fast;
        queue.configureBatching(true);check(queue.batchingEnabled()==batch,"AMR opt-in controls batching");
        bool badStep=false;try{queue.writeImplicit(0,input);}catch(const std::invalid_argument&){badStep=true;}
        check(badStep && !std::filesystem::exists("output"),"invalid completed step prepared a frame");
        const auto first=queue.writeImplicit(2,input);
        check(first.compute_s>=0. && queue.pending()==(batch?1u:0u),"implicit frame reaches original submission mode");
        check(std::filesystem::exists("output/meta/packed_mask_diag_2.csv"),"mask diagnostic precedes queued field write");
        phi.initialize(mesh,3.);queue.writeImplicit(4,input);
        phi.initialize(mesh,99.);rho.initialize(mesh,99.);electric.initialize(mesh,{99.,99.});amounts[0].initialize(mesh,99.);
        if(batch){bool rejected=false;try{queue.refreshCoordinates(mesh);}catch(const std::logic_error&){rejected=true;}
            check(rejected && queue.pending()==2 && !std::filesystem::exists("output/fields/fields_000002.dat"),"pending frames prevent coordinate changes before flush");}
        queue.flush();
        for(int step:{2,4}) {const auto display=readDisplay(std::string("output/fields/fields_00000")+std::to_string(step)+".dat");
            check(display.names==std::vector<std::string>{"X","Y","muFaceProxy_SI_Fpm","Ex_SI_Vpm","Ey_SI_Vpm","phi_SI_V","rho_SI_Cpm3","n_electron_SI_1pm3","Te_e_eV","psi_composite_SI_m"},"implicit selection order and aliases copied into owner");
            check(display.elements==4 && display.nodes.size()==(cells?16u:9u),"vertex/cell layout selected without changing physical mesh");
            bool differingDensity=false;
            for(const auto& row:display.nodes){check(row[0]>=0. && row[0]<=2. && row[1]>=0. && row[1]<=2.,"SI coordinates owned by output");
                // Packed fields retain the original default stream precision.
                check(row[2]==1.77084e-11 && row[3]==12. && row[4]==-6.
                    && row[5]==(step==2?12.:18.) && row[6]==6.64064e-11 && row[9]==0.,"queued physical sample, SI values and absent composite are preserved");
                if(fast)check(std::abs(row[7]-7.)<1e-10,"fast density consumes same-stage deposited amount");
                differingDensity|=std::abs(row[7]-7.)>1e-6;}
            check(fast || differingDensity,"CSR density uses particles instead of supplied fast amount");
        }
        check(species[0].particle[0].trackId==10 && species[0].particle[0].rngState==20,"output does not mutate particle identity or RNG");
        for(int step:{6,8,10})queue.writeImplicit(step,input);
        check(queue.pending()==0 && std::filesystem::exists("output/fields/fields_000010.dat"),"three prepared outputs trigger FIFO threshold flush");
        QuadtreeMesh other(settings);bool rejected=false;
        try{queue.refreshCoordinates(other);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"coordinate refresh rejects another mesh owner");
        DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
        check(mesh.cellIndexNum()==16,"output stage real AMR");graph.build(mesh);mask.assign(16,1);metricCells.resize(16);
        phi.initialize(mesh,2.);rho.initialize(mesh,5.);electric.initialize(mesh,{4.,-2.});amounts[0].initialize(mesh,7.);
        mu.assign(16,7.);kappa.assign(16,3.);segments.assign(graph.segments().size(),3.);coupling.refreshAfterMeshChange();
        for(auto& p:species[0].particle)p.inCellL=mesh.findCellLabel(p.position);
        queue.refreshCoordinates(mesh);queue.writeImplicit(12,input);queue.flush();
        const auto adapted=readDisplay("output/fields/fields_000012.dat");
        check(adapted.elements==16 && adapted.nodes.size()==(cells?64u:25u),"post-AMR frame and coordinate sizes refresh together");
        for(const auto& row:adapted.nodes)check(row[0]>=0. && row[0]<=2. && row[1]>=0. && row[1]<=2. && row[5]==12.,"post-AMR physical display");
        // Gates precede any dereference of optional physical state.
        auto gated=input;gated.normalization=nullptr;
        Queue noUnits(mesh,nullptr,{"phi"},Queue::Options{1});check(noUnits.writeImplicit(20,gated).compute_s==0.,"null units gate");
        Queue noSelection(mesh,units.get(),{},Queue::Options{1});check(noSelection.writeImplicit(20,gated).compute_s==0.,"empty selection gate");
        Queue unknown(mesh,units.get(),{"unknown"},Queue::Options{1,true,false,true,true});
        check(unknown.writeImplicit(21,input).compute_s==0. && std::filesystem::exists("output/meta/packed_mask_diag_21.csv")
            && !std::filesystem::exists("output/fields/fields_000021.dat"),"unknown selection preserves preceding mask output but submits no frame");
    }
    check(logs.out.str().find("batching disabled (AMR active)")!=std::string::npos
        && logs.out.str().find("batching forced with AMR")!=std::string::npos,"output policy retains mode reports");
}
struct Frame {
    std::vector<std::array<double,4>> nodes;
    std::vector<std::array<int,4>> cells;
};
Frame read(const std::filesystem::path& path,const std::string& title) {
    std::ifstream file(path);std::string line;
    std::getline(file,line);check(line=="Title= \""+title+"\"","title/path convention changed");
    std::getline(file,line);check(line=="Variables= \"X\" \"Y\" \"phi\" \"rho\" ","variable order changed");
    std::getline(file,line);std::smatch match;
    check(std::regex_match(line,match,std::regex("ZONE N = ([0-9]+), E = ([0-9]+), F = FEPOINT, ET = QUADRILATERAL")),"legacy zone format");
    Frame frame;frame.nodes.resize(std::stoi(match[1]));frame.cells.resize(std::stoi(match[2]));
    for(auto& row:frame.nodes)for(auto& value:row)check(bool(file>>value),"truncated node row");
    for(auto& row:frame.cells)for(auto& value:row) {
        check(bool(file>>value),"truncated connectivity");
        check(value>0 && value<=static_cast<int>(frame.nodes.size()),"one-based connectivity outside output nodes");
    }
    std::string extra;check(!(file>>extra),"unexpected output records");return frame;
}
}
int main() try {
    const auto root=std::filesystem::current_path()/
        ("pic_field_output_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);std::filesystem::current_path(root);
    checkSpeciesSeries();
    checkImplicitParticleSummary();
    checkLiveParticleSamples();
    checkDiagnosticSample(false);
    checkDiagnosticSample(true);
    checkBoundarySeries();
    checkStepDiagnostics();
    checkBoundaryWork();
    checkDiagnosticWindow();
    checkEnergyState();
    checkChargeBalance();
    checkImplicitSpawnBalance();
    checkCurrentBalanceSeries();
    checkParticleManagementSeries();
    checkEmbeddedGeometry();
    checkDomainOwnership();
    checkGeometrySetup();
    checkPowerSeries();
    checkStepSeries();
    checkPerformanceSeries();
    checkImplicitConsole();
    checkExplicitConsole();
    checkAmrSeries();
    checkCutcellStateOutput();
    checkImplicitDiagnosticSeries();
    checkExplicitEmbeddedBand();
    checkExplicitReports();
    checkGeometryReports();
    checkOutputManifest();
    checkOutputDirectory();
    checkPhysicsDiagnostics();
    checkPhysicsDiagnosticStage();
    checkDetailSeries();
    checkWaveformSeries();
    checkImplicitOutputStage();
    checkExplicitOutputStage();
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    DictionaryMesh::RefineBox box;box.min={0.,.5};box.max={1.,1.};box.level=1;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()==10,"real coarse/fine test mesh");
    checkBoundaryDiagnostics(mesh);
    checkElectrodeCurrent(mesh);
    checkBoundaryFieldSampling(mesh);
    checkReconstruction(mesh);
    checkFieldSelection(mesh);
    checkFieldFrame();
    checkOutputCoupling(mesh);
    checkProfiles(mesh);
    checkGeometryFields(mesh);
    checkMaskDiagnostics(mesh);
    checkOutputPreparation(mesh);
    checkParticleFields(mesh);
    std::vector<double> x(mesh.vertIndexNum()),y(x.size());
    std::vector<VertField<double>> vertices;vertices.emplace_back(mesh,0.);vertices.emplace_back(mesh,0.);
    std::map<std::pair<double,double>,int> vertexIndices;
    for(int vi=0;vi<mesh.vertIndexNum();++vi) {
        const auto point=mesh.verts(mesh.vertIndexToLabel[vi]).coord();
        x[vi]=3.+2.*point.x;y[vi]=7.+4.*point.y;
        vertices[0].value[vi]=10.+x[vi]-2.*y[vi];vertices[1].value[vi]=-vertices[0].value[vi];
        vertexIndices[{x[vi],y[vi]}]=vi;
    }
    std::vector<CellField<double>> cells;cells.emplace_back(mesh,0.);cells.emplace_back(mesh,0.);
    std::vector<unsigned char> mask(mesh.cellIndexNum(),0);
    std::set<std::pair<double,double>> expectedCenters;
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const auto point=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord();
        const double px=3.+2.*point.x,py=7.+4.*point.y;
        cells[0].value[ci]=10.+px-2.*py;cells[1].value[ci]=-cells[0].value[ci];
        if(point.x<.5) {mask[ci]=1;expectedCenters.insert({px,py});}
    }
    check(expectedCenters.size()==5,"masked coarse/fine cells");
    const std::vector<std::string> names{"phi","rho"};
    auto write=[&](int layout,const std::vector<unsigned char>* chosen) {
        if(layout==0)PicFieldOutput::writeExplicitVertices(mesh,"explicit",x,y,names,vertices,chosen);
        else if(layout==1)PicFieldOutput::writeImplicitVertices(mesh,"nested/implicit",x,y,names,vertices,chosen);
        else PicFieldOutput::writeImplicitCells(mesh,"nested/cells",x,y,names,cells,chosen);
    };
    const std::array<std::string,3> titles{"explicit","nested/implicit","nested/cells"};
    const std::array<std::filesystem::path,3> paths{"output/fields/explicit.dat","output/nested/implicit.dat","output/nested/cells.dat"};
    for(int layout=0;layout<3;++layout) {
        write(layout,&mask);const auto frame=read(paths[layout],titles[layout]);
        check(frame.cells.size()==expectedCenters.size(),"mask not applied to output cells");
        if(layout==2)check(frame.nodes.size()==4*frame.cells.size(),"cell-constant vertices must remain separate");
        else {
            check(frame.nodes.size()<4*frame.cells.size(),"vertex output must share nodes");
            int previous=-1;std::pair<double,double> previousCoordinate{-1.,-1.};
            for(const auto& node:frame.nodes) {
                const int vi=vertexIndices.at({node[0],node[1]});
                check(node[2]==10.+node[0]-2.*node[1] && node[3]==-node[2],"vertex fields or supplied coordinates changed");
                if(layout==0)check(vi>previous,"explicit mesh-index vertex order");
                else check(std::make_pair(node[1],node[0])>previousCoordinate,"implicit coordinate vertex order");
                previous=vi;previousCoordinate={node[1],node[0]};
            }
        }
        std::set<std::pair<double,double>> actualCenters;
        for(const auto& connectivity:frame.cells) {
            double cx=0.,cy=0.;
            for(int label:connectivity) {cx+=frame.nodes[label-1][0]/4.;cy+=frame.nodes[label-1][1]/4.;}
            check(actualCenters.insert({cx,cy}).second,"duplicated output element");
            if(layout==2)for(int label:connectivity) {
                const auto& node=frame.nodes[label-1];
                check(node[2]==10.+cx-2.*cy && node[3]==-node[2],"cell constants were averaged or reordered");
            }
        }
        check(actualCenters==expectedCenters,"connectivity does not represent selected physical cells");
        write(layout,nullptr);const auto all=bytes(paths[layout]);
        const std::vector<unsigned char> wrongSize{0};write(layout,&wrongSize);
        check(bytes(paths[layout])==all,"legacy wrong-size mask must be ignored");
        const std::vector<unsigned char> empty(mesh.cellIndexNum(),0);write(layout,&empty);
        const auto none=read(paths[layout],titles[layout]);check(none.nodes.empty() && none.cells.empty(),"empty mask output");
    }
    PicFieldOutput::writeImplicitVertices(mesh,"queue/vertices",x,y,names,vertices,&mask);
    PicFieldOutput::writeImplicitCells(mesh,"queue/cells",x,y,names,cells,&mask);
    const auto vertexReference=bytes("output/queue/vertices.dat");
    const auto cellReference=bytes("output/queue/cells.dat");
    std::filesystem::rename("output/queue/vertices.dat","output/queue/vertices.reference");
    std::filesystem::rename("output/queue/cells.dat","output/queue/cells.reference");
    PicFieldOutputQueue queue(mesh,x,y,2);
    const auto savedX=x,savedY=y;x.assign(x.size(),1000.);y.assign(y.size(),-1000.);
    check(queue.pending()==0 && queue.flush()==0.,"empty queue flush");
    auto ownedVertices=vertices;auto ownedNames=names;auto changedMask=mask;
    check(queue.submitVertices("queue/vertices",std::move(ownedNames),std::move(ownedVertices),changedMask,true)==0.
        && queue.pending()==1 && !std::filesystem::exists("output/queue/vertices.dat"),"first frame must remain queued");
    changedMask.assign(mask.size(),0);
    check(queue.submitCells("queue/cells",std::vector<std::string>(names),
        std::vector<CellField<double>>(cells),mask,true)>=0. && queue.pending()==0,"mixed frame threshold flush");
    check(bytes("output/queue/vertices.dat")==vertexReference && bytes("output/queue/cells.dat")==cellReference,
        "queued field ownership, copied mask or FIFO serialization changed");
    x=savedX;y=savedY;
    check(queue.submitVertices("queue/remainder",std::vector<std::string>(names),
        std::vector<VertField<double>>(vertices),mask,true)==0. && queue.pending()==1,"remaining frame queued");
    check(!std::filesystem::exists("output/queue/remainder.dat") && queue.flush()>=0. && queue.pending()==0,
        "explicit final flush writes remaining frames");
    check(read("output/queue/remainder.dat","queue/remainder").cells.size()==expectedCenters.size(),"remaining frame content");
    auto immediateFields=vertices;auto immediateNames=names;
    check(queue.submitVertices("queue/immediate",std::move(immediateNames),std::move(immediateFields),mask,false)>=0.
        && queue.pending()==0 && immediateFields.size()==vertices.size() && immediateNames==names,
        "immediate path must not take ownership or queue a frame");
    check(read("output/queue/immediate.dat","queue/immediate").cells.size()==expectedCenters.size(),"immediate frame content");
    {
        PicFieldOutputQueue abandoned(mesh,x,y,2);
        abandoned.submitVertices("queue/abandoned",std::vector<std::string>(names),
            std::vector<VertField<double>>(vertices),mask,true);
        check(abandoned.pending()==1,"abandoned frame was queued");
    }
    check(!std::filesystem::exists("output/queue/abandoned.dat"),"queue destruction must not write pending frames");
    std::filesystem::create_directory("output/queue/failure_second.dat");
    {
        PicFieldOutputQueue failed(mesh,x,y,2);
        failed.submitVertices("queue/failure_first",std::vector<std::string>(names),
            std::vector<VertField<double>>(vertices),mask,true);
        bool writeFailed=false;
        try {failed.submitCells("queue/failure_second",std::vector<std::string>(names),
                std::vector<CellField<double>>(cells),mask,true);}
        catch(const std::runtime_error&) {writeFailed=true;}
        check(writeFailed && failed.pending()==2 && std::filesystem::exists("output/queue/failure_first.dat"),
            "failed flush must preserve queued jobs and earlier writes");
        std::ofstream marker("output/queue/failure_first.dat");marker<<"do not retry";
    }
    check(bytes("output/queue/failure_first.dat")=="do not retry","destruction retried partially failed IO");

    const auto before=bytes(paths[0]);bool rejected=false;
    try {PicFieldOutput::writeExplicitVertices(mesh,"explicit",{},y,names,vertices);}
    catch(const std::runtime_error&) {rejected=true;}
    check(rejected && bytes(paths[0])==before,"coordinate preflight must preserve existing file");
    rejected=false;
    try {PicFieldOutput::writeImplicitVertices(mesh,"nested/implicit",x,y,{},vertices);}
    catch(const std::runtime_error&) {rejected=true;}
    check(rejected,"vertex name/field mismatch rejection");
    rejected=false;
    try {PicFieldOutput::writeImplicitCells(mesh,"nested/cells",x,y,{},cells);}
    catch(const std::runtime_error&) {rejected=true;}
    check(rejected,"cell name/field mismatch rejection");
    rejected=false;
    try {PicFieldOutput::writeExplicitVertices(mesh,"missing_parent/file",x,y,names,vertices);}
    catch(const std::runtime_error&) {rejected=true;}
    check(rejected,"explicit writer must retain its original parent-directory policy");
    std::filesystem::create_directory("output/blocked.dat");rejected=false;
    try {PicFieldOutput::writeImplicitCells(mesh,"blocked",x,y,names,cells);}
    catch(const std::runtime_error&) {rejected=true;}
    check(rejected,"output-open failure must propagate");
    std::cout<<"PIC reconstruction, periodic vertices, field formats, AMR masks/connectivity, queued ownership, flush ordering and IO failures passed\n";
    return 0;
}
catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
