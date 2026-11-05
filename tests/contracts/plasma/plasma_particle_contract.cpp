#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Diagnostics/PicParticleTrajectories.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include "Modules/Plasma/Physics/PicChargeDeposition.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include "Modules/Plasma/Physics/PicParticleSources.h"
#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Modules/Plasma/ExplicitPIC/Field/MatchedGather.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Physics/PicRefinementIndicators.h"
#include "Modules/Plasma/Physics/PicMeshRefinement.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicPolarization.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicFieldSolve.h"
#include "Modules/Plasma/ImplicitPIC/Orbit/ImplicitPicOrbit.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include "Modules/Plasma/Boundary/PicBoundarySetup.h"
#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/ConsoleControl.h"
#include "Core/IO/DictionaryUnits.h"
#include "Modules/Plasma/IO/DictionaryNormalization.h"
#include "Modules/Plasma/IO/PicCaseInputLoader.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include <chrono>
#include <cstdlib>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
void check(bool valid,const char* reason) { if (!valid) throw std::runtime_error(reason); }
void near(double actual,double expected,const char* reason)
{ check(std::isfinite(actual) && std::abs(actual-expected)<1e-11,reason); }
struct BoundaryLogs {
    std::ostringstream output,errors;
    std::streambuf* oldOutput=std::cout.rdbuf(output.rdbuf());
    std::streambuf* oldErrors=std::cerr.rdbuf(errors.rdbuf());
    ~BoundaryLogs() { std::cout.rdbuf(oldOutput);std::cerr.rdbuf(oldErrors); }
};
void normalizationInputs(const std::filesystem::path& directory)
{
    BoundaryLogs logs;
    const auto close = [](double value,double reference,const char* message) {
        check(std::isfinite(value) && std::abs(value/reference-1.) < 5e-14,message);
    };
    const auto unitsPath=directory/"units.dict", normPath=directory/"normalization.dict";
    { std::ofstream out(unitsPath);out<<"Units\n{\nmode SI\nL0 2 cm\nT0 3 ms\nPhi0 7 V\nn0 1e12 1/m^3\nv0 99 m/s\nE0 88 V/m\nrho0 77 C/m^3\nwriteSI yes\n}\n"; }
    DictionaryUnits units;check(units.readDict(unitsPath.string()),"Core unit reader");
    check(units.context().mode()==UnitsMode::kSI && units.writeSIOutputs(),"unit mode/output flag");
    close(units.context().scales().L0_m,.02,"length reader SI conversion");
    const auto parse=[&](const std::string& text) {
        {std::ofstream out(normPath);out<<"Normalization\n{\ndimension 2\nslabThickness 2 mm\n"<<text<<"\n}\n";}
        DictionaryNormalization dict;check(dict.readDict(normPath.string()),"Plasma normalization reader");
        check(dict.hasDimension() && dict.hasSlabThickness() && dict.hasStrategy() && dict.hasAny(),"normalization presence flags");
        close(dict.config().slabThickness_m,.002,"normalization slab SI conversion");
        return dict.config();
    };
    auto config=parse("strategy manual");std::string warning;
    auto scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    check(warning.empty() && scales.T0_s==.003 && scales.Phi0_V==7. && scales.v0_mps==99.
        && scales.E0_Vpm==88. && scales.rho0_Cpm3==77.,"manual preserves explicit derived references");
    config=parse("strategy plasmaFrequency\nrefSpecies electron");
    scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    close(scales.T0_s,1.7725907105982083e-8,"electron plasma reference time");
    check(scales.Phi0_V==7. && scales.v0_mps==-1. && scales.E0_Vpm==-1. && scales.rho0_Cpm3==-1.,"plasma derived reset sentinels");
    const auto electronTime=scales.T0_s;
    config=parse("strategy wp\nrefSpecies Hplus");scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    close(scales.T0_s,7.595613604411247e-7,"proton plasma reference time");
    config=parse("strategy thermal\nTe0 3 eV\nrefSpecies electron");
    scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    close(scales.Phi0_V,3.,"thermal voltage");
    close(scales.T0_s,2.7533325513071663e-8,"thermal reference time");
    const auto finalized=std::make_shared<UnitsContext>(UnitsContext::MakeSI(scales));
    close(finalized->scales().v0_mps,scales.L0_m/scales.T0_s,"thermal velocity finalized by common units");
    check(finalized->scales().rho0_Cpm3==1.,"legacy reset charge scale fallback remains one");
    const auto norm=NormalizationContext::MakeFromUnits(finalized,config);
    close(norm.VolumeCodeToSI(5.),5.*.02*.02*.002,"normalized slab volume closure");
    config=parse("strategy te\nTe0 34813.55436465025 K\nrefCharge 2 e\nrefMass 4 amu");
    close(config.Te0_eV,3.,"Kelvin temperature reader");
    close(config.refMass_kg,4.*PhysicalConstants::amu_kg,"mass reader");
    scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    close(scales.Phi0_V,1.5,"explicit positive reference charge has priority");
    config.strategy=NormalizationStrategy::kPlasmaFrequency;config.refSpecies="unknown";warning.clear();
    resolvePlasmaReferenceScales(units.context(),config,warning);check(warning.empty(),"complete positive reference overrides unknown species");
    config.refCharge_C=-PhysicalConstants::qe_C;warning.clear();
    scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    close(scales.T0_s,electronTime,"legacy negative explicit reference falls back to species");
    check(warning.find("fallback to electron")!=std::string::npos,"unknown species warning");
    warning.clear();config.strategy=NormalizationStrategy::kThermal;config.Te0_eV=0.;
    scales=resolvePlasmaReferenceScales(units.context(),config,warning);
    check(scales.T0_s==.003 && scales.Phi0_V==7. && warning.find("Te0 missing")!=std::string::npos,"invalid thermal fallback/warning precedence");
    warning.clear();config.Te0_eV=3.;
    scales=resolvePlasmaReferenceScales(UnitsContext::MakeCodeUnits(),config,warning);
    check(scales.T0_s==1. && warning=="Units mode is CodeUnits; normalization strategy ignored.","code units ignore physical strategy");
    auto emptyDensity=units.context().scales();emptyDensity.n0_m3=0.;config.strategy=NormalizationStrategy::kPlasmaFrequency;
    scales=resolvePlasmaReferenceScales(UnitsContext::MakeSI(emptyDensity),config,warning);
    check(scales.T0_s==.003 && warning.find("invalid n0/q/m")!=std::string::npos,"invalid plasma density fallback");
    auto overflow=units.context().scales();overflow.n0_m3=std::numeric_limits<double>::infinity();
    scales=resolvePlasmaReferenceScales(UnitsContext::MakeSI(overflow),config,warning);
    check(scales.T0_s==.003 && warning.find("omega_p invalid")!=std::string::npos,"nonfinite plasma frequency fallback");
    check(units.context().scales().T0_s==.003 && config.refSpecies=="unknown","resolution does not mutate inputs");
    std::filesystem::remove(unitsPath);std::filesystem::remove(normPath);
}
void particleBoundaryOperations(QuadtreeMesh& mesh, ImmersedGeometry<QuadtreeMesh>& geometry,
    DictionaryGeometry& dictionary, PicEmbeddedGeometry& embedded)
{
    ParticleType original("energy",2.*PhysicalConstants::eV_J,1.);
    const std::array<double,5> speeds{0.,1.,2.,std::numeric_limits<double>::quiet_NaN(),3.};
    const std::array<double,5> weights{2.,0.,3.,5.,7.};
    for(std::size_t j=0;j<speeds.size();++j){Particle p;p.velocity={speeds[j],0.};p.w=weights[j];p.trackId=j+1;original.particle.push_back(p);}
    Particle energy;energy.velocity={3.,4.};energy.vz=12.;
    near(ParticleManager::kineticEnergyEV(energy,2.*PhysicalConstants::eV_J,2.),676.,"cutoff uses all three SI velocity components");
    near(ParticleManager::kineticEnergyEV(energy,0.,2.),0.,"invalid mass energy gate");
    std::vector<ParticleType> species{original,original};std::vector<SpeciesConfig> config(1);config[0].energyCutoff_eV=1.;
    std::vector<long long> macro{10};std::vector<double> physical{20.};
    check(ParticleManager::pruneByEnergyCutoff(species,config,{4.},1.,&macro,&physical)==2,"cutoff removes cold and nonfinite velocities");
    check(species[0].number==3 && species[1].particle.size()==5 && macro[0]==12,"cutoff uses common configuration prefix and accumulates counts");
    near(physical[0],27.,"cutoff physical loss preserves before/after subtraction and fallback weights");
    check(species[0].particle[0].trackId==2 && species[0].particle[1].trackId==3 && species[0].particle[2].trackId==5,"cutoff retains equal threshold and stable survivor order");
    auto disabled=original;disabled.number=-1;
    check(ParticleManager::pruneBelowEnergyCutoff(disabled,0.,1.)==0 && disabled.number==5,"disabled cutoff repairs number without filtering");
    check(ParticleManager::pruneBelowEnergyCutoff(disabled,1.,0.)==0,"invalid velocity scale does not filter");
    check(ParticleManager::pruneBelowEnergyCutoff(disabled,std::numeric_limits<double>::quiet_NaN(),1.)==0,"legacy nonfinite cutoff gate");

    ParticleType seed("cull",1.,1.);
    for(double x:{.1,.9,.2,.8,2.}){Particle p;p.position={x,.125};p.w=0.;p.trackId=seed.particle.size()+1;p.inCellL=mesh.findCellLabel({.125,.125});p.inBackgroundCellL=-1;seed.particle.push_back(p);}
    const int oldHost=seed.particle[0].inCellL;
    for(bool refresh:{false,true}) {
        species={seed,seed};species[1].number=99;macro={5,6};physical={7.,8.};const std::vector<unsigned char> mask{1};
        const auto report=ParticleManager::cullEmbedded(mesh,&geometry,embedded,species,{3.},refresh,&mask,&macro,&physical);
        check(report.removed==2 && report.kept==3 && species[0].number==3,"solid cull reports selected population");
        check(species[1].number==99 && species[1].particle.size()==5 && macro==std::vector<long long>({7,6}),"short touched mask leaves unselected species intact");
        near(physical[0],13.,"solid cull adds removed fallback weights");near(physical[1],8.,"unselected physical accounting unchanged");
        check(species[0].particle[0].trackId==2 && species[0].particle[1].trackId==4 && species[0].particle[2].trackId==5,"solid cull is stable and does not remove outside-domain particles");
        const auto& survivor=species[0].particle[0];const auto& outside=species[0].particle[2];
        check(refresh?(survivor.inCellL==mesh.findCellLabel(survivor.position) && survivor.inBackgroundCellL>=0 && outside.inCellL==-1 && outside.inBackgroundCellL==-1)
                     :(survivor.inCellL==oldHost && survivor.inBackgroundCellL==-1 && outside.inCellL==oldHost),"only implicit solid cull refreshes surviving hosts");
    }
    PicEmbeddedGeometry absent(PicEmbeddedGeometry::Mode::Explicit);species={seed};species[0].number=99;
    check(ParticleManager::cullEmbedded(mesh,&geometry,absent,species,{},true).removed==0 && species[0].number==99,"inactive geometry is a complete no-op");
    check(ParticleManager::cullEmbedded(mesh,nullptr,embedded,species,{},true).removed==0,"missing geometry is a complete no-op");

    using Emission=PicSecondaryEmission;using Side=ExplicitPIC::Pusher::DomainSide;
    ExplicitPIC::Pusher::EmbeddedImpact impact;impact.pos={.625,.125};impact.cellLabel=mesh.findCellLabel(impact.pos);
    Vector2D normal;double magnitude=0.;
    check(Emission::embeddedNormal(mesh,&geometry,dictionary,embedded,1,impact,normal,magnitude),"real embedded normal resolved");
    near(normal.x,1.,"normal points out of negative solid");near(normal.y,0.,"wall normal tangential component");
    check(!Emission::embeddedNormal(mesh,nullptr,dictionary,embedded,1,impact,normal,magnitude),"normal requires geometry");
    auto invalid=dictionary;invalid.boundaryById[1].levelSetIndex=99;
    check(!Emission::embeddedNormal(mesh,&geometry,invalid,embedded,1,impact,normal,magnitude),"boundary override validates level-set index");
    SEEConfig settings;settings.emittedEnergyModel="mono";settings.emittedEnergy_eV=3.;settings.angular="fixed";
    Emission emission(settings,Emission::Mode::Explicit);Emission::Engine rng(75),reference(75);Emission::Uniform uniform(0.,1.);
    const Emission::Kinematics kinematics{2.*PhysicalConstants::eV_J,7.,2.,.03};
    const Emission::DomainBox box{{0.,0.},{1.,1.},{.01,.02}};
    const std::array<Side,5> sides{Side::kLeft,Side::kRight,Side::kBottom,Side::kTop,Side::kUnknown};
    const std::array<Vector2D,5> expected{{{.03,.98},{.97,.98},{.01,.03},{.01,.97},{.01,.98}}};
    for(std::size_t j=0;j<sides.size();++j){ExplicitPIC::Pusher::DomainImpact domain;domain.side=sides[j];domain.particle.position={-2.,2.};
        const auto created=emission.spawnDomain(mesh,domain,box,kinematics,rng,uniform);check(created.has_value(),"domain emission clamps into a valid host");
        near(created->particle.position.x,expected[j].x,"domain emission x clamp");near(created->particle.position.y,expected[j].y,"domain emission y clamp");
        near(created->energyJ/PhysicalConstants::eV_J,3.,"emission carries sampled SI energy");near(created->particle.w,7.,"emission weight");
        check(created->particle.inCellL>=0 && created->particle.inBackgroundCellL>=0 && created->particle.rngState==0 && created->particle.trackId==0,"emission fills hosts and retains fresh particle metadata");}
    check(rng==reference,"mono fixed placement consumes no RNG");
    settings.emittedEnergyModel="maxwell";settings.emittedTe_eV=2.;settings.angular="cosine";
    for(auto mode:{Emission::Mode::Explicit,Emission::Mode::Implicit}) {
        Emission actual(settings,mode),expectedModel(settings,mode);Emission::Engine actualRng(63),expectedRng(63);Emission::Uniform actualUniform(0.,1.),expectedUniform(0.,1.);
        for(int step=0;step<2;++step){actual.beginStep();expectedModel.beginStep();
            for(int j=0;j<3;++j){impact.pos={j==1?2.:.625,.125};impact.particle.inCellL=-1;
                const double energyJ=expectedModel.sampleEnergyJ(expectedRng,expectedUniform);const auto direction=expectedModel.sampleDirection({1.,0.,0.},expectedRng,expectedUniform);
                const auto created=actual.spawnEmbedded(mesh,impact,{1.,0.},kinematics,actualRng,actualUniform);
                check(created.has_value()==(j!=1) && actualRng==expectedRng,"host rejection follows full energy/direction sampling");
                if(created){near(created->energyJ/PhysicalConstants::eV_J,energyJ/PhysicalConstants::eV_J,"placement preserves Gaussian state");near(created->particle.position.x,.655,"embedded injection distance");
                    const double speed=std::sqrt(2.*energyJ/kinematics.massKg)/kinematics.velocityScale;
                    near(created->particle.velocity.x,speed*direction[0],"embedded velocity x");near(created->particle.velocity.y,speed*direction[1],"embedded velocity y");near(created->particle.vz,speed*direction[2],"embedded velocity z");}
            }
        }
    }
}
void emissionStage(QuadtreeMesh& mesh, ImmersedGeometry<QuadtreeMesh>& geometry,
    DictionaryGeometry& dictionary, PicEmbeddedGeometry& embedded)
{
    using Model=PicSecondaryEmission;using Mode=Model::Mode;using Side=ExplicitPIC::Pusher::DomainSide;
    static_assert(std::is_same_v<decltype(std::declval<const Model&>().stats()),const Model::Stats&>);
    const Model::BoundaryGeometry context{mesh,&geometry,dictionary,embedded};
    const Model::DomainBox box{{0.,0.},{1.,1.},{1e-6,1e-6}};
    SEEConfig config;config.enabled=true;config.gamma=2.;config.emittedSpecies="electron";
    config.cathodeSide="LEFT";config.emittedEnergyModel="mono";config.emittedEnergy_eV=5.;
    config.angular="fixed";config.injectEps_m=.06;config.maxEmitMacro=3;
    Particle ion;ion.position={-.1,.125};ion.velocity={3.,4.};ion.vz=12.;ion.w=0.;
    ion.inCellL=mesh.findCellLabel({.125,.125});
    ExplicitPIC::Pusher::DomainImpactBuffer domain;domain.Push(ion,Side::kLeft);domain.Push(ion,Side::kRight);domain.Push(ion,Side::kUnknown);
    ExplicitPIC::Pusher::EmbeddedImpactBuffer walls;walls.Push(ion,{.625,.125},mesh.findCellLabel({.625,.125}),false);
    const auto& impact=walls.impacts().front();
    const std::vector<ParticleType> species{{"ion",2.*PhysicalConstants::eV_J,1.},{"electron",2.*PhysicalConstants::eV_J,-1.}};
    const std::unordered_map<std::string,std::size_t> index{{"ion",0},{"electron",1}};
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        Model model(config,mode);std::unordered_set<int> ids{1};std::vector<double> weights{6.,2.};
        model.initialize(ids,species,weights,index,2.,2.);ids.clear();weights[1]=9.;
        check(model.enabled() && model.emittedIndex()==1,"SEE resolves its target after restore and retains prepared configuration");
        const auto& stats=model.stats();model.beginStep(7);
        std::vector<Particle> pending;ParticleSpawnBuffer named;
        Model::Engine rng(912),reference(912);Model::Uniform uniform(0.,1.),refUniform(0.,1.);
        if(mode==Mode::Explicit) {
            model.processExplicit(context,species,1,6.,domain,&walls,box,{1},pending,rng,uniform);
            check(pending.empty() && rng==reference,"negative incident charge bypasses explicit emission and RNG");
            model.processExplicit(context,species,0,6.,domain,&walls,box,{1},pending,rng,uniform);
            (void)refUniform(reference);(void)refUniform(reference);
            check(pending.size()==6 && stats.ionImpactMacroCathode==2 && stats.emittedMacro==6 && rng==reference,"domain then embedded emission retain capped counts and exact RNG");
            near(stats.ionImpactPhysCathode,12.,"explicit impact weights");near(stats.ionImpactEnergyJCathode/PhysicalConstants::eV_J,8112.,"three-dimensional incident SI energy");
            near(stats.emittedPhys,24.,"cap changes macro weight rather than emitted physical population");
            near(stats.emittedEnergyJ/PhysicalConstants::eV_J,120.,"pre-thinning emitted energy");
            model.accountRetained(pending,species,weights);
            near(stats.injectedCharge[0],-12.,"domain injection charge excludes embedded emission");near(stats.injectedCharge[1],0.,"unselected domain side stays empty");
            for(std::size_t j=0;j<pending.size();++j) {
                near(pending[j].w,4.,"emission uses initialized target weight despite caller changes");
                near(pending[j].position.x,j<3?.03:.655,"mode scale and domain/embedded placement order");
                near(pending[j].velocity.x,std::sqrt(5.)/2.,"emitted SI speed conversion");
            }
            PicSpawnManagement manager(PicSpawnManagement::Mode::Explicit,{},index,weights);
            long long retained=0;double leftCharge=0.;
            for(std::size_t i=0;i<pending.size();++i)if(refUniform(reference)<=.5){++retained;if(i<3)leftCharge-=8.;}
            const auto thinned=model.thinExplicit(manager,pending,true,.5,rng,uniform);
            model.accountRetained(pending,species,weights);
            check(retained>0 && retained<6 && stats.emittedMacro==retained && thinned.kept==retained && rng==reference,"retained report follows actual thinning without changing RNG");
            near(stats.emittedPhys,retained*8.,"compensated retained weight");
            near(stats.emittedEnergyJ/PhysicalConstants::eV_J,retained*40.,"retained velocity-based energy");
            near(stats.injectedCharge[0],leftCharge,"retained domain charge follows compacted origin tags");
            for(const auto& particle:pending)(void)refUniform(reference);
            model.thinExplicit(manager,pending,true,0.,rng,uniform);model.accountRetained(pending,species,weights);
            check(pending.empty() && rng==reference && stats.emittedMacro==0 && stats.emittedPhys==0. && stats.emittedEnergyJ==0. && stats.emittedEnergyJ2==0.,"all-discarded emission contributes zero count, weight and energy moments");
            check(stats.injectedCharge==std::array<double,4>{} && stats.ionImpactMacroCathode==2 && stats.ionImpactPhysCathode==12.,"empty emission clears charge while preserving incident-ion totals");

        } else {
            model.processImplicitEmbedded(context,impact,0,species,0,6.,676.*PhysicalConstants::eV_J,named,rng,uniform);
            model.processImplicitEmbedded(context,impact,1,species,1,6.,676.*PhysicalConstants::eV_J,named,rng,uniform);
            check(named.empty() && rng==reference && stats.ionImpactMacroCathode==0,"implicit boundary and charge gates precede impact accounting");
            model.processImplicitEmbedded(context,impact,1,species,0,6.,676.*PhysicalConstants::eV_J,named,rng,uniform);
            (void)refUniform(reference);
            check(named.size()==3 && rng==reference && stats.ionImpactMacroCathode==1,"implicit embedded stage preserves count draw");
            near(stats.ionImpactPhysCathode,6.,"implicit uses caller-sanitized weight");near(stats.ionImpactEnergyJCathode/PhysicalConstants::eV_J,4056.,"implicit uses caller-computed impact energy");
            for(const auto& item:named.items())check(item.species_name=="electron" && item.particle.w==4.,"implicit pending buffer keeps live target name and compensated weights");
        }
        bool rejected=false;try{model.initialize({},species,weights,index,2.,2.);}catch(const std::logic_error&){rejected=true;}
        check(rejected && model.enabled() && stats.ionImpactMacroCathode>0,"duplicate initialization preserves existing configuration and statistics");
        model.beginStep(8);check(&model.stats()==&stats && stats.ionImpactMacroCathode==0 && stats.emittedMacro==0,"next step clears reports without replacing the borrowed report object");
        auto invalidDictionary=dictionary;invalidDictionary.boundaryById[1].levelSetIndex=99;
        Model::BoundaryGeometry missing{mesh,&geometry,invalidDictionary,embedded};const auto before=rng;
        if(mode==Mode::Explicit)model.processExplicit(missing,species,0,6.,{},&walls,box,{1},pending,rng,uniform);
        else model.processImplicitEmbedded(missing,impact,1,species,0,6.,676.*PhysicalConstants::eV_J,named,rng,uniform);
        check(rng==before && stats.ionImpactMacroCathode==1,"failed normal resolution follows impact accounting and precedes count sampling");
    }
    {
        auto settings=config;settings.cathodeSides={"left","right","bottom","top"};
        const std::vector<double> weights{6.,2.};Model model(settings,Mode::Explicit);
        model.initialize({1},species,weights,index,2.,2.);
        PicSpawnManagement manager(PicSpawnManagement::Mode::Explicit,{},index,weights);
        Model::Engine rng(833),reference(833);Model::Uniform uniform(0.,1.),refUniform(0.,1.);
        std::vector<Particle> pending;
        for(double probability:{1.,.5,0.}) {
            model.beginStep();pending.clear();ExplicitPIC::Pusher::DomainImpactBuffer allSides;
            for(auto side:{Side::kLeft,Side::kRight,Side::kBottom,Side::kTop})allSides.Push(ion,side);
            model.processExplicit(context,species,0,6.,allSides,&walls,box,{1},pending,rng,uniform);
            check(pending.size()==15,"four domain sides plus embedded emission");
            for(int i=0;i<5;++i)(void)refUniform(reference);
            std::array<double,4> expected{};long long retained=0;
            for(int i=0;i<15;++i)if(probability==1. || refUniform(reference)<=probability){++retained;if(i<12)expected[i/3]-=4./probability;}
            model.thinExplicit(manager,pending,true,probability,rng,uniform);model.accountRetained(pending,species,weights);
            check(rng==reference && model.stats().emittedMacro==retained,"tag tracking preserves all/no/partial-retention RNG");
            for(std::size_t side=0;side<4;++side)near(model.stats().injectedCharge[side],expected[side],"each side receives only its retained compensated charge");
            const double physical=retained?retained*4./probability:0.;
            near(model.stats().emittedPhys,physical,"domain plus embedded retained population");
            near(model.stats().emittedEnergyJ2/(PhysicalConstants::eV_J*PhysicalConstants::eV_J),25.*physical,"retained energy second moment");
            const auto charge=model.stats().injectedCharge;model.accountRetained(pending,species,weights);
            check(model.stats().injectedCharge==charge,"repeated accounting replaces rather than accumulates charge");
        }
        pending.push_back(ion);const auto saved=rng;bool rejected=false;
        try{model.thinExplicit(manager,pending,true,.5,rng,uniform);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected && rng==saved && pending.size()==1,"out-of-band pending mutation rejected before thinning");
        rejected=false;try{model.accountRetained(pending,species,weights);}catch(const std::logic_error&){rejected=true;}
        check(rejected,"accounting rejects missing boundary provenance");
    }
    for(Mode mode:{Mode::Explicit,Mode::Implicit})for(int gate=0;gate<4;++gate) {
        auto settings=config;auto targets=species;auto weights=std::vector<double>{6.,2.};auto names=index;
        if(gate==0)settings.enabled=false;
        if(gate==1)names.clear();
        if(gate==2)targets[1]=ParticleType("electron",0.,-1.);
        if(gate==3)weights[1]=0.;
        BoundaryLogs logs;Model model(settings,mode);model.initialize({1},targets,weights,names,2.,2.);model.beginStep();
        Model::Engine rng(77),reference(77);Model::Uniform uniform(0.,1.),refUniform(0.,1.);
        std::vector<Particle> pending;ParticleSpawnBuffer named;
        if(mode==Mode::Explicit)model.processExplicit(context,targets,0,6.,domain,nullptr,box,{},pending,rng,uniform);
        else model.processImplicitEmbedded(context,impact,1,targets,0,6.,1.,named,rng,uniform);
        if(mode==Mode::Explicit && gate==2)(void)refUniform(reference);
        check(pending.empty() && named.empty() && rng==reference,"eligibility retains mode-specific mass gate and RNG consumption");
        check(model.stats().ionImpactMacroCathode==(mode==Mode::Explicit?1:0),"disabled/unusable emission retains each mode's impact tally gate");
    }
    {
        auto settings=config;settings.emittedEnergyModel="maxwell";settings.emittedTe_eV=2.;settings.angular="cosine";
        Model model(settings,Mode::Explicit),kernel(settings,Mode::Explicit);
        model.initialize({1},species,{6.,2.},index,2.,2.);model.beginStep();kernel.beginStep();
        Model::Engine rng(327),reference(327);Model::Uniform uniform(0.,1.),refUniform(0.,1.);
        const int count=kernel.sampleCount(6.,reference,refUniform);
        for(int i=0;i<count;++i){(void)kernel.sampleEnergyJ(reference,refUniform);(void)kernel.sampleDirection(Side::kLeft,reference,refUniform);}
        std::vector<Particle> pending;
        model.processExplicit(context,species,0,6.,domain,nullptr,{{2.,2.},{3.,3.},{1e-6,1e-6}},{},pending,rng,uniform);
        check(pending.empty() && rng==reference && model.stats().ionImpactMacroCathode==1 && model.stats().emittedMacro==0,"failed placement keeps all count/energy/direction draws without crediting emission");
    }
    {
        BoundaryLogs logs;Model model(config,Mode::Implicit);model.initialize({1},species,{6.,2.},index,2.,2.,true);
        Model::Engine rng(53);Model::Uniform uniform(0.,1.);ParticleSpawnBuffer pending;
        for(int step=1;step<=10;++step){model.beginStep(step);model.processImplicitEmbedded(context,impact,1,species,0,6.,1.,pending,rng,uniform);}
        const auto log=logs.output.str();std::size_t pos=0,count=0;
        while((pos=log.find("[ImplicitPIC][SEE][",pos))!=std::string::npos){++count;++pos;}
        check(count==25 && log.find("step=8 ")!=std::string::npos && log.find("step=9 ")==std::string::npos,"implicit debug budget spans steps and includes original impact/yield/emit stages");
    }
}
void resamplingStage()
{
    using Step=ParticleManager::Step;using Mode=Step::Mode;
    static_assert(!std::is_copy_constructible_v<Step> && !std::is_move_constructible_v<Step>);
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph faces;faces.build(mesh);
    PicEmbeddedGeometry empty(PicEmbeddedGeometry::Mode::Explicit);
    ParticleManagementConfig config;config.enabled=true;config.perCellMax=2;config.splitEnabled=false;
    ParticleType seed("test",1.,1.);
    for(int j=0;j<5;++j){Particle p;p.position={.125,.125};p.velocity={2.,3.};p.vz=4.;p.w=2.;p.inCellL=mesh.findCellLabel(p.position);p.inBackgroundCellL=mesh.findBackgroundCellLabel(p.position);seed.particle.push_back(p);}
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        Step step(config,mesh,faces,empty,mode);check(!step.skipMergeCells(),"no cut metrics means no protection mask");
        ParticleType particles=seed,referenceParticles=seed;std::mt19937_64 rng(135),reference(135);ParticleManagementStats legacy;
        ParticleManager::Apply(config,mesh,referenceParticles,0,referenceParticles.mass(),1.,reference,legacy,false,nullptr);
        const auto report=step.apply(particles,0,1.,1.,rng,true);
        check(report.macroBefore==5 && report.macroAfter==2 && report.counts.mergedPairs==3,"resampling stage performs real pair merges");
        near(report.physBefore,10.,"pre-merge physical population");near(report.physAfter,10.,"merge conserves physical population");
        double px=0.,py=0.,pz=0.;for(const auto& p:particles.particle){px+=p.w*p.velocity.x;py+=p.w*p.velocity.y;pz+=p.w*p.vz;}
        near(px,20.,"merge momentum x");near(py,30.,"merge momentum y");near(pz,40.,"merge momentum z");check(rng==reference,"merge stage preserves original kernel RNG consumption");
        check(particles.particle.size()==referenceParticles.particle.size(),"merge stage preserves original kernel count");
        for(std::size_t j=0;j<particles.particle.size();++j){near(particles.particle[j].position.x,referenceParticles.particle[j].position.x,"original merge position x");near(particles.particle[j].position.y,referenceParticles.particle[j].position.y,"original merge position y");near(particles.particle[j].w,referenceParticles.particle[j].w,"original merge weight");}
        check(report.weights.max>=report.weights.p50 && report.weights.top1PhysFrac>0.,"merged weight-tail report");
    }
    config.mergeEnabled=false;config.splitEnabled=true;config.splitWMax=3.;config.splitMaxPerCell=8;
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        ParticleType particles("split",1.,1.);particles.particle.assign(2,seed.particle[0]);for(auto& p:particles.particle)p.w=8.;
        ParticleType reference=particles;std::mt19937_64 rng(921),refRng(921);ParticleManagementStats legacy;
        ParticleManager::Apply(config,mesh,reference,0,reference.mass(),1.,refRng,legacy,false,nullptr);
        Step step(config,mesh,faces,empty,mode);const auto report=step.apply(particles,0,1.,1.,rng,false);
        check(report.counts.splitCount>0 && report.counts.splitEmergencyCount>0 && report.macroAfter>report.macroBefore,"weight-triggered split executed");
        near(report.physAfter,16.,"splitting conserves physical population");check(rng==refRng && particles.particle.size()==reference.particle.size(),"stage preserves direct-kernel RNG and count");
        for(std::size_t j=0;j<particles.particle.size();++j){near(particles.particle[j].position.x,reference.particle[j].position.x,"original split jitter x");near(particles.particle[j].position.y,reference.particle[j].position.y,"original split jitter y");near(particles.particle[j].w,reference.particle[j].w,"original split weight");}
        check(mode==Mode::Explicit?report.weights.max>0.:report.weights.max==0.,"weight diagnostics keep mode-specific gate");
    }
    config.splitEnabled=false;ParticleType tails("tail",1.,1.);
    for(double w:{0.,1.,3.,5.,9.}){Particle p=seed.particle[0];p.w=w;tails.particle.push_back(p);}
    std::mt19937_64 rng(14);Step tailStep(config,mesh,faces,empty,Mode::Explicit);
    const auto tail=tailStep.apply(tails,0,2.,1.,rng,false);near(tail.physAfter,20.,"report-only default-weight fallback");
    near(tail.weights.p50,3.,"lower-index median");near(tail.weights.p90,5.,"lower-index upper quantile");near(tail.weights.p99,5.,"upper quantile is not interpolated");near(tail.weights.max,9.,"maximum macro weight");near(tail.weights.top1PhysFrac,.45,"top one percent includes at least one macro");
    DictionaryGeometry dict;dict.levelSetName={"wall"};dict.innerPoints={{.1,.5}};dict.moving={false};
    dict.linesOnInterface={{{0,false}}};dict.phaseName={"solid"};dict.levelSetOperation={{LevelSetOperation::non}};
    dict.levelSetIndexA={{{0,false}}};dict.levelSetIndexB={{}};dict.boundaryLines={{{0,false}}};
    Line line;line.points[0]={.625,0.};line.points[1]={.625,1.};dict.lines={line};
    DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.kind=GeometricBoundary::fixed;boundary.levelSetIndex=0;boundary.phase=0;boundary.lineIds={0};dict.boundaries={boundary};dict.boundaryById[1]=boundary;
    ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh);geometry.levelSets[0].initialize(mesh,0.);geometry.LS_dx(0).initialize(mesh,1.);geometry.LS_dy(0).initialize(mesh,0.);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry.levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.625;
    PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);embedded.resolveSelection(dict,1);embedded.rebuildComposite(mesh,&geometry);embedded.rebuildMetrics(mesh,&geometry,faces,{1});
    check(embedded.metricsReady(),"resampling test has real cut metrics");
    particleBoundaryOperations(mesh,geometry,dict,embedded);
    emissionStage(mesh,geometry,dict,embedded);
    config.mergeEnabled=true;Step explicitCut(config,mesh,faces,embedded,Mode::Explicit),implicitCut(config,mesh,faces,embedded,Mode::Implicit,true);
    check(explicitCut.skipMergeCells() && implicitCut.skipMergeCells(),"cut masks belong to each stage");
    const int farLeft=mesh.cells(mesh.findCellLabel({.125,.125})).index();
    const int wallCell=mesh.cells(mesh.findCellLabel({.625,.125})).index();
    check((*explicitCut.skipMergeCells())[farLeft]==0 && (*implicitCut.skipMergeCells())[farLeft]!=0,"explicit cut-neighbor band differs from implicit solid/face propagation");
    check((*explicitCut.skipMergeCells())[wallCell]!=0 && (*implicitCut.skipMergeCells())[wallCell]!=0,"wall cell is protected in both modes");
    ParticleType protectedParticles=seed;const auto protectedReport=implicitCut.apply(protectedParticles,0,1.,1.,rng,true);
    check(protectedReport.macroAfter==5 && protectedReport.counts.mergedPairs==0,"mask suppresses population changes in protected cell");
    Step disabledCut(config,mesh,faces,embedded,Mode::Implicit,false);check(!disabledCut.skipMergeCells(),"implicit cut-metric gate retained");
}
void populationStage()
{
    using Step=ParticleManager::Step;using Mode=Step::Mode;
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph faces;faces.build(mesh);
    PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);
    std::vector<ParticleType> seed;
    for(const char* name:{"merge","split","empty","fallback"})seed.emplace_back(name,1e-30,1.);
    for(std::size_t si=0;si<seed.size();++si) {
        const int count=si==0?5:si==2?0:2;
        for(int j=0;j<count;++j) {
            Particle p;p.position={.125,.125};p.velocity={double(j+1),2.};p.vz=3.;
            p.w=si==1?8.:si==3||j==0?0.:2.;p.rngState=j%2?100+j:0;
            p.trackId=1000*si+j;p.emitChance=.2;p.levelSetValues={.1,.3};p.weight[1]=.7;
            p.inCellL=mesh.findCellLabel(p.position);p.inBackgroundCellL=mesh.findBackgroundCellLabel(p.position);
            seed[si].particle.push_back(p);
        }
        seed[si].number=count;
    }
    const auto sameSpecies=[](const ParticleType& a,const ParticleType& b) {
        check(a.name()==b.name() && a.number==b.number && a.particle.size()==b.particle.size(),"population species identity/count");
        for(std::size_t j=0;j<a.particle.size();++j) {
            const auto& p=a.particle[j];const auto& q=b.particle[j];
            check(p.position.x==q.position.x && p.position.y==q.position.y && p.velocity.x==q.velocity.x && p.velocity.y==q.velocity.y && p.vz==q.vz && p.w==q.w,"population position/velocity/weight exactly match single-species reference");
            check(p.inCellL==q.inCellL && p.inBackgroundCellL==q.inBackgroundCellL && p.discreteVertL==q.discreteVertL && p.rngState==q.rngState && p.trackId==q.trackId && p.emitChance==q.emitChance && p.levelSetValues==q.levelSetValues,"population particle caches and identities exactly match reference");
            for(int k=0;k<4;++k)check(p.weight[k]==q.weight[k],"population interpolation weights");
        }
    };
    const auto sameCounts=[](const ParticleManagementStats& a,const ParticleManagementStats& b) {
        check(a.mergedPairs==b.mergedPairs && a.splitCount==b.splitCount && a.splitEmergencyCount==b.splitEmergencyCount && a.splitCountBasedCount==b.splitCountBasedCount,"population count aggregation");
    };
    for(Mode mode:{Mode::Explicit,Mode::Implicit})for(bool reportWeights:{false,true})
    for(int operation=0;operation<4;++operation)for(int arraySize:{0,1,4}) {
        ParticleManagementConfig config;config.enabled=true;config.perCellMax=2;config.perCellMin=6;
        config.mergeEnabled=operation==0;config.splitEnabled=operation==1||operation==2;
        config.splitWMax=operation==1?3.:0.;config.splitByCountEnabled=operation==2;
        config.mergeMaxDeltaE_bySpecies_eV={1.,0.};
        Step stage(config,mesh,faces,embedded,mode);
        const std::vector<double> weights{-3.,8.}; // Exercise signed physical changes and missing weights.
        auto actual=seed,reference=seed;std::mt19937_64 rng(1481),refRng(1481);
        std::vector<ParticleManagementReport> reports;ParticleManagementStats expected;
        for(std::size_t si=0;si<reference.size();++si) {
            reports.push_back(stage.apply(reference[si],si,si<weights.size()?weights[si]:1.,1.,refRng,reportWeights));
            const auto& counts=reports.back().counts;
            expected.mergedPairs+=counts.mergedPairs;expected.splitCount+=counts.splitCount;
            expected.splitEmergencyCount+=counts.splitEmergencyCount;expected.splitCountBasedCount+=counts.splitCountBasedCount;
        }
        check(operation!=0||expected.mergedPairs>0,"population test performs merges");
        check(operation!=1||expected.splitEmergencyCount>0,"population test performs weight splits");
        check(operation!=2||expected.splitCountBasedCount>0,"population test performs count splits");
        std::vector<long long> macro(arraySize,17);std::vector<double> physical(4-arraySize,23.);
        auto expectedMacro=macro;auto expectedPhysical=physical;std::size_t observed=0;
        const auto totals=stage.applyAll(actual,weights,1.,rng,reportWeights,
            [&](const ParticleType& sp,const ParticleManagementReport& report) {
                const std::size_t si=observed++;const auto& r=reports[si];
                check(&sp==&actual[si],"observer sees live species in storage order, including empty species");
                sameSpecies(sp,reference[si]);
                for(std::size_t later=si+1;later<actual.size();++later)sameSpecies(actual[later],seed[later]);
                check(report.macroBefore==r.macroBefore && report.macroAfter==r.macroAfter && report.physBefore==r.physBefore && report.physAfter==r.physAfter,"population before/after report");
                sameCounts(report.counts,r.counts);
                check(report.weights.p50==r.weights.p50 && report.weights.p90==r.weights.p90 && report.weights.p99==r.weights.p99 && report.weights.max==r.weights.max && report.weights.top1PhysFrac==r.weights.top1PhysFrac,"mode-specific weight-tail gate");
                if(mode==Mode::Implicit) {
                    if(si<expectedMacro.size())expectedMacro[si]+=r.macroBefore-r.macroAfter;
                    if(si<expectedPhysical.size())expectedPhysical[si]+=r.physBefore-r.physAfter;
                }
                check(macro==expectedMacro && physical==expectedPhysical,"signed implicit losses are accumulated before observer; explicit losses stay untouched");
            },&macro,&physical);
        check(observed==seed.size() && rng==refRng,"population observer count and global RNG sequence");sameCounts(totals,expected);
        if(mode==Mode::Implicit && arraySize==4 && (operation==1||operation==2))check(macro[operation==1?1:0]<17,"splits retain negative macro loss");
        if(mode==Mode::Implicit && arraySize==0 && operation==0)check(physical[0]<23.,"negative physical loss is not clamped");
        auto silent=seed;std::mt19937_64 silentRng(1481);
        sameCounts(stage.applyAll(silent,weights,1.,silentRng,reportWeights),expected);
        check(silentRng==refRng,"absent observer and loss arrays preserve RNG");
        for(std::size_t si=0;si<silent.size();++si)sameSpecies(silent[si],reference[si]);
        std::vector<ParticleType> noSpecies;auto emptyRng=rng;
        sameCounts(stage.applyAll(noSpecies,{},1.,rng,true,[&](const auto&,const auto&){check(false,"empty population must not call observer");}),{});
        check(rng==emptyRng,"empty population consumes no RNG");

        auto interrupted=seed;std::mt19937_64 interruptedRng(1481),firstRng(1481);
        auto first=seed[0];stage.apply(first,0,weights[0],1.,firstRng,reportWeights);
        std::vector<long long> partialLoss(4,17);bool thrown=false;
        try { stage.applyAll(interrupted,weights,1.,interruptedRng,reportWeights,
            [](const auto&,const auto&){throw std::runtime_error("observer failure");},&partialLoss); }
        catch(const std::runtime_error& e){thrown=std::string(e.what())=="observer failure";}
        check(thrown && interruptedRng==firstRng,"observer failure propagates after current mutation without consuming later RNG");
        sameSpecies(interrupted[0],first);
        check(partialLoss[0]==17+(mode==Mode::Implicit?reports[0].macroBefore-reports[0].macroAfter:0),"observer failure retains completed accounting");
        for(std::size_t si=1;si<seed.size();++si){sameSpecies(interrupted[si],seed[si]);check(partialLoss[si]==17,"later accounting untouched on failure");}
    }
}
void spawnManagement()
{
    using Manager=PicSpawnManagement;using Mode=Manager::Mode;
    const std::unordered_map<std::string,std::size_t> names{{"electron",0},{"ion",1}};
    const std::vector<double> weights{2.,5.};ParticleManagementConfig config;
    config.softCapTotal=10;config.hardCapTotal=15;config.minKeepProb=.3;
    std::vector<ParticleType> species;species.emplace_back("electron",1.,-1.);species.emplace_back("ion",2.,1.);
    species[0].particle.resize(12);
    Manager explicitManager(Mode::Explicit,config,names,weights),implicitManager(Mode::Implicit,config,names,weights);
    near(explicitManager.keepProbability(species,8,true),.5,"spawn cap includes pending and live population");
    near(explicitManager.keepProbability(species,100,true),.3,"minimum thinning probability");
    near(explicitManager.keepProbability(species,100,false),1.,"inactive thinning ignores caps");
    config.softCapTotal=1;near(explicitManager.keepProbability(species,8,true),.5,"population settings are a snapshot");
    config.softCapTotal=0;config.hardCapTotal=5;config.minKeepProb=1.25;
    Manager unusual(Mode::Implicit,config,names,weights);near(unusual.keepProbability(species,0,true),1.25,"legacy minimum above one is not clamped again");
    for(const auto* manager:{&explicitManager,&implicitManager}) {
        ParticleSpawnBuffer pending;std::vector<Particle> expected;double before=0.,after=0.;
        std::mt19937_64 rng(2347),reference(2347);std::uniform_real_distribution<double> uniform(0.,1.),refUniform(0.,1.);
        for(int j=0;j<21;++j) {
            Particle p;p.position={double(j),.5};p.w=j%2?0.:3.;
            const std::string name=j%3==0?"electron":j%3==1?"ion":"unknown";
            const double weight=p.w>0.?p.w:j%3==0?2.:j%3==1?5.:1.;before+=weight;
            pending.Push(name,p);
            if(refUniform(reference)<=.4){p.w=weight/.4;expected.push_back(p);after+=p.w;}
        }
        const auto stats=manager->thin(pending,true,.4,rng,uniform);
        check(stats.in==21 && stats.kept==static_cast<long long>(expected.size()) && stats.drop==21-stats.kept && stats.drop>0,"stochastic thinning count");
        near(stats.physIn,before,"source-name default weights including unknown fallback");near(stats.physKept,after,"inverse-probability weight compensation");
        near(stats.reportKeepProb(),.4,"per-call reported probability");check(rng==reference,"one uniform draw per candidate");
        for(std::size_t j=0;j<expected.size();++j) {
            near(pending.items()[j].particle.position.x,expected[j].position.x,"thinning keeps original order");
            near(pending.items()[j].particle.w,expected[j].w,"retained weight");
        }
        const auto untouched=manager->thin(pending,false,.1,rng,uniform);
        check(rng==reference && untouched.drop==0,"disabled group consumes no random numbers");near(untouched.reportKeepProb(),1.,"disabled group reports one");
        Manager::Stats combined;combined.accumulate(stats);combined.accumulate(untouched);
        near(combined.reportKeepProb(),(.4*21.+double(untouched.in))/double(21+untouched.in),"implicit half-step probability weighted by input counts");
        std::vector<Particle> see(13);for(std::size_t j=0;j<see.size();++j){see[j].position={double(j),.5};see[j].w=0.;}
        std::size_t keep=0;for(std::size_t j=0;j<see.size();++j)if(refUniform(reference)<=.5)++keep;
        const auto seeStats=manager->thin(see,1,true,.5,rng,uniform);
        check(see.size()==keep && rng==reference,"single-species SEE shares thinning RNG kernel");near(seeStats.physIn,65.,"SEE fallback weight");
        for(const auto& p:see)near(p.w,10.,"SEE weight compensation");
        std::vector<Particle> tagged(20);std::vector<int> tags,expectedTags;
        for(int i=0;i<20;++i){tagged[i].position.x=i;tagged[i].w=2.;tags.push_back(100+i);if(refUniform(reference)<=.4)expectedTags.push_back(100+i);}
        manager->thin(tagged,0,true,.4,rng,uniform,&tags);
        check(tags==expectedTags && tags.size()==tagged.size() && rng==reference,"optional tags use the same stable selection and RNG");
        for(std::size_t i=0;i<tags.size();++i){near(tagged[i].position.x,tags[i]-100,"tag follows surviving particle");near(tagged[i].w,5.,"tagged thinning preserves weight compensation");}
        manager->thin(tagged,0,false,.1,rng,uniform,&tags);manager->thin(tagged,0,true,1.,rng,uniform,&tags);
        check(tags==expectedTags && rng==reference,"disabled and unit-probability thinning preserve tags without draws");
        auto badTags=tags;badTags.push_back(-1);bool badSize=false;
        try{manager->thin(tagged,0,true,.4,rng,uniform,&badTags);}catch(const std::invalid_argument&){badSize=true;}
        check(badSize && rng==reference && tagged.size()==tags.size(),"tag mismatch rejected before particle mutation or RNG");
        for(std::size_t i=0;i<tagged.size();++i)(void)refUniform(reference);
        manager->thin(tagged,0,true,0.,rng,uniform,&tags);manager->thin(tagged,0,true,.5,rng,uniform,&tags);
        check(tagged.empty() && tags.empty() && rng==reference,"all discarded tags clear and empty batch draws nothing");
        ParticleSpawnBuffer empty;const auto noStats=manager->thin(empty,true,.2,rng,uniform);
        check(noStats.in==0 && noStats.reportKeepProb()==1. && rng==reference,"empty batch draws nothing");
        Particle infinite;infinite.w=std::numeric_limits<double>::infinity();empty.Push("electron",infinite);
        const auto nonfinite=manager->thin(empty,true,0.,rng,uniform);
        check(nonfinite.kept==0,"zero probability drops candidate");
        check(manager==&explicitManager?nonfinite.physKept==0.:std::isnan(nonfinite.physKept),"retain original nonfinite tally arithmetic");
    }
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    ParticleSpawnBuffer pending;Particle valid;valid.position={.2,.2};valid.w=0.;valid.inCellL=-1;valid.inBackgroundCellL=-1;
    Particle outside=valid;outside.position={2.,.2};Particle weighted=valid;weighted.w=7.;
    pending.Push("electron",valid);pending.Push("unknown",valid);pending.Push("ion",outside);pending.Push("electron",weighted);
    species[0].particle.clear();species[0].number=0;species[1].number=99;std::vector<unsigned char> touched(2,0);
    explicitManager.commit(pending,mesh,species,&touched);
    check(species[0].particle.size()==2 && species[1].particle.size()==1 && touched==std::vector<unsigned char>({1,1}),"explicit appends known species and marks touched");
    check(species[0].particle[0].inCellL==-1 && species[0].particle[0].w==0.,"explicit leaves host/weight preparation to original caller stage");
    check(species[0].number==2 && species[1].number==1 && pending.size()==4,"commit retains diagnostic buffer");
    for(auto& sp:species){sp.particle.clear();sp.number=99;}
    {BoundaryLogs logs;implicitManager.commit(pending,mesh,species,nullptr,true);
     check(logs.errors.str().find("spawned species 'unknown'")!=std::string::npos,"implicit MCC-only missing-species warning");}
    check(species[0].particle.size()==2 && species[1].particle.empty() && species[1].number==0,"implicit repairs metadata for every species and drops outside host");
    check(species[0].particle[0].inCellL>=0 && species[0].particle[0].inBackgroundCellL>=0,"implicit repairs incoming host labels");
    near(species[0].particle[0].w,2.,"implicit fills missing weight");near(species[0].particle[1].w,7.,"implicit preserves set weight");
    check(pending.items()[0].particle.inCellL==-1,"implicit commits a copy and keeps buffer intact");
    ParticleSpawnBuffer unknown;unknown.Push("unknown",valid);
    {BoundaryLogs logs;implicitManager.commit(unknown,mesh,species);check(logs.errors.str().empty(),"SEE/source unknown targets are silent");}
    species[0].particle[0].inCellL=-1;bool rejected=false;
    {BoundaryLogs logs;try{implicitManager.commit(unknown,mesh,species,nullptr,true);}catch(const std::runtime_error&){rejected=true;}}
    check(rejected,"MCC still checks all existing particles after commit");
    std::vector<Particle> see{weighted};touched={0,0};explicitManager.commitExplicitSee(see,1,species,touched);
    check(species[1].number==1 && touched==std::vector<unsigned char>({0,1}) && see.size()==1,"explicit SEE append and retained buffer");
}
void particleSources()
{
    using Sources=PicParticleSources;using Mode=Sources::Mode;
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryGeometry dict;
    Line line;line.points[0]={0.,1.};line.points[1]={0.,0.};dict.lines={line};
    DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.lineIds={0};dict.boundaryById[1]=boundary;
    PicEmbeddedGeometry noGeometry(PicEmbeddedGeometry::Mode::Explicit);
    SpeciesConfig config;config.name="electron";config.q=-2.;config.m=PhysicalConstants::eV_J;config.weight=2.;
    config.source.enabled=true;config.source.boundaryId=1;config.source.numberDensity_m3=8.;
    config.source.every=2;config.drift={2.,0.};
    const std::unordered_map<std::string,std::size_t> names{{"electron",0}};
    std::mt19937_64 rng(97),reference(97);std::uniform_real_distribution<double> uniform(0.,1.),refUniform(0.,1.);
    Sources density(Mode::Explicit,{config},names,&dict,mesh,noGeometry,nullptr);
    ParticleSpawnBuffer pending;
    const std::vector<ParticleType> species{{"electron",PhysicalConstants::eV_J,-2.}};
    const std::vector<double> weights{2.};
    density.emit({0,.5,1.,.5},mesh,noGeometry,nullptr,rng,uniform,pending);
    check(pending.empty() && rng==reference,"source cadence skips all draws");
    config.drift={99.,99.}; // Prepared settings are independent of later input edits.
    density.emit({1,.5,1.,.5},mesh,noGeometry,nullptr,rng,uniform,pending);
    check(pending.size()==2,"density * area * drift * dt / weight");
    const auto densityReport=density.accountRetained(pending,species,weights,names,1.,{1,2,3,4});
    near(densityReport.domainChargeC[0],-8.,"retained source charge");
    near(densityReport.energyJ/PhysicalConstants::eV_J,8.,"retained source kinetic energy");
    for(const auto& item:pending.items()) {
        near(item.particle.velocity.x,2.,"source owns drift snapshot");near(item.particle.vz,0.,"cold source");
        check(item.particle.inCellL>=0 && item.particle.position.x>0.,"valid injected host");
    }
    {
        auto geometryInput=dict;
        for(int id=2;id<=5;++id){auto extra=boundary;extra.id=id;geometryInput.boundaryById[id]=extra;}
        std::vector<SpeciesConfig> emitters;
        for(int id=1;id<=5;++id){auto emitter=config;emitter.drift={2.,0.};emitter.source.every=1;
            emitter.source.boundaryId=id;emitter.source.current_A=16.;emitter.source.maxEmitMacro=4;emitters.push_back(emitter);}
        Sources source(Mode::Explicit,emitters,names,&geometryInput,mesh,noGeometry,nullptr);
        PicSpawnManagement manager(PicSpawnManagement::Mode::Explicit,{},names,weights);
        for(double probability:{1.,.5,0.}) {
            pending.Clear();source.emit({0,1.,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);
            check(pending.size()==20,"same species can have multiple source boundaries");
            reference=rng;std::array<double,4> expected{};long long keep=0;
            for(int i=0;i<20;++i)if(probability==1. || refUniform(reference)<=probability){++keep;if(i<16)expected[i/4]-=4./probability;}
            const auto stats=source.thinExplicit(manager,pending,true,probability,rng,uniform);
            const auto report=source.accountRetained(pending,species,weights,names,2.,{1,2,3,4});
            check(rng==reference && stats.kept==keep && pending.size()==std::size_t(keep),"source thinning retains RNG/order and compensated population");
            for(std::size_t side=0;side<4;++side)near(report.domainChargeC[side],expected[side],"source boundary tags survive compaction");
            near(report.energyJ/PhysicalConstants::eV_J,keep?16.*keep/probability:0.,"all retained sources contribute SI energy including embedded source");
            const auto again=source.accountRetained(pending,species,weights,names,2.,{1,2,3,4});
            check(again.domainChargeC==report.domainChargeC && again.energyJ==report.energyJ,"source accounting is repeatable");
        }
        // Unknown targets retain the original explicit energy/commit gate.
        pending.Clear();source.emit({0,1.,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);
        const auto unknown=source.accountRetained(pending,species,weights,{},1.,{1,2,3,4});
        check(unknown.energyJ==0. && unknown.domainChargeC==std::array<double,4>{},"unknown source targets contribute nothing");
        pending.items_mut().pop_back();const auto saved=rng;bool rejected=false;
        try{source.thinExplicit(manager,pending,true,.5,rng,uniform);}catch(const std::logic_error&){rejected=true;}
        check(rejected && rng==saved && pending.size()==19,"source provenance mismatch rejected before thinning");
        rejected=false;try{source.accountRetained(pending,species,weights,names,1.,{1,2,3,4});}catch(const std::logic_error&){rejected=true;}
        check(rejected,"source accounting rejects missing provenance");
        pending.Clear();source.emit({0,1.,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);
        source.emit({0,1.,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);
        check(pending.size()==40,"explicit source batch supports append without losing provenance");
        near(source.accountRetained(pending,species,weights,names,1.,{1,1,0,0}).domainChargeC[0],-32.,"duplicate domain ids preserve first-side precedence");
    }
    config.drift={-2.,0.};config.Te_eV=1.;config.source.every=1;config.source.current_A=12.;
    Sources thermal(Mode::Implicit,{config},names,&dict,mesh,noGeometry,nullptr);
    rng.seed(97);reference.seed(97);pending.Clear();
    for(int step=0;step<3;++step) {
        const auto oldSize=pending.size();thermal.emit({step,1.,1.,0.},mesh,noGeometry,nullptr,rng,uniform,pending);
        check(pending.size()==oldSize+3,"current source overrides density and appends buffer");
        refUniform(reference);std::normal_distribution<double> normal(0.,1.);
        for(int j=0;j<3;++j) {
            refUniform(reference);const double t=refUniform(reference);
            const double vx=-2.+normal(reference),vy=normal(reference),vz=normal(reference);
            const auto& p=pending.items()[oldSize+j].particle;
            near(p.position.y,std::min(std::max(1.-t,1e-12+1e-6),1.-1e-12-1e-6),"segment sampling");
            near(p.velocity.x,std::abs(vx),"incoming thermal velocity reflection");near(p.velocity.y,vy,"thermal tangential velocity");near(p.vz,vz,"three velocity components");
        }
        check(rng==reference,"Gaussian cache discarded per emitter per call, including odd count");
    }
    auto fluxConfig=config;fluxConfig.source.current_A=0.;fluxConfig.source.numberDensity_m3=4.;
    fluxConfig.Te_eV=8.*3.14159265358979323846;fluxConfig.drift={0.,0.};
    Sources flux(Mode::Implicit,{fluxConfig},names,&dict,mesh,noGeometry,nullptr);pending.Clear();
    flux.emit({0,.5,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);
    check(pending.size()==2,"thermal half-space flux contributes sqrt(kT / (2 pi m))");
    config.source.maxEmitMacro=1;Sources capped(Mode::Explicit,{config},names,&dict,mesh,noGeometry,nullptr);
    pending.Clear();capped.emit({0,1.,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);check(pending.size()==1,"source macro cap");
    // On a reversed domain edge, preserve explicit's line normal and implicit's
    // domain-inward fallback; this distinction affects reflected velocity.
    std::swap(dict.lines[0].points[0],dict.lines[0].points[1]);config.Te_eV=0.;config.drift={2.,0.};
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        Sources source(mode,{config},names,&dict,mesh,noGeometry,nullptr);pending.Clear();
        source.emit({0,1.,1.,1.},mesh,noGeometry,nullptr,rng,uniform,pending);
        check(pending.size()==1,"domain source placed after clamping");near(pending.items()[0].particle.velocity.x,mode==Mode::Explicit?-2.:2.,"mode-specific domain orientation");
    }
    {
        BoundaryLogs logs;auto bad=config;bad.source.boundaryId=0;
        Sources missing(Mode::Explicit,{bad},names,&dict,mesh,noGeometry,nullptr);
        Sources absent(Mode::Implicit,{config},names,nullptr,mesh,noGeometry,nullptr);
        check(missing.empty() && absent.empty() && logs.errors.str().find("missing sourceBoundaryId")!=std::string::npos,"invalid/no geometry sources skipped");
    }
    dict.lines[0].points[0]={.625,0.};dict.lines[0].points[1]={.625,1.};
    dict.levelSetName={"solid"};dict.innerPoints={{.1,.5}};dict.moving={false};dict.linesOnInterface={{{0,false}}};
    dict.phaseName={"solid"};dict.levelSetOperation={{LevelSetOperation::non}};dict.levelSetIndexA={{{0,false}}};
    dict.levelSetIndexB={{}};dict.boundaryLines={{{0,false}}};boundary.kind=GeometricBoundary::fixed;
    boundary.levelSetIndex=0;boundary.phase=0;dict.boundaries={boundary};dict.boundaryById[1]=boundary;
    ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh);geometry.levelSets[0].initialize(mesh,0.);
    geometry.LS_dx(0).initialize(mesh,1.);geometry.LS_dy(0).initialize(mesh,0.);
    for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry.levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.625;
    PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);embedded.resolveSelection(dict,1);embedded.rebuildComposite(mesh,&geometry);
    config.Te_eV=1.;
    Sources explicitSource(Mode::Explicit,{config},names,&dict,mesh,embedded,&geometry);
    Sources implicitSource(Mode::Implicit,{config},names,&dict,mesh,embedded,&geometry);
    for(auto* source:{&explicitSource,&implicitSource}) {
        pending.Clear();source->emit({0,1.,1.,1.},mesh,embedded,&geometry,rng,uniform,pending);
        check(pending.size()==1 && pending.items()[0].particle.position.x>.625 && pending.items()[0].particle.velocity.x>0.,"embedded source selects plasma side");
    }
    for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry.levelSets[0].value[vi]=-1.;
    pending.Clear();reference=rng;refUniform(reference);refUniform(reference);refUniform(reference);
    explicitSource.emit({0,1.,1.,1.},mesh,embedded,&geometry,rng,uniform,pending);
    check(pending.empty() && rng==reference,"explicit current-solid rejection precedes velocity draws");
    const auto rejectedReport=explicitSource.accountRetained(pending,species,weights,names,1.,{1,2,3,4});
    check(rejectedReport.energyJ==0. && rejectedReport.domainChargeC==std::array<double,4>{},"rejected placements contribute no energy or charge");
    implicitSource.emit({0,1.,1.,1.},mesh,embedded,&geometry,rng,uniform,pending);
    check(pending.size()==1,"implicit retains prepared normal and original placement policy");
}
void secondaryEmission()
{
    using Model=PicSecondaryEmission;using Mode=Model::Mode;using Side=ExplicitPIC::Pusher::DomainSide;
    static_assert(!std::is_copy_constructible_v<Model> && !std::is_move_constructible_v<Model>);
    SEEConfig config;config.gamma=2.;Model constant(config,Mode::Explicit);config.gamma=9.;
    near(constant.yield(123.),2.,"emission owns configuration snapshot");
    config.gammaModel="PoWeRlAw";config.gamma0=.5;config.Eref_eV=10.;config.alpha=2.;config.Emin_eV=5.;config.gammaMax=3.;
    Model power(config,Mode::Implicit);near(power.yield(5.),0.,"powerlaw threshold");near(power.yield(15.),.5,"powerlaw scaling");near(power.yield(100.),3.,"yield ceiling");
    config.gammaModel="piecewise";config.piecewiseE_eV={0.,10.,20.};config.piecewiseGamma={0.,1.,3.};
    Model piecewise(config,Mode::Explicit);near(piecewise.yield(-1.),0.,"lower table endpoint");near(piecewise.yield(15.),2.,"table interpolation");near(piecewise.yield(100.),3.,"upper table endpoint");
    config.piecewiseGamma.push_back(7.);Model explicitTail(config,Mode::Explicit),implicitTail(config,Mode::Implicit);
    near(explicitTail.yield(100.),7.,"legacy explicit unmatched table tail");near(implicitTail.yield(100.),3.,"legacy implicit common-prefix table tail");
    config.piecewiseE_eV.clear();config.gamma=.4;Model fallback(config,Mode::Explicit);near(fallback.yield(12.),.4,"short table fallback");
    config=SEEConfig{};config.emittedEnergyModel="maxwell";config.emittedTe_eV=2.;config.angular="cosine";config.maxEmitMacro=0;
    for(Mode mode:{Mode::Explicit,Mode::Implicit}) {
        Model model(config,mode);Model::Engine rng(1789),reference(1789);Model::Uniform uniform(0.,1.),refUniform(0.,1.);std::normal_distribution<double> normal(0.,1.);
        for(int step=0;step<5;++step) {
            model.beginStep();if(mode==Mode::Explicit)normal.reset();
            const int expected=3+(refUniform(reference)<.25?1:0);
            check(model.sampleCount(3.25,rng,uniform)==expected,"stochastic count matches original draw");
            const double a=normal(reference),b=normal(reference),c=normal(reference);
            near(model.sampleEnergyJ(rng,uniform)/PhysicalConstants::eV_J,a*a+b*b+c*c,"Gaussian spare-variate lifetime matches each integrator");
            const double mu=std::sqrt(refUniform(reference));refUniform(reference);
            const auto direction=model.sampleDirection({3.,4.,0.},rng,uniform);
            near(.6*direction[0]+.8*direction[1],mu,"cosine angular sample about inward normal");
            near(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2],1.,"unit emission direction");
            check(rng==reference,"energy and direction keep complete RNG state");
            check(uniform(rng)==refUniform(reference),"source/diagnostic draw interleaving preserved");
        }
    }
    config.emittedEnergyModel="mono";config.emittedEnergy_eV=3.;config.angular="fixed";config.maxEmitMacro=2;
    Model mono(config,Mode::Explicit);Model::Engine rng(91),reference(91);Model::Uniform uniform(0.,1.),refUniform(0.,1.);
    near(mono.sampleEnergyJ(rng,uniform)/PhysicalConstants::eV_J,3.,"mono energy");
    const std::array<Side,5> sides{Side::kLeft,Side::kRight,Side::kBottom,Side::kTop,Side::kUnknown};
    const std::array<std::array<double,3>,5> normals{{{1.,0.,0.},{-1.,0.,0.},{0.,1.,0.},{0.,-1.,0.},{1.,0.,0.}}};
    for(std::size_t j=0;j<sides.size();++j)check(mono.sampleDirection(sides[j],rng,uniform)==normals[j],"domain inward emission directions");
    check(mono.sampleDirection({0.,0.,0.},rng,uniform)==normals[0] && rng==reference,"fixed and mono branches draw no random numbers");
    check(mono.sampleCount(0.,rng,uniform)==0 && rng==reference,"zero expected count draws nothing");
    refUniform(reference);check(mono.sampleCount(8.,rng,uniform)==2 && rng==reference,"integral count still draws before cap");
    config.emittedEnergyModel="exponential";config.emittedEmean_eV=4.;Model exponential(config,Mode::Implicit);
    const double expected=-4.*std::log(std::max(1e-12,refUniform(reference)));
    near(exponential.sampleEnergyJ(rng,uniform)/PhysicalConstants::eV_J,expected,"exponential inverse CDF");check(rng==reference,"exponential consumes one uniform sample");
    config.emittedEnergyModel="maxwell";config.emittedTe_eV=-1.;Model zero(config,Mode::Implicit);
    near(zero.sampleEnergyJ(rng,uniform),0.,"nonpositive thermal energy");check(rng==reference,"nonpositive thermal energy consumes no normal sample");
    // Distribution moments provide an independent physics check beyond the
    // exact engine-sequence contract above (2D normal, 3D velocity space).
    config.emittedTe_eV=2.;config.angular="isotropic";Model isotropic(config,Mode::Implicit);
    double energy=0.,cosine=0.;constexpr int samples=30000;
    for(int j=0;j<samples;++j){energy+=isotropic.sampleEnergyJ(rng,uniform)/PhysicalConstants::eV_J;cosine+=isotropic.sampleDirection({1.,0.,0.},rng,uniform)[0];}
    check(std::abs(energy/samples-3.)<.08 && std::abs(cosine/samples-.5)<.015,"Maxwell energy and isotropic hemisphere moments");
}
struct ScopedEnvironment {
    std::string name,oldValue;bool present;
    static void set(const char* name,const char* value) {
#ifdef _WIN32
        check(_putenv_s(name,value?value:"")==0,"set test environment");
#else
        check((value?setenv(name,value,1):unsetenv(name))==0,"set test environment");
#endif
    }
    ScopedEnvironment(const char* key,const char* value):name(key),present(std::getenv(key)!=nullptr) {
        if(present)oldValue=std::getenv(key);set(key,value);
    }
    ~ScopedEnvironment(){set(name.c_str(),present?oldValue.c_str():nullptr);}
};
void implicitEnvironmentInput()
{
    ScopedEnvironment finite("LAMP_IMPLICIT_CONSOLE_DQREL_WARN",nullptr),
        positive("LAMP_IMPLICIT_AMR_REINIT_ACTIVE_MOVING_MAX_STREAK",nullptr),
        integer("LAMP_IMPLICIT_SEGMENTER_MAXLEVEL",nullptr),
        bounded("LAMP_IMPLICIT_SEGMENTER_MAXSEG",nullptr),
        boolean("LAMP_IMPLICIT_ORBITAVG_FUSED",nullptr),
        debug("LAMP_IMPLICIT_SEE_DEBUG",nullptr),
        path("LAMP_IMPLICIT_INIT_GEOM_SNAPSHOT_PATH",nullptr),
        mode("LAMP_IMPLICIT_AMR_GEOM_REFRESH_MODE",nullptr),
        gate("LAMP_IMPLICIT_AMR_GATE_EFIELD_REL",nullptr),
        omega("LAMP_IMPLICIT_POL_CG_TEMPORAL_OMEGA",nullptr),
        trace("LAMP_SEGMENTER_TRACE",nullptr),hint("LAMP_SEGMENTER_NOHINT",nullptr);
    auto absent=readImplicitPicOverrides();
    check(!absent.consoleDqrelWarn && !absent.amrReinitActiveMovingMaxStreak && !absent.segmenterMaxLevel
        && !absent.segmenterMaxSegments && !absent.orbitAvgFused && !absent.initGeomSnapshotPath
        && !absent.amrGeomRefreshMode && !absent.amrGateEFieldRel && !absent.polarization.temporalOmega
        && !absent.segmenterTrace && !absent.segmenterIgnoreHostHints,"absent overrides retain dictionary/default selection");
    const auto set=ScopedEnvironment::set;
    set(finite.name.c_str(),"0.25trailing");set(positive.name.c_str(),"-3");set(integer.name.c_str(),"-2rest");
    set(bounded.name.c_str(),"0");set(boolean.name.c_str(),"oFf");set(debug.name.c_str(),"unexpected");
    set(path.name.c_str(),"relative snapshot");set(mode.name.c_str()," Skip ");
    set(gate.name.c_str(),"inf");set(omega.name.c_str(),"bad");
    set(trace.name.c_str(),"0");set(hint.name.c_str(),"false");
    const auto captured=readImplicitPicOverrides();
    check(captured.consoleDqrelWarn==.25 && captured.amrReinitActiveMovingMaxStreak==1
        && captured.segmenterMaxLevel==-2 && captured.segmenterMaxSegments==0,
        "legacy numeric prefixes and positive-only/raw-integer distinction");
    check(captured.orbitAvgFused==false && captured.seeDebug==true && captured.segmenterTrace
        && captured.segmenterIgnoreHostHints,"boolean, permissive SEE and presence-only segmenter flags remain distinct");
    check(captured.initGeomSnapshotPath=="relative snapshot" && captured.amrGeomRefreshMode==" Skip ",
        "owned strings preserve spaces; path stays relative to execution directory");
    check(std::isinf(*captured.amrGateEFieldRel) && captured.polarization.temporalOmega==.8,
        "legacy gate infinity and malformed omega fallback");
    for (const char* invalid:{"nan","inf","bad","1e999"}) {
        set(finite.name.c_str(),invalid);check(!readImplicitPicOverrides().consoleDqrelWarn,"nonfinite/invalid ordinary real retains fallback");
    }
    set(integer.name.c_str(),"999999999999999999999");set(bounded.name.c_str(),"bad");
    set(boolean.name.c_str()," true ");set(debug.name.c_str(),"NO");set(gate.name.c_str(),"nan");
    set(path.name.c_str(),nullptr);set(mode.name.c_str(),"adaptive");
    auto changed=readImplicitPicOverrides();
    check(!changed.segmenterMaxLevel && !changed.segmenterMaxSegments && !changed.orbitAvgFused
        && changed.seeDebug==false && changed.amrGateEFieldRel==0. && !changed.initGeomSnapshotPath,
        "malformed integers, untrimmed bool and NaN gate keep separate old policies");
    check(captured.initGeomSnapshotPath=="relative snapshot" && captured.amrGeomRefreshMode==" Skip "
        && captured.consoleDqrelWarn==.25,"captured values remain owned after process changes");
    PicCaseInput prepared;prepared.implicitOverrides=captured;
    auto copied=prepared;prepared.implicitOverrides.initGeomSnapshotPath="changed";
    check(copied.implicitOverrides.initGeomSnapshotPath=="relative snapshot","prepared input copy owns override storage");

    // Exercise the numerical segmenter with process flags still enabled. Only
    // explicit options may request trace/no-hint; a new environment cannot do so.
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=16;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    using Segmenter=ImplicitPIC::Orbit::PredictorOrbitSegmenter;
    std::vector<Vector2D> from{{.125,.25}},to{{.875,.25}};
    std::vector<int> offsets,hints{mesh.findCellLabel(from[0])};
    std::vector<ImplicitPIC::Orbit::OrbitSegment> segments;
    BoundaryLogs logs;Segmenter::Options options;
    Segmenter::segmentLinearOrbits(mesh,from,to,&hints,1.,offsets,segments,options);
    check(logs.output.str().empty() && logs.errors.str().empty() && segments.size()==2,
        "segmenter ignores live debug environment and crosses actual cells");
    const auto baseline=segments;
    Segmenter::segmentLinearOrbits(mesh,from,to,nullptr,1.,offsets,segments,options);
    check(segments.size()==baseline.size(),"ordinary lookup without an initial hint preserves segmentation");
    for(std::size_t i=0;i<segments.size();++i)
        check(segments[i].cellLabel==baseline[i].cellLabel && segments[i].t0==baseline[i].t0 && segments[i].t1==baseline[i].t1,
            "ordinary lookup preserves cells and crossing times");
}

void explicitAmrSetup()
{
    PicSpeciesState particles(2,123);particles.species().emplace_back("negative_ion",1.,-1.);
    particles.species().emplace_back("hot_ELECTRON",1.,1.);const auto rng=particles.rng();
    const auto rejected=[](auto action){try{action();}catch(const std::logic_error&){return true;}return false;};
    for(int present=0;present<4;++present)for(bool physics:{false,true})for(int flags=0;flags<4;++flags)for(int names=0;names<4;++names) {
        PicCaseInput input;auto& cfg=input.pic.amr;auto& mesh=input.mesh;
        cfg.enabled=true;cfg.debye.enabled=(flags&1)!=0;cfg.indicator.useEField=(flags&2)!=0;
        cfg.adaptEverySteps=13;cfg.indicatorEverySteps=11;cfg.persistenceSteps=9;cfg.hysteresis=.3;cfg.debye.Te_eV=0.;
        cfg.report.enabled=true;cfg.report.every=7;cfg.report.path="owned.csv";cfg.particleRelabel.mode="original";
        input.pic.species.resize(2);input.pic.species[0].Te_eV=3.;input.pic.species[1].Te_eV=7.;
        mesh.geometryRefinement.present=(present&1)!=0;mesh.physicsRefinement.present=(present&2)!=0;
        auto& p=mesh.physicsRefinement;p.enabled=physics;p.adaptEverySteps=names==0?-1:4;
        p.indicatorEverySteps=names==1?0:5;p.persistenceSteps=names==2?-3:6;p.hysteresis=names==3?-.2:0.;p.reportEvery=names==0?-7:2;
        if(names==1)p.indicators={"dEbYeLeNgTh","ELECTRICFIELD"};
        if(names==2)p.indicators={" DebyeLength ","unknown"};
        if(names==3)p.indicators={"ElectricField"};
        ExplicitPicAmrSetup setup(input);const auto& c=setup.config();const auto& refiner=setup.refinerMesh();
        const bool mapped=(present&2)!=0,active=present!=0 && physics;
        const bool debye=(flags&1)!=0 && (!mapped || (active && (names==0 || names==1)));
        const bool electric=(flags&2)!=0 && (!mapped || (active && names!=2));
        check(c.debye.enabled==debye && c.indicator.useEField==electric,"explicit indicator allow-list, presence and enablement gates");
        check(c.adaptEverySteps==(mapped && names!=0?4:13) && c.indicatorEverySteps==(mapped && names!=1?5:11)
            && c.persistenceSteps==(mapped && names!=2?6:9) && c.hysteresis==(mapped && names!=3?0.:.3),"explicit positive/zero/negative mesh overrides");
        check(c.report.every==(mapped?(names==0?0:2):7) && c.report.enabled==(!mapped || names!=0),"explicit report override");
        check(refiner.physicsRefinement.enabled==(present==0 || physics),"legacy no-block refiner policy");
        if(present==0){std::vector<std::string> expected;if(flags&1)expected.push_back("DebyeLength");if(flags&2)expected.push_back("ElectricField");
            check(refiner.physicsRefinement.indicators==expected && refiner.physicsRefinement.hysteresis==.3,"legacy refiner indicators/order");}
        else check(refiner.physicsRefinement.indicators==p.indicators && refiner.physicsRefinement.hysteresis==p.hysteresis,"declared mesh policy retained");
        input.pic.amr.report.path="changed";input.pic.amr.adaptEverySteps=999;input.mesh.physicsRefinement.enabled=!physics;
        check(c.report.path=="owned.csv" && c.particleRelabel.mode=="original" && c.adaptEverySteps!=999,"explicit setup must own copied inputs");
        std::ostringstream errors;check(rejected([&]{setup.selectElectron(particles,errors);}) && errors.str().empty(),"particle gate requires mesh activation");
        setup.activateForMesh(2,false,errors);check(setup.enabled()==active,"AMR mesh stage before electron/indicator gate");
        check(rejected([&]{setup.temperatureEV(input.pic.species);}),"temperature cannot precede electron selection");
        setup.selectElectron(particles,errors);
        check(setup.enabled()==(active && (debye||electric)) && setup.electronIndex()==(active && debye?1:-1),"live named electron precedence and inactive policy");
        check(setup.temperatureEV(input.pic.species)==(active && debye?7.:0.),"explicit configured electron temperature fallback");
        check(rejected([&]{setup.activateForMesh(0,true,errors);}) && rejected([&]{setup.selectElectron(particles,errors);}),"explicit setup phases cannot repeat");
    }
    PicSpeciesState none(1,4);none.species().emplace_back("ion",1.,1.);
    PicCaseInput input;input.mesh.physicsRefinement.present=true;input.mesh.physicsRefinement.enabled=true;
    input.pic.amr.indicator.useEField=false;input.pic.amr.debye.enabled=true;input.pic.amr.debye.Te_eV=-1.;
    std::ostringstream errors;ExplicitPicAmrSetup noLevels(input);noLevels.activateForMesh(0,true,errors);noLevels.selectElectron(particles,errors);
    check(!noLevels.enabled() && noLevels.electronIndex()==-1 && noLevels.temperatureEV({})==-1.
        && errors.str()=="[ExplicitPIC][AMR] WARNING: AMR enabled but mesh.maxLevel=0. Check mesh.dict adaptive settings (maxLevel_local). Disabling AMR.\n","original zero-level warning and inactive fallback");
    for(bool geometry:{false,true}) {
        errors.str("");ExplicitPicAmrSetup missing(input);missing.activateForMesh(2,geometry,errors);missing.selectElectron(none,errors);
        const std::string expected="[ExplicitPIC][AMR] WARNING: no electron species found; disabling Debye indicator.\n";
        check(errors.str()==expected+(geometry?"":"[ExplicitPIC][AMR] WARNING: no active refinement indicators; disabling AMR.\n")
            && missing.enabled()==geometry && !missing.config().debye.enabled,"missing electron warnings and geometry-only survival");
    }
    for(double temperature:{-1.,0.,2.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        input.pic.amr.debye.Te_eV=temperature;ExplicitPicAmrSetup setup(input);setup.activateForMesh(2,false,errors);setup.selectElectron(particles,errors);
        const double actual=setup.temperatureEV({});check(std::isnan(temperature)?std::isnan(actual):actual==(temperature<=0.?2.:temperature),"missing configured species and original nonfinite temperature behavior");
        std::vector<SpeciesConfig> values(2);values[1].Te_eV=std::numeric_limits<double>::quiet_NaN();
        if(temperature<=0.)check(std::isnan(setup.temperatureEV(values)),"configured NaN temperature must not be silently repaired");
    }
    check(particles.rng()==rng && particles.species()[0].charge()==-1. && particles.species()[1].name()=="hot_ELECTRON","AMR setup must not modify particles or RNG");
}
void implicitAmrSetup()
{
    using Mode=PicAmrGeometry::RefreshMode;
    ScopedEnvironment liveMode("LAMP_IMPLICIT_AMR_GEOM_REFRESH_MODE","skip"),
        liveRemap("LAMP_IMPLICIT_AMR_FAST_REMAP","1");
    PicCaseInput input;input.pic.amr.debye.Te_eV=0.;input.pic.amr.debye.CD=0.;
    auto setup=resolveImplicitPicAmr(input,-1,false);
    check(!setup.indicators.enabled && !setup.indicators.wantDebye && !setup.indicators.wantEField
        && setup.adaptEvery==0 && setup.indicators.indicatorEvery==1 && setup.indicators.persistenceSteps==1,
        "implicit AMR absent blocks and cadence defaults");
    check(setup.indicators.temperatureEV==2. && setup.indicators.debyeScale==.5
        && setup.indicators.electricPercentile==95. && !setup.fastRemap && !setup.resetTransients
        && setup.geometry.mode==Mode::Full,"prepared AMR setup ignores live environment");
    std::ostringstream output,errors;setup.writeReport(output,errors,true);
    check(errors.str().empty() && output.str()=="[ImplicitPIC][AMR] geometry refresh mode=full relTh=0.2 allowWithMovingEB=0 movingEbRefreshMode=full reinitActiveMovingOnly=0 hybridSensitiveInactive=0 hybridUseDelta=1 reinitActiveMovingMaxStreak=0 (dict: built-in defaults via pic.dict/AMR/GeometryRefresh, env override supported)\n",
        "default AMR report exact format/source");
    for(int present=0;present<4;++present) {
        input.mesh.geometryRefinement.present=(present&1)!=0;input.mesh.physicsRefinement.present=(present&2)!=0;
        input.mesh.physicsRefinement.enabled=true;input.mesh.physicsRefinement.indicators={"DEBYELENGTH","ElectricField","ignored"};
        setup=resolveImplicitPicAmr(input,-1,false);
        check(setup.indicators.enabled==(present!=0) && setup.indicators.wantDebye==(present!=0)
            && setup.indicators.wantEField==(present!=0),"legacy combined mesh-block presence gate");
    }
    input.mesh.physicsRefinement.enabled=false;setup=resolveImplicitPicAmr(input,9,true);
    check(setup.indicators.enabled && !setup.indicators.wantDebye && !setup.indicators.wantEField
        && setup.indicators.electronIndex==9 && setup.indicators.temperatureEV==2.,"geometry-only activation and invalid electron fallback");
    input.mesh.physicsRefinement.enabled=true;input.mesh.physicsRefinement.adaptEverySteps=4;
    input.mesh.physicsRefinement.indicatorEverySteps=-1;input.mesh.physicsRefinement.persistenceSteps=-7;
    input.mesh.physicsRefinement.indicators={" DebyeLength ","ElectricField"};
    input.pic.species.resize(2);input.pic.species[0].Te_eV=0.;input.pic.species[1].Te_eV=7.;
    setup=resolveImplicitPicAmr(input,1,false);
    check(setup.adaptEvery==4 && setup.indicators.indicatorEvery==4 && setup.indicators.persistenceSteps==1
        && !setup.indicators.wantDebye && setup.indicators.wantEField && setup.indicators.temperatureEV==7.,
        "cadence fallback, species temperature and untrimmed indicator names");
    check(resolveImplicitPicAmr(input,0,false).indicators.temperatureEV==1e-6,"electron temperature floor");
    input.mesh.physicsRefinement.indicatorEverySteps=3;input.mesh.physicsRefinement.persistenceSteps=2;
    input.pic.amr.debye.Te_eV=11.;input.pic.amr.debye.CD=.8;input.pic.amr.indicator.eFieldPercentile=101.;
    setup=resolveImplicitPicAmr(input,0,false);
    check(setup.indicators.temperatureEV==11. && setup.indicators.debyeScale==.8 && setup.indicators.electricPercentile==101.
        && setup.indicators.indicatorEvery==3 && setup.indicators.persistenceSteps==2,"explicit reference values retain old bounds");
    for(double invalid:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        input.pic.amr.indicator.eFieldPercentile=invalid;
        check(resolveImplicitPicAmr(input,0,false).indicators.electricPercentile==95.,"nonfinite/negative percentile fallback");
    }
    input.pic.amr.indicator.eFieldPercentile=0.;
    check(resolveImplicitPicAmr(input,0,false).indicators.electricPercentile==0.,"zero percentile is valid");

    auto& mesh=input.mesh.geometryRefresh;auto& pic=input.pic.amr.geometryRefresh;
    pic.present=true;pic.mode="SKIP";pic.changeRel=.7;pic.verbose=true;pic.allowWithMovingEb=true;pic.movingEbRefreshMode="EXACT_REFRESH";
    setup=resolveImplicitPicAmr(input,0,false);output.str("");errors.str("");setup.writeReport(output,errors,true);
    check(setup.geometry.mode==Mode::Skip && setup.geometry.amrGeomRefreshChangeRelTh==.7
        && setup.geometry.amrGeomRefreshVerbose && setup.geometry.amrGeomRefreshAllowWithMovingEb,
        "deprecated PIC-only geometry policy remains effective");
    check(output.str().find("move this block to mesh.dict/GeometryRefresh.")!=std::string::npos
        && output.str().find("(dict: pic.dict/AMR/GeometryRefresh,")!=std::string::npos
        && errors.str().find("requires LAMP_IMPLICIT_AMR_ALLOW_EXACT_MOVINGEB_REFRESH=1.")!=std::string::npos,
        "deprecated source and exact-refresh compatibility warning");
    mesh.present=true;mesh.mode="adaptive";mesh.changeRel=.3;
    setup=resolveImplicitPicAmr(input,0,false);output.str("");errors.str("");setup.writeReport(output,errors,true);
    check(setup.geometry.mode==Mode::Adaptive && setup.geometry.amrGeomRefreshChangeRelTh==.3
        && !setup.geometry.amrGeomRefreshVerbose && !setup.geometry.amrGeomRefreshAllowWithMovingEb,
        "mesh dictionary wins every geometry policy field");
    check(errors.str()=="[ImplicitPIC][AMR] WARNING: GeometryRefresh is defined in both mesh.dict and pic.dict with different values; runtime is using mesh.dict/GeometryRefresh.\n"
        && output.str().find("runtime GeometryRefresh now follows mesh.dict/GeometryRefresh.")!=std::string::npos,
        "conflict warning exact format and discarded PIC exact request");
    mesh.mode=" A D A P T I V E ";pic.mode="adaptive";pic.changeRel=.3+5e-13;
    pic.verbose=false;pic.allowWithMovingEb=false;pic.movingEbRefreshMode="full";
    setup=resolveImplicitPicAmr(input,0,false);output.str("");errors.str("");setup.writeReport(output,errors,true);
    check(setup.geometry.mode==Mode::Full && errors.str().empty(),"conflict whitespace normalization differs from runtime token matching");
    pic.changeRel=.3+2e-12;setup=resolveImplicitPicAmr(input,0,false);
    output.str("");errors.str("");setup.writeReport(output,errors,false);
    check(output.str().empty() && !errors.str().empty(),"conflict tolerance and warnings independent of info level");
    auto& overrides=input.implicitOverrides;
    overrides.amrFastRemap=true;overrides.resetTransientsOnAmr=true;overrides.amrMaxChanges=-8;
    overrides.amrGateEFieldRel=-.25;overrides.amrGateDebyeRel=std::numeric_limits<double>::infinity();
    overrides.amrGeomRefreshMode="SKIP";overrides.amrGeomRefreshChangeRel=-1.;
    overrides.amrGeomRefreshVerbose=true;overrides.amrGeomRefreshAllowWithMovingEb=true;
    overrides.amrMovingEbRefreshMode="EXACT_REFRESH";overrides.amrAllowExactMovingEbRefresh=true;
    overrides.amrReinitActiveMovingOnly=true;overrides.amrReinitHybridSensitiveInactive=true;
    overrides.amrReinitHybridUseDelta=false;overrides.amrReinitHybridIncludeMaxIter=false;
    overrides.amrReinitHybridNearZeroMin=-2;overrides.amrReinitHybridNarrowBandMin=3;
    overrides.amrReinitHybridNearZeroDeltaMin=4;overrides.amrReinitHybridNarrowBandDeltaMin=-5;
    overrides.amrReinitHybridMaxIterThreshold=0;overrides.amrReinitActiveMovingMaxStreak=-1;
    const auto captured=resolveImplicitPicAmr(input,1,false);const auto& geometry=captured.geometry;
    check(captured.fastRemap && captured.resetTransients && captured.indicators.maxChanges==0
        && captured.indicators.gateElectricRelative==-.25 && std::isinf(captured.indicators.gateDebyeRelative),
        "typed gates retain caller values while maximum changes clamps");
    check(geometry.mode==Mode::Skip && geometry.amrGeomRefreshChangeRelTh==0. && geometry.amrGeomRefreshVerbose
        && geometry.amrGeomRefreshAllowWithMovingEb && geometry.amrReinitActiveMovingOnly && geometry.amrReinitHybridSensitiveInactive
        && !geometry.amrReinitHybridUseDelta && !geometry.amrReinitHybridIncludeMaxIter
        && geometry.amrReinitHybridNearZeroMin==0 && geometry.amrReinitHybridNarrowBandMin==3
        && geometry.amrReinitHybridNearZeroDeltaMin==4 && geometry.amrReinitHybridNarrowBandDeltaMin==0
        && geometry.amrReinitHybridMaxIterThreshold==1 && geometry.amrReinitActiveMovingMaxStreak==0,
        "all explicit geometry overrides and independent lower bounds");
    input.implicitOverrides={};mesh.present=pic.present=false;mesh.mode="full";pic.changeRel=mesh.changeRel;
    output.str("");errors.str("");captured.writeReport(output,errors,true);
    check(output.str().find("fast remap enabled")!=std::string::npos && output.str().find("movingEbRefreshMode=exact_refresh")!=std::string::npos
        && output.str().find("(dict: mesh.dict/GeometryRefresh,")!=std::string::npos
        && errors.str().find("requires LAMP_IMPLICIT") == std::string::npos,
        "resolved setup owns its input/provenance snapshot and exact permission affects report only");
    output.str("");errors.str("");captured.writeReport(output,errors,false);
    check(output.str()=="[ImplicitPIC][AMR] fast remap enabled (LAMP_IMPLICIT_AMR_FAST_REMAP=1)\n"
        && !errors.str().empty(),"fast-remap and conflict messages retain independent info gates");
}
void segmenterLookup()
{
    using Segmenter=ImplicitPIC::Orbit::PredictorOrbitSegmenter;
    using Segment=ImplicitPIC::Orbit::OrbitSegment;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=256;input.maxInitialIterations=0;
    input.maxLevel_local={2};QuadtreeMesh mesh(input);
    std::vector<Vector2D> from{{.125,.2},{.875,.2},{.2,.125},{.2,.875},{.125,.125},{.3,.2},{-.1,.2}},
                          to{{.875,.2},{.125,.2},{.2,.875},{.2,.125},{.875,.875},{1.2,.2},{.3,.2}};
    Segmenter::Options options;options.maxSegmentsPerParticle=16;
    const auto verify=[&](const Segmenter::Options& opt) {
        std::vector<int> expectedOffsets,hints;
        for(auto p:from)hints.push_back(mesh.findCellLabel(p));
        std::vector<Segment> expected;
        Segmenter::segmentLinearOrbits(mesh,from,to,&hints,1.,expectedOffsets,expected,opt);
        for(int mode=0;mode<5;++mode) {
            // Correct, wrong-but-live, retired/unused, out-of-range and absent.
            const int invalid=mode==1?mesh.findCellLabel({.8,.8}):mode==2?255:9999;
            if(mode>0)hints.assign(from.size(),invalid);
            auto independent=opt;independent.ignoreHostHints=true;
            std::vector<int> offsets;std::vector<Segment> actual;
            Segmenter::segmentLinearOrbits(mesh,from,to,mode==4?nullptr:&hints,1.,offsets,actual,independent);
            check(offsets==expectedOffsets && actual.size()==expected.size(),"position lookup changed orbit segmentation");
            for(std::size_t i=0;i<actual.size();++i) {
                const auto& a=actual[i];const auto& b=expected[i];
                check(a.cellLabel==b.cellLabel && a.t0==b.t0 && a.t1==b.t1
                    && a.x0.x==b.x0.x && a.x0.y==b.x0.y && a.x1.x==b.x1.x && a.x1.y==b.x1.y,
                    "ignored hints changed segment cells, positions or times");
                check(a.t1>a.t0 && mesh.isInsideCell(a.cellLabel,(a.x0+a.x1)*.5,0.),
                    "segment interior is outside its reported host");
            }
        }
        return expected;
    };
    auto result=verify(options);
    check(result[0].cellLabel==0 && result[1].cellLabel==1,"analytic horizontal host sequence");
    near(result[0].t1,.5,"analytic first crossing time");near(result[1].t0,.5,"analytic second segment starts at crossing");
    near(result[1].t1,1.,"complete orbit duration");
    options.periodicX=options.periodicY=true;verify(options); // Existing wrapped path and bisection policy.
    options.periodicX=options.periodicY=false;
    DictionaryMesh::RefineBox box;box.min={.51,0.};box.max={1.,1.};box.level=2;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()>4,"lookup test needs coarse/fine interfaces");
    verify(options);options.segmentCellMaxLevel=0;verify(options);options.segmentCellMaxLevel=-1;
    // A planar wall in the refined half must stop the path at x=0.7.
    VertField<double> phi(mesh,0.);
    for(int i=0;i<phi.size();++i)phi.value[i]=.7-mesh.verts(phi.label(i)).coord().x;
    std::vector<double> alpha(mesh.cellIndexNum());
    for(int i=0;i<mesh.cellIndexNum();++i) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);const auto width=cell.size().x;
        alpha[i]=std::clamp((.7-(cell.centerCoord().x-.5*width))/width,0.,1.);
    }
    from={{.125,.2},{.9,.2}};to={{.875,.2},{.95,.2}};
    options.enableEmbeddedIntercept=true;options.embeddedLevelSetVert=&phi;options.alphaByCellIndex=&alpha;
    result=verify(options);check(!result.empty(),"wall interception emitted no fluid segment");
    near(result.back().x1.x,.7,"analytic embedded-wall endpoint");
    near(result.back().t1,(.7-.125)/.75,"analytic embedded-wall hit time");
    options={};options.ignoreHostHints=true;
    from.assign(96,{.125,.2});to.assign(96,{.875,.2});
    // Compare the OpenMP chunk path with independent one-particle calls.
    std::vector<int> offsets;std::vector<Segment> batch,single;
    Segmenter::segmentLinearOrbits(mesh,from,to,nullptr,1.,offsets,batch,options);
    check(offsets.size()==97 && offsets.back()==static_cast<int>(batch.size()),"parallel orbit CSR shape");
    std::vector<int> one;
    Segmenter::segmentLinearOrbits(mesh,{from[0]},{to[0]},nullptr,1.,one,single,options);
    for(int p=0;p<96;++p) {
        check(offsets[p+1]-offsets[p]==static_cast<int>(single.size()),"parallel lookup dropped an orbit");
        for(std::size_t i=0;i<single.size();++i) {
            const auto& a=batch[offsets[p]+i];const auto& b=single[i];
            check(a.cellLabel==b.cellLabel && a.t0==b.t0 && a.t1==b.t1,"parallel lookup changed crossing sequence");
        }
    }
    options.maxSegmentsPerParticle=1;bool rejected=false;
    try { Segmenter::segmentLinearOrbits(mesh,{from[0]},{to[0]},nullptr,1.,one,single,options); }
    catch(const std::runtime_error& e){rejected=std::string(e.what()).find("exceeded maxSegmentsPerParticle")!=std::string::npos;}
    check(rejected,"position lookup bypassed the real segment budget");
}
void speciesRunStartup(const std::filesystem::path& directory)
{
    using State=PicSpeciesState;using Mode=PicDomain::Mode;
    struct ConsoleScope {ConsoleControl::Config previous=ConsoleControl::g_config;
        ~ConsoleScope(){ConsoleControl::g_config=previous;}} console;
    ConsoleControl::setConfig("debug",1);BoundaryLogs logs;
    const auto bytes=[](const std::filesystem::path& p){std::ifstream f(p,std::ios::binary);
        check(bool(f),"particle startup cache missing");return std::string(std::istreambuf_iterator<char>(f),{});};
    const auto rejected=[](const auto& action){try{action();}catch(const std::logic_error&){return true;}return false;};
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=64;input.maxInitialIterations=0;
    NormalizationScales scales;scales.L0_m=2.;scales.T0_s=.5;scales.Phi0_V=3.;scales.v0_mps=4.;
    std::vector<SpeciesConfig> config(3);
    config[0].name="electron";config[0].q=-2.*PhysicalConstants::eV_J;config[0].m=PhysicalConstants::eV_J;
    config[0].weight=3.;config[0].numParticles=31;config[0].Te_eV=8.;
    config[1]=config[0];config[1].name="ion";config[1].q=-config[0].q;config[1].weight=5.;config[1].numParticles=17;config[1].Te_eV=0.;
    config[2]=config[1];config[2].name="neutral";config[2].q=0.;config[2].numParticles=0;
    for(Mode mode:{Mode::Explicit,Mode::Implicit})for(bool wall:{false,true}) {
        const bool implicit=mode==Mode::Implicit;
        const auto dir=directory/(std::string(implicit?"species_implicit":"species_explicit")+(wall?"_wall":"_plain"));
        std::filesystem::create_directories(dir);
        PicDomain domain(input,input,mode);auto& mesh=domain.mesh();
        DictionaryGeometry dict;Line line;line.points[0]={.5,0.};line.points[1]={.5,1.};line.setType(LineType::straight);
        dict.lines={line};dict.levelSetName={"wall"};dict.innerPoints={{.25,.5}};dict.moving={false};dict.linesOnInterface={{{0,false}}};
        bool eb=wall,refine=false;check(domain.initializeGeometry(wall?&dict:nullptr,eb,refine,{}),"particle startup domain");
        if(wall) {
            auto* geometry=domain.geometry();geometry->levelSets[0].initialize(mesh,0.);
            for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry->levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.5;
            geometry->refreshLevelSetDerivatives(mesh);domain.embedded().rebuildComposite(mesh,geometry);
        }
        domain.rebuildFaceGraph();
        if(wall)domain.embedded().rebuildMetrics(mesh,domain.geometry(),domain.faceGraph(),{1});
        const auto same=[&](State& actual,State& expected) {
            check(actual.rng()==expected.rng() && actual.weights()==expected.weights() && actual.ions()==expected.ions()
                && actual.electrons()==expected.electrons() && actual.neutrals()==expected.neutrals(),"startup preserves master RNG and configured metadata");
            for(const auto& pair:{std::make_pair(&actual,"actual"),std::make_pair(&expected,"expected")}) {
                const auto path=dir/pair.second;std::filesystem::create_directories(path);
                for(const auto& sp:pair.first->species())sp.writeCache(PicCheckpointStore::pathPrefix(path));
            }
            for(std::size_t si=0;si<config.size();++si) {
                const auto file="particle_"+config[si].name+"_cache.bin";
                check(bytes(dir/"actual"/file)==bytes(dir/"expected"/file),"startup exact particle payload, metadata and host labels");
                const auto& a=actual.species()[si];const auto& b=expected.species()[si];
                check(a.number==b.number && a.particle.size()==b.particle.size(),"startup species counts");
                for(std::size_t pi=0;pi<a.particle.size();++pi)check(a.particle[pi].trackId==b.particle[pi].trackId,"startup/cull preserves runtime IDs");
            }
        };
        const auto original=[&](State& state,const PicCheckpointReader* restart,bool restoreRng) {
            state.initialize(config,scales,mesh,restart!=nullptr,[&](const Vector2D& position,int cellLabel) {
                return (implicit || (domain.embedded().ready() && domain.embedded().absorbing())) &&
                    domain.embedded().containsPoint(mesh,domain.geometry(),position,cellLabel);
            });
            if(restart && !restart->readParticles(mesh,state.species(),state.rng(),restoreRng))return false;
            ParticleManager::cullEmbedded(mesh,domain.geometry(),domain.embedded(),state.species(),state.weights(),implicit);
            return true;
        };
        for(bool enabled:{false,true})for(bool crossing:{false,true})for(bool periodic:{false,true}) {
            State policy(0,79);const auto rng=policy.rng();
            check(rejected([&]{(void)policy.pushOptions();}),"unconfigured push policy cannot be consumed");
            policy.configurePush(domain,enabled,periodic,!periodic,crossing);
            const auto& options=policy.pushOptions();const auto& embedded=domain.embedded();
            const bool absorb=(implicit?enabled:true) && embedded.ready() && embedded.absorbing();
            const bool cross=(implicit?(enabled && embedded.ready()):absorb) && crossing;
            const State::PushOptions defaults;
            const auto indices=implicit || absorb?embedded.levelSetIndices():defaults.embeddedLevelSetIds;
            check(options.absorbAtDomain && options.absorbAtEmbedded==absorb && options.absorbEmbeddedCrossing==cross
                && options.periodicX==periodic && options.periodicY==!periodic,
                "owned policy preserves explicit/implicit readiness, enable, periodic and crossing gates");
            check(options.embeddedInsideIsNegative==(implicit || absorb?embedded.insideIsNegative():defaults.embeddedInsideIsNegative)
                && options.embeddedLevelSetIds==indices,"policy retains mode-specific embedded selection and sign");
            check(policy.species().empty() && policy.rng()==rng,"configuring push never samples or changes particles/RNG");
            check(rejected([&]{policy.configurePush(defaults);}) && &policy.pushOptions()==&options
                && options.periodicX==periodic,"reconfiguration rejects without replacing the owned snapshot");
        }
        {
            State copied(0,17);State::PushOptions options;options.embeddedLevelSetIds={3,7};options.periodicX=true;
            copied.configurePush(options);options.embeddedLevelSetIds={9};options.periodicX=false;
            check(copied.pushOptions().embeddedLevelSetIds==std::vector<int>({3,7}) && copied.pushOptions().periodicX,
                "custom push policy owns its vector and scalar inputs");
            State missing(0,19);missing.initialize({},scales,mesh,false);const auto rng=missing.rng();
            CellField<Vector2D> field(mesh,{0.,0.});
            check(rejected([&]{missing.beginPush();}) && rejected([&]{missing.initializeBackwardHalfStep(mesh,{},field,0.,nullptr);})
                && missing.rng()==rng,"unconfigured push and half-kick reject before advancing state");
            missing.configurePush(options);missing.initializeBackwardHalfStep(mesh,{},field,0.,nullptr);
            check(rejected([&]{missing.initializeBackwardHalfStep(mesh,{},field,0.,nullptr);}),
                "missing-policy rejection preserves the single startup half-kick for later configuration");
        }

        State actual(3,77),expected(3,77);
        const auto first=actual.initializeForRun(config,scales,domain,nullptr,true);
        check(first && *first==0 && original(expected,nullptr,true),"fresh startup returns step zero");same(actual,expected);
        check(actual.species()[0].number>0 && actual.species()[0].number==(wall?expected.species()[0].number:31)
            && (!wall || actual.species()[0].number<31),"initial exclusion rejects samples without resampling");
        check(actual.index().empty() && expected.index().empty(),"startup leaves name indexing at its later driver stage");
        actual.rebuildIndex();check(actual.index().at("neutral")==2 && actual.electronIndex()==0,"later name consumers use current species");
        const auto savedRng=actual.rng();
        check(rejected([&]{actual.initializeForRun(config,scales,domain,nullptr,true);}) && actual.rng()==savedRng,"repeated startup rejects before resampling");
        std::vector<CouplingWeightsCache> emptyWeights(3);CellField<Vector2D> zeroField(mesh,{0.,0.});
        for(std::size_t si=0;si<3;++si){const auto n=actual.species()[si].particle.size();emptyWeights[si].assignCSR(n,std::vector<int>(n+1,0),{},{});}
        State::PushOptions options;options.absorbAtDomain=false;options.absorbAtEmbedded=false;
        actual.configurePush(options);
        actual.initializeBackwardHalfStep(mesh,emptyWeights,zeroField,0.,nullptr);
        check(rejected([&]{actual.initializeBackwardHalfStep(mesh,emptyWeights,zeroField,0.,nullptr);}),"fresh run retains exactly one backward half-step");

        State seed(3,9123);seed.initialize(config,scales,mesh,true);
        for(std::size_t si=0;si<2;++si)for(int p=0;p<4;++p) {
            Particle value;value.position={p%2?.75:.25,.25+.5*(p/2)};value.velocity={1.+p,-2.};value.vz=3.;
            value.w=p%2?0.:7.;value.inCellL=-1;value.inBackgroundCellL=-1;value.rngState=100+p;value.trackId=900+p;
            seed.species()[si].particle.push_back(value);seed.species()[si].number=4;
        }
        seed.rng().discard(13);
        for(int number:{2,3,4}) {
            const auto checkpoint=dir/("pic."+std::to_string(number));std::filesystem::create_directories(checkpoint);
            check(PicCheckpointStore::writeMetadata(checkpoint,number,.1*number,.1),"particle checkpoint metadata");
            if(number!=3)check(PicCheckpointStore::writeRng(checkpoint,seed.rng()),"particle checkpoint RNG");
            for(std::size_t si=0;si<(number==4?1u:3u);++si)seed.species()[si].writeCache(PicCheckpointStore::pathPrefix(checkpoint));
        }
        PicCheckpointStore store(dir,{});double dt=.1;
        for(int number:{2,3})for(bool restoreRng:{false,true}) {
            auto reader=PicCheckpointReader::select(store,number,dt,implicit?"ImplicitPIC":"ExplicitPIC");check(reader.has_value(),"particle checkpoint selected");
            State restored(3,77),reference(3,77);auto step=restored.initializeForRun(config,scales,domain,&*reader,restoreRng);
            check(step && *step==number && original(reference,&*reader,restoreRng),"restored loop step published after all particle reads");same(restored,reference);
            check(restored.rng()==(restoreRng && number==2?seed.rng():State::Engine(77)),"RNG restore/disabled/missing-file policies");
            check(restored.species()[0].number==(wall?2:4) && restored.index().empty(),"restored cull runs before late indexing");
            check(rejected([&]{restored.initializeBackwardHalfStep(mesh,{},zeroField,0.,nullptr);}),"restored particles cannot receive startup half-kick");
            for(const auto& sp:restored.species())for(const auto& p:sp.particle)
                check(p.trackId==0 && p.inCellL==mesh.findCellLabel(p.position),"restore repairs hosts without assigning new trajectory IDs");
        }
        auto incomplete=PicCheckpointReader::select(store,4,dt,"PIC");check(incomplete.has_value(),"legacy specified checkpoint selection");
        State failed(3,77);check(!failed.initializeForRun(config,scales,domain,&*incomplete,true),"missing later species fails startup");
        check(failed.rng()==seed.rng() && failed.species()[0].number==4 && failed.species()[1].particle.empty()
            && failed.species()[0].particle[0].position.x==.25,"failed restore preserves earlier RNG/particle reads and does not prematurely cull");
        check(rejected([&]{failed.initializeForRun(config,scales,domain,nullptr,false);}),"partial restore cannot be initialized again");
        if(wall)for(bool selected:{false,true}) {
            State culled(3,9),manual(3,9);culled.initialize(config,scales,mesh,true);manual.initialize(config,scales,mesh,true);
            culled.species()=seed.species();manual.species()=seed.species();
            std::vector<unsigned char> mask{1};std::vector<long long> loss{10};std::vector<double> physical{20.};
            auto expectedLoss=loss;auto expectedPhysical=physical;
            const auto result=ParticleManager::cullEmbedded(mesh,domain.geometry(),domain.embedded(),manual.species(),manual.weights(),implicit,
                selected?&mask:nullptr,&expectedLoss,&expectedPhysical);
            const auto removed=culled.cullEmbedded(domain,"contract",selected?&mask:nullptr,&loss,&physical);
            check(removed==result.removed && removed==(selected?2u:4u) && loss==expectedLoss && physical==expectedPhysical,"runtime cull keeps masks and bounded loss accumulation");
            same(culled,manual);
            check(implicit?culled.species()[0].particle[0].inCellL>=0:culled.species()[0].particle[0].inCellL==-1,"only implicit cull refreshes survivor hosts");
        }
        State uninitialized(0,9);check(rejected([&]{uninitialized.cullEmbedded(domain,"invalid");}),"uninitialized cull rejected");
    }
    check(logs.output.str().find("[ExplicitPIC][EB] Culled")!=std::string::npos
        && logs.output.str().find("[ImplicitPIC][EB] Culled")!=std::string::npos
        && logs.errors.str().find("RNG state not found")!=std::string::npos,"startup retains mode-specific cull and missing RNG reports");
}
void speciesInitializationAndRelabel()
{
    using State=PicSpeciesState;using Pusher=ExplicitPIC::Pusher::ParticlePusherLeapfrog;
    const auto same=[](const Particle& a,const Particle& b) {
        check(a.position.x==b.position.x && a.position.y==b.position.y && a.velocity.x==b.velocity.x
            && a.velocity.y==b.velocity.y && a.vz==b.vz && a.w==b.w && a.inCellL==b.inCellL
            && a.inBackgroundCellL==b.inBackgroundCellL && a.discreteVertL==b.discreteVertL
            && a.emitChance==b.emitChance && a.rngState==b.rngState && a.trackId==b.trackId
            && a.levelSetValues==b.levelSetValues,"particle transition differs from original operation");
        for(int i=0;i<4;++i)check(a.weight[i]==b.weight[i],"particle interpolation history changed");
    };
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=256;input.maxInitialIterations=0;
    const auto norm=NormalizationContext::MakeLegacy();std::vector<SpeciesConfig> config(3);
    for(int i=0;i<3;++i) {config[i].name="initial_"+std::to_string(i);config[i].m=1.;config[i].q=i?0.:1.;config[i].numParticles=0;}
    for(int variant=0;variant<3;++variant) {
        QuadtreeMesh mesh(input);CellField<Vector2D> field(mesh,Vector2D(.5,-.25));
        DictionaryGeometry dict;dict.points={{.5,0.},{.5,1.}};Line wall;
        wall.points[0]={.5,0.};wall.points[1]={.5,1.};wall.setType(LineType::straight);wall.levelSetIndex=0;
        dict.lines={wall};dict.levelSetName={"wall"};dict.innerPoints={{.75,.5}};
        dict.linesOnInterface={{{0,false}}};dict.moving={false};
        ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh);
        for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry.levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.5;
        geometry.refreshLevelSetDerivatives(mesh);
        auto* solid=variant==1?&geometry:nullptr;
        State state(3,77);state.initialize(config,norm.scales(),mesh,false);
        auto& species=state.species();species[0].set_kappa_s(2.);species[1].set_kappa_s(0.);
        const int count=variant==2?20016:5;
        const std::vector<Vector2D> points{{.25,.2},{.75,.7},{-.1,.2},{1.1,.8},{.4,.3}};
        for(int i=0;i<count;++i) {
            Particle p;p.position=variant==2?Vector2D(i%2?-.1:.25,.2):points[i];p.velocity={1.,-2.};
            p.vz=3.;p.w=2.;p.rngState=100+i;p.trackId=200+i;p.discreteVertL=4;
            p.emitChance=.25;p.weight[0]=.4;p.levelSetValues={-2.,1.};
            p.inCellL=i%2?-1:mesh.findCellLabel(p.position);p.inBackgroundCellL=-1;
            species[0].particle.push_back(p);
        }
        species[1].particle.push_back(species[0].particle[0]);species[0].number=-1;
        std::vector<CouplingWeightsCache> caches(3);
        for(std::size_t si=0;si<species.size();++si) {
            std::vector<int> offsets,indices;std::vector<double> weights;
            for(std::size_t p=0;p<species[si].particle.size();++p) {
                offsets.push_back(static_cast<int>(indices.size()));indices.push_back(0);indices.push_back(1);
                weights.push_back(.25);weights.push_back(.75);
            }
            offsets.push_back(static_cast<int>(indices.size()));caches[si].assignCSR(species[si].particle.size(),offsets,indices,weights);
        }
        const auto initial=species;auto expected=species;const auto rng=state.rng();Pusher::Options options;
        options.absorbAtEmbedded=true;options.embeddedInsideIsNegative=false;options.embeddedLevelSetIds={0};
        options.periodicX=true; // Original half-kick absorbs outside; it does not perform the full push's wrap.
        state.configurePush(options);
        auto invalid=caches;invalid[1].assignCSR(2,{0,1,2},{0,1},{1.,1.});bool rejected=false;
        try {state.initializeBackwardHalfStep(mesh,invalid,field,.1,solid);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected && species[0].number==-1,"late-species mismatch rejects before the first kick");
        same(species[0].particle[0],initial[0].particle[0]);
        rejected=false;try {state.initializeBackwardHalfStep(mesh,{},field,.1,solid);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"missing half-step species caches rejected");
        CellField<Vector2D> stale;
        rejected=false;try {state.initializeBackwardHalfStep(mesh,caches,stale,.1,solid);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"half-step field shape rejected");
        for(std::size_t si=0;si<expected.size();++si) {
            if(expected[si].particle.empty())continue;std::vector<Vector2D> gathered;
            ExplicitPIC::Field::MatchedGather::gatherVector2D(caches[si],field,gathered,false);
            Pusher::initializeHalfStep(mesh,expected[si],gathered,-.1,options,solid);
        }
        state.initializeBackwardHalfStep(mesh,caches,field,.1,solid);
        for(std::size_t si=0;si<species.size();++si) {
            check(species[si].number==expected[si].number && species[si].particle.size()==expected[si].particle.size(),"half-kick survivor count");
            for(std::size_t pi=0;pi<species[si].particle.size();++pi)same(species[si].particle[pi],expected[si].particle[pi]);
        }
        check(species[0].particle.size()==(variant==2?10008u:(variant==1?2u:3u)),"half-step domain/embedded compaction exercised");
        near(species[0].particle[0].velocity.x,.95,"analytic backward x half-kick");
        near(species[0].particle[0].velocity.y,-1.975,"analytic backward y half-kick");
        check(species[0].particle[0].position.x==.25 && species[0].particle[0].vz==3. && state.rng()==rng,"half-step does not drift or draw RNG");
        rejected=false;try {state.initializeBackwardHalfStep(mesh,caches,field,.1,solid);}catch(const std::logic_error&){rejected=true;}
        check(rejected,"duplicate half-kick rejected");same(species[0].particle[0],expected[0].particle[0]);
        State restored(3,77);restored.initialize(config,norm.scales(),mesh,true);restored.species()=initial;
        rejected=false;try {restored.initializeBackwardHalfStep(mesh,caches,field,.1,solid);}catch(const std::logic_error&){rejected=true;}
        check(rejected,"restored particles cannot receive startup half-kick");same(restored.species()[0].particle[0],initial[0].particle[0]);
        if(variant==0) {
            State advanced(3,77);advanced.initialize(config,norm.scales(),mesh,false);advanced.species()=initial;
            advanced.configurePush(options);
            advanced.pushAll(mesh,caches,field,0.,nullptr);rejected=false;
            try {advanced.initializeBackwardHalfStep(mesh,caches,field,.1,nullptr);}catch(const std::logic_error&){rejected=true;}
            check(rejected,"starting a permanent push closes startup initialization");
        }
    }
    using Options=State::RelabelOptions;
    check(Options::fromExplicitConfig("LOCALMAPTHENFINDCELL","FINDCELLLABEL").useHint
        && Options::fromExplicitConfig("localMap","findCell").useFallback,"legacy relabel aliases are case insensitive");
    check(!Options::fromExplicitConfig("other","none").useHint && !Options::fromExplicitConfig("other","none").useFallback,"unknown relabel choices retain disabled legacy gates");
    QuadtreeMesh mesh(input);State state(3,99);state.initialize(config,norm.scales(),mesh,false);
    auto& species=state.species();
    for(int i=0;i<5;++i) {
        Particle p;p.position={i==4?1.1:.1+.15*i,.2};p.velocity={3.,4.};p.vz=5.;p.rngState=300+i;p.trackId=400+i;
        p.inCellL=i==1?-1:mesh.findCellLabel(p.position);species[0].particle.push_back(p);
    }
    const auto coarse=species;const auto rng=state.rng();
    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
    check(mesh.cellIndexNum()==64,"relabel actual refinement");
    auto compareRelabel=[&](const std::vector<ParticleType>& initial,bool hint,bool fallback) {
        species=initial;auto expected=initial;species[0].number=-8;
        for(auto& sp:expected) {
            for(auto& p:sp.particle) {
                int cell=p.inCellL;if(hint)cell=mesh.findCellLabelFromHint(cell,p.position);
                if(cell<0 && fallback)cell=mesh.findCellLabel(p.position);
                p.inCellL=cell;p.inBackgroundCellL=cell>=0?mesh.findBackgroundCellLabel(p.position):-1;
            }
            sp.number=static_cast<int>(sp.particle.size());
        }
        state.relabelAfterMeshChange(mesh,{hint,fallback});
        for(std::size_t si=0;si<species.size();++si) {
            check(species[si].number==expected[si].number && species[si].particle.size()==expected[si].particle.size(),"relabel updates counts without culling");
            for(std::size_t pi=0;pi<species[si].particle.size();++pi)same(species[si].particle[pi],expected[si].particle[pi]);
        }
        check(state.rng()==rng,"relabel draws no RNG");
        if(hint || (fallback && initial[0].particle[0].inCellL<0))
            check(species[0].particle[0].inCellL==mesh.findCellLabel(species[0].particle[0].position),"relabel locates current leaf");
        else check(species[0].particle[0].inCellL==initial[0].particle[0].inCellL,"disabled hint retains old nonnegative label");
    };
    for(bool hint:{false,true})for(bool fallback:{false,true})compareRelabel(coarse,hint,fallback);
    const auto refined=species;
    CellField<double> gradient(mesh,0.),target(mesh,0.);std::vector<VertField<double>> levelSets;
    std::vector<VertField<double>*> vertices;std::vector<HEdgeField<double>*> horizontal;
    std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> cells{&target};
    mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&target);
    check(mesh.cellIndexNum()==16,"relabel actual coarsening");compareRelabel(refined,true,true);
    State uninitialized(1,1);bool rejected=false;
    try {uninitialized.relabelAfterMeshChange(mesh,{true,true});}catch(const std::logic_error&){rejected=true;}
    check(rejected,"uninitialized relabel rejected");
    CellField<Vector2D> field(mesh,Vector2D(0.,0.));rejected=false;
    try {uninitialized.initializeBackwardHalfStep(mesh,{},field,.1,nullptr);}catch(const std::logic_error&){rejected=true;}
    check(rejected,"uninitialized half-kick rejected");
}

void speciesPush()
{
    using State=PicSpeciesState;using Pusher=ExplicitPIC::Pusher::ParticlePusherLeapfrog;
    static_assert(std::is_same_v<decltype(std::declval<const State&>().pushResult(0)),const State::PushResult&>);
    auto equalParticle=[](const Particle& a,const Particle& b){
        check(a.position.x==b.position.x && a.position.y==b.position.y && a.velocity.x==b.velocity.x
            && a.velocity.y==b.velocity.y && a.vz==b.vz && a.w==b.w && a.inCellL==b.inCellL
            && a.inBackgroundCellL==b.inBackgroundCellL && a.discreteVertL==b.discreteVertL
            && a.emitChance==b.emitChance && a.rngState==b.rngState && a.trackId==b.trackId
            && a.levelSetValues==b.levelSetValues,"pushed particle state/order differs from original kernel");
        for(int i=0;i<4;++i)check(a.weight[i]==b.weight[i],"particle interpolation history changed");};
    for(bool allSpecies:{false,true})for(int variant=0;variant<4;++variant) {
        DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);CellField<Vector2D> field(mesh,Vector2D(.5,-.25));
        DictionaryGeometry dict;dict.points={{.5,0.},{.5,1.}};Line wall;
        wall.points[0]={.5,0.};wall.points[1]={.5,1.};wall.setType(LineType::straight);wall.levelSetIndex=0;
        dict.lines={wall};dict.levelSetName={"wall"};
        dict.innerPoints={{.75,.5}};dict.linesOnInterface={{{0,false}}};dict.moving={false};
        ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh);geometry.lines[0].levelSetIndex=0;
        for(int vi=0;vi<mesh.vertIndexNum();++vi)geometry.levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.5;
        geometry.refreshLevelSetDerivatives(mesh);
        auto* geometryPtr=variant>=2?&geometry:nullptr;
        std::vector<SpeciesConfig> config(3);for(int i=0;i<3;++i){config[i].name="push_"+std::to_string(i);config[i].m=1.;config[i].q=i==1?0.:1.;config[i].weight=2.;}
        State owner(3,12345);const auto norm=NormalizationContext::MakeLegacy();owner.initialize(config,norm.scales(),mesh,true);
        auto& species=owner.species();species[0].set_kappa_s(0.);species[1].set_kappa_s(2.);
        const std::vector<Vector2D> positions{{.25,.25},{.02,.2},{.98,.2},{.2,.02},{.2,.98},{.6,.5},{.45,.5},{.55,.7}};
        const std::vector<Vector2D> velocities{{.1,.2},{-1.,0.},{1.,0.},{0.,-1.},{0.,1.},{0.,0.},{1.,0.},{-1.,0.}};
        for(std::size_t i=0;i<positions.size();++i){Particle p;p.position=positions[i];p.velocity=velocities[i];p.vz=3.;p.w=i+1.;
            p.rngState=123+i;p.trackId=20+i;p.inCellL=mesh.findCellLabel(p.position);p.inBackgroundCellL=mesh.findBackgroundCellLabel(p.position);species[0].particle.push_back(p);}
        Particle analytic;analytic.position={.1,.1};analytic.rngState=222;analytic.trackId=333;
        analytic.inCellL=mesh.findCellLabel(analytic.position);analytic.inBackgroundCellL=mesh.findBackgroundCellLabel(analytic.position);
        species[1].particle.push_back(analytic);auto reference=species;const auto rng=owner.rng();
        std::vector<CouplingWeightsCache> caches(3);
        auto build=[&]{for(std::size_t si=0;si<species.size();++si){std::vector<int> offsets,indices;std::vector<double> weights;
            for(const auto& p:species[si].particle){offsets.push_back(static_cast<int>(indices.size()));indices.push_back(mesh.cells(mesh.findCellLabel(p.position)).index());weights.push_back(1.);}
            offsets.push_back(static_cast<int>(indices.size()));caches[si].assignCSR(species[si].particle.size(),offsets,indices,weights);}};build();
        Pusher::Options options;options.periodicX=options.periodicY=variant==1;
        options.absorbAtEmbedded=variant>=2;options.absorbEmbeddedCrossing=variant==3;
        options.embeddedInsideIsNegative=false;options.embeddedLevelSetIds={0};
        owner.configurePush(options);
        const bool recordEmbedded=allSpecies || variant==3;
        auto step=[&]{
            std::vector<State::PushResult> expected(3);
            for(std::size_t si=0;si<reference.size();++si){if(reference[si].particle.empty())continue;std::vector<Vector2D> gathered;
                ExplicitPIC::Field::MatchedGather::gatherVector2D(caches[si],field,gathered,false);
                Pusher::advance(mesh,reference[si],gathered,.1,options,geometryPtr,&expected[si].domain,
                    recordEmbedded?&expected[si].embedded:nullptr);}
            if(allSpecies)owner.pushAll(mesh,caches,field,.1,geometryPtr);
            else {owner.beginPush();for(std::size_t si=0;si<species.size();++si)
                if(!species[si].particle.empty())owner.pushSpecies(si,mesh,caches[si],field,.1,geometryPtr,recordEmbedded);}
            for(std::size_t si=0;si<species.size();++si){
                check(species[si].number==reference[si].number && species[si].particle.size()==reference[si].particle.size(),"push survivor count");
                for(std::size_t p=0;p<species[si].particle.size();++p)equalParticle(species[si].particle[p],reference[si].particle[p]);
                const auto& a=owner.pushResult(si);const auto& b=expected[si];
                check(a.domain.size()==b.domain.size() && a.embedded.size()==b.embedded.size(),"owned impact sizes match original optional capture");
                for(std::size_t j=0;j<a.domain.size();++j){equalParticle(a.domain.impacts()[j].particle,b.domain.impacts()[j].particle);
                    check(a.domain.impacts()[j].side==b.domain.impacts()[j].side,"domain side and merge order");}
                for(std::size_t j=0;j<a.embedded.size();++j){const auto& x=a.embedded.impacts()[j];const auto& y=b.embedded.impacts()[j];equalParticle(x.particle,y.particle);
                    check(x.pos.x==y.pos.x && x.pos.y==y.pos.y && x.cellLabel==y.cellLabel && x.crossing==y.crossing
                        && x.boundaryId==y.boundaryId && x.levelSetIndex==y.levelSetIndex && x.lineIndex==y.lineIndex
                        && x.normal.x==y.normal.x && x.normal.y==y.normal.y,"embedded impact snapshot and classification order");}
            }
            check(owner.rng()==rng,"push consumes no master RNG");
        };
        step();near(species[1].particle[0].velocity.x,.1,"analytic gather kick");near(species[1].particle[0].velocity.y,-.05,"analytic y kick");
        near(species[1].particle[0].position.x,.11,"analytic drift uses kicked velocity");near(species[1].particle[0].position.y,.095,"analytic y drift");
        if(variant==1)check(owner.pushResult(0).domain.empty() && species[0].particle.size()==8,"periodic wrap precedes absorption");
        else check(owner.pushResult(0).domain.size()==4,"all four domain sides captured");
        if(variant>=2 && recordEmbedded)check(!owner.pushResult(0).embedded.empty(),"fixture exercises actual embedded absorption");
        if(variant==3){const bool crossed=std::any_of(owner.pushResult(0).embedded.impacts().begin(),owner.pushResult(0).embedded.impacts().end(),
            [](const auto& impact){return impact.crossing;});
            if(!crossed){std::cerr<<"push fixture all="<<allSpecies<<"\n";
                for(const auto& impact:owner.pushResult(0).embedded.impacts())std::cerr<<"impact id="<<impact.particle.trackId<<" x="<<impact.pos.x<<" crossing="<<impact.crossing<<'\n';
                for(const auto& p:species[0].particle)std::cerr<<"kept id="<<p.trackId<<" x="<<p.position.x<<'\n';}
            check(crossed,"crossing flag reaches owned report");}
        for(const auto& impact:owner.pushResult(0).embedded.impacts())
            check(impact.boundaryId==0 && impact.levelSetIndex==-1 && impact.lineIndex==-1,
                "crossing constructor must not write the flag into boundary identity metadata");
        const auto impacts=owner.pushResult(0).domain.size(),count=species[0].particle.size();bool rejected=false;
        try{owner.pushAll(mesh,{},field,.1,geometryPtr);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected && species[0].particle.size()==count && owner.pushResult(0).domain.size()==impacts,"short weights rejected before clearing records or pushing");
        if(variant==0){DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
            check(mesh.cellIndexNum()==64,"push real AMR fixture");field.initialize(mesh,Vector2D(.5,-.25));build();step();}
        species[0].particle.clear();species[0].number=0;reference[0].particle.clear();reference[0].number=0;
        build();step();check(owner.pushResult(0).domain.empty() && owner.pushResult(0).embedded.empty(),"empty species does not retain previous step impacts");
        owner.beginPush();for(std::size_t si=0;si<3;++si)check(owner.pushResult(si).domain.empty() && owner.pushResult(si).embedded.empty(),"new push stage clears all results");
    }
    {
        DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=32;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);CellField<Vector2D> field(mesh,Vector2D(0.,0.));
        SpeciesConfig config;config.name="parallel";config.m=1.;config.q=0.;
        State owner(1,7);owner.initialize({config},NormalizationContext::MakeLegacy().scales(),mesh,true);
        auto& species=owner.species()[0];std::vector<int> offsets,indices;std::vector<double> weights;
        for(int i=0;i<20016;++i){Particle p;p.position={i%2?.5:.05,.5};p.velocity={-1.,0.};p.trackId=i;p.rngState=900+i;
            p.inCellL=mesh.findCellLabel(p.position);p.inBackgroundCellL=mesh.findBackgroundCellLabel(p.position);species.particle.push_back(p);
            offsets.push_back(i);indices.push_back(mesh.cells(p.inCellL).index());weights.push_back(1.);}
        offsets.push_back(20016);std::vector<CouplingWeightsCache> cache(1);cache[0].assignCSR(20016,offsets,indices,weights);
        auto reference=species;std::vector<Vector2D> gathered;ExplicitPIC::Field::MatchedGather::gatherVector2D(cache[0],field,gathered,false);
        State::PushResult expected;Pusher::Options options;options.absorbAtEmbedded=false;
        Pusher::advance<QuadtreeMesh>(mesh,reference,gathered,.1,options,nullptr,&expected.domain,&expected.embedded);
        owner.configurePush(options);
        owner.pushAll(mesh,cache,field,.1,nullptr);
        check(species.particle.size()==10008 && owner.pushResult(0).domain.size()==10008,"parallel threshold fixture survivors and domain impacts");
        for(std::size_t i=0;i<10008;++i){equalParticle(species.particle[i],reference.particle[i]);
            equalParticle(owner.pushResult(0).domain.impacts()[i].particle,expected.domain.impacts()[i].particle);
            check(species.particle[i].trackId==2*i+1 && owner.pushResult(0).domain.impacts()[i].particle.trackId==2*i,"parallel compaction and impact merge order");}
    }
    State uninitialized(1,9);bool rejected=false;try{uninitialized.beginPush();}catch(const std::logic_error&){rejected=true;}
    check(rejected,"push requires initialized species");
}
void implicitFieldRebuild(const std::filesystem::path& directory)
{
    using Solver=PicPoissonFieldSolve;using Policy=PicCurrentWeights::Policy;
    const auto bytes=[](const std::filesystem::path& path) {
        std::ifstream in(path,std::ios::binary);check(bool(in),"field rebuild fixture file");
        return std::string(std::istreambuf_iterator<char>(in),{});
    };
    for(bool periodic:{false,true}) {
        const auto dir=directory/(periodic?"rebuild_periodic":"rebuild_domain");
        const auto snapshot=[&](const ImplicitPicFieldHistory& history,const std::string& name) {
            const auto path=dir/name;std::filesystem::create_directories(path);history.writeCheckpoint(path);
            return bytes(path/PicCheckpointStore::fieldName);
        };
        DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input),otherMesh(input);FaceSegmentGraph graph;graph.build(mesh);
        PicFieldState actual(mesh),reference(mesh),foreign(mesh);
        ImplicitPicFieldHistory history(actual),historyRef(reference),historyForeign(foreign);
        Solver::Options options;options.cacheEnabled=false;options.cacheLogOnMiss=false;
        Solver solver(actual,options),solverRef(reference,options),solverForeign(foreign,options);
        PicCurrentWeights weights(mesh),weightsRef(mesh),wrongWeights(otherMesh);
        DepositionConfig config;Policy coupling(mesh,config,Policy::Mode::Implicit,periodic,false);
        Policy wrongPolicy(otherMesh,config,Policy::Mode::Implicit,false,false);
        Policy explicitPolicy(mesh,config,Policy::Mode::Explicit,periodic,false);
        PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Implicit);
        const auto norm=NormalizationContext::MakeLegacy();std::vector<ParticleType> species;
        species.emplace_back("charge",1.,.02*norm.scales().rho0_Cpm3);species.emplace_back("empty",1.,0.);
        for(int i=0;i<3;++i) {Particle p;p.position={.2+.2*i,.3};p.w=i;p.rngState=700+i;p.trackId=800+i;
            p.inCellL=mesh.findCellLabel(p.position);species[0].particle.push_back(p);}
        const std::vector<double> macroWeight{3.};
        const Solver::ImplicitRebuildInput particles{species,macroWeight,norm,coupling,embedded,false};
        PicDomainBC bc;bc.right.value=1.;bc.top.type=bc.bottom.type=PicBCType::Neumann;
        if(periodic) {bc.left.type=bc.right.type=PicBCType::Periodic;bc.top.type=bc.bottom.type=PicBCType::Dirichlet;bc.top.value=1.;}
        const Solver::Input fieldInput{graph,bc};
        const auto original=[&](bool seed) {
            weightsRef.rebuild(species,coupling,embedded,false);
            reference.depositFullCellCharge(species,macroWeight,weightsRef.caches(),norm);
            reference.copyChargeDensityToRhs();solverRef.solveForRebuild(fieldInput);
            if(seed)historyRef.seedFromCurrent();
        };
        for(int step=0;step<3;++step) {
            if(step==1) {species[0].particle[0].position.y+=.1;species[0].particle[0].inCellL=mesh.findCellLabel(species[0].particle[0].position);bc.right.value=2.;}
            if(step==2 && !periodic) {
                const auto old=weights.caches()[0].weights();DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;
                mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()==64,"implicit rebuild actual AMR");graph.build(mesh);
                bool rejected=false;try {solver.rebuildImplicit(fieldInput,particles,weights,&history);}catch(const std::logic_error&){rejected=true;}
                check(rejected && weights.caches()[0].weights()==old,"stale fields rejected before rebuilding weights");
                CellField<double> mapped(mesh,0.);for(int ci=0;ci<mesh.cellIndexNum();++ci)mapped.value[ci]=2.*mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x;
                actual.resetAfterMeshChange(&mapped);reference.resetAfterMeshChange(&mapped);rejected=false;
                try {solver.reseedImplicitFromPotential(fieldInput,history);}catch(const std::logic_error&){rejected=true;}
                check(rejected && actual.electricField().value[0].x==0.,"stale history rejected before remap E evaluation");
                history.resetAfterMeshChange();historyRef.resetAfterMeshChange();coupling.refreshAfterMeshChange();
                solver.reseedImplicitFromPotential(fieldInput,history);
                solverRef.evaluateCurrent(fieldInput,Solver::BoundaryPolicy::LegacyDefault);historyRef.seedFromCurrent();
                check(snapshot(history,"fast_actual")==snapshot(historyRef,"fast_reference"),"fast remap keeps phi/charge and seeds exact E history");
                for(int ci=0;ci<64;++ci) {near(actual.electricField().value[ci].x,-2.,"analytic remapped E");near(actual.chargeAmount().value[ci],0.,"fast remap does not deposit charge");}
                for(auto& p:species[0].particle)p.inCellL=mesh.findCellLabel(p.position);
            }
            solver.rebuildImplicit(fieldInput,particles,weights,step==0?nullptr:&history);original(step!=0);
            check(snapshot(history,"actual")==snapshot(historyRef,"reference"),"particle rebuild preserves all seven fields and history stages");
            for(std::size_t si=0;si<species.size();++si) {
                const auto& a=weights.caches()[si];const auto& b=weightsRef.caches()[si];
                check(a.offsets()==b.offsets() && a.cellIndices()==b.cellIndices() && a.weights()==b.weights(),"implicit rebuild exact current-position CSR");
            }
            if(step==0)near(history.nextElectricField().value[0].x,0.,"initial rebuild does not seed later-created history");
            for(std::size_t pi=0;pi<species[0].particle.size();++pi)
                check(species[0].particle[pi].rngState==700+pi && species[0].particle[pi].trackId==800+pi,"field rebuild preserves particle RNG/identity");
        }
        const auto saved=weights.caches()[0].weights();const auto before=snapshot(history,"before_rejection");
        auto reject=[&](PicCurrentWeights& cache,const Policy& policy,ImplicitPicFieldHistory* h) {
            bool caught=false;try {solver.rebuildImplicit(fieldInput,{species,macroWeight,norm,policy,embedded,false},cache,h);}
            catch(const std::invalid_argument&){caught=true;}
            check(caught && weights.caches()[0].weights()==saved && snapshot(history,"after_rejection")==before,"foreign rebuild owner rejected before mutation");
        };
        reject(wrongWeights,coupling,&history);reject(weights,wrongPolicy,&history);
        reject(weights,explicitPolicy,&history);reject(weights,coupling,&historyForeign);
        // Preserve distinct current/next/half layers; a valid restore must not
        // replace this time history with a freshly solved, uniformly seeded E.
        PicDomainBC nextBC=bc;if(periodic)nextBC.top.value=3.;else nextBC.right.value=3.;
        Solver::Report report;solver.solveExplicit({graph,nextBC},3,report);
        solver.evaluateNext({graph,nextBC},history);history.formHalfStep();
        check(history.nextElectricField().value[0].x!=actual.electricField().value[0].x
            || history.nextElectricField().value[0].y!=actual.electricField().value[0].y,"restore fixture has distinct E time layers");
        const auto encoded=snapshot(history,"valid");weights.clear();
        check(solver.restoreImplicit(dir/"valid",fieldInput,particles,weights,history) && weights.caches().empty(),"valid field restore skips all particle rebuild work");
        check(snapshot(history,"valid_roundtrip")==encoded,"IPF1 field layers survive exact restoration");
        check(solverForeign.restoreImplicit(dir/"valid",fieldInput,particles,weights,historyForeign),"valid reader may initialize stale field/history storage");
        check(snapshot(historyForeign,"valid_foreign")==encoded,"stale owner restored to current checkpoint dimensions");
        bool rejected=false;
        try {solver.restoreImplicit(dir/"valid",fieldInput,particles,weights,historyForeign);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected,"restore rejects foreign field history before IO");
        for(const auto& fixture:{std::string("missing"),std::string("bad_header"),std::string("truncated")}) {
            const auto path=dir/fixture;std::filesystem::create_directories(path);
            if(fixture!="missing") {std::ofstream stream(path/PicCheckpointStore::fieldName,std::ios::binary);
                const auto content=fixture=="bad_header"?std::string("invalid IPF1"):encoded.substr(0,encoded.size()-1);
                stream.write(content.data(),static_cast<std::streamsize>(content.size()));}
            BoundaryLogs logs;check(!solver.restoreImplicit(path,fieldInput,particles,weights,history),"rejected field checkpoint invokes rebuild");
            check(!historyRef.readCheckpoint(path),"original reader rejects same fixture");original(true);
            check(snapshot(history,"fallback_actual")==snapshot(historyRef,"fallback_reference"),"fallback replaces partial read and seeds all E layers exactly");
        }
        const auto electric=actual.electricField().value[0];const auto next=history.nextElectricField().value[0];
        species[0].particle[0].w=7.;FaceSegmentGraph empty;empty.build(mesh,[](const QuadtreeMesh&,const Cell&){return false;});
        {BoundaryLogs logs;bool failed=false;try {solver.rebuildImplicit({empty,bc},particles,weights,&history);}
            catch(const std::runtime_error&){failed=true;}check(failed,"rebuild propagates numerical solver failure");}
        check(actual.electricField().value[0].x==electric.x && actual.electricField().value[0].y==electric.y
            && history.nextElectricField().value[0].x==next.x && history.nextElectricField().value[0].y==next.y,"failed solve does not evaluate or seed E history");
    }
}

void implicitFieldAdvance()
{
    const auto norm=NormalizationContext::MakeLegacy();const double q=norm.scales().rho0_Cpm3;
    const double dt=std::sqrt(PhysicalConstants::epsilon0_Fpm);
    for(int closure=0;closure<4;++closure) {
        DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
        PicFieldState actual(mesh),reference(mesh),foreign(mesh);
        PicPoissonFieldSolve::Options options;options.cacheEnabled=false;
        if(closure==3){options.useIC=false;options.maxIter=1;options.relTol=1e-14;}
        PicPoissonFieldSolve eval(actual,options),evalRef(reference,options),evalForeign(foreign,options);
        ImplicitPicFieldHistory history(actual),historyRef(reference),historyForeign(foreign);
        ImplicitPicPolarization pol(mesh,2,graph.segments().size()),polRef(mesh,2,graph.segments().size());
        PicConfig::PolarizationPoissonConfig config;config.mode=closure==3?"auto":"direct";
        ImplicitPicFieldSolve solver(actual,options,config,false,false),solverRef(reference,options,config,false,false);
        PicDomainBC initial;initial.right.value=1.;initial.bottom.type=initial.top.type=PicBCType::Neumann;
        auto seed=[&]{eval.solveForRebuild({graph,initial});evalRef.solveForRebuild({graph,initial});
            history.seedFromCurrent();historyRef.seedFromCurrent();};seed();
        PicDomainBC bc=initial;const bool closed=closure==1 || closure==2,periodic=closure==2;
        if(closed)for(auto* side:{&bc.left,&bc.right,&bc.bottom,&bc.top}){side->type=periodic?PicBCType::Periodic:PicBCType::Neumann;side->value=0.;}
        PicEmbeddedBC embedded;embedded.enabled=true;embedded.type=PicBCType::Dirichlet;embedded.g=7.;
        // Even a configured embedded Dirichlet does not change the legacy
        // domain-only neutrality decision; cut-geometry cases run separately.
        const ImplicitPicFieldSolve::Input fieldInput{graph,bc,nullptr,&embedded};
        std::vector<ParticleType> species;species.emplace_back("positive",q*q,q);species.emplace_back("empty",1.,0.);
        std::vector<double> defaults{2.,3.};std::vector<CouplingWeightsCache> caches(2);
        auto fill=[&]{const int n=mesh.cellIndexNum();species[0].particle.resize(n);
            std::vector<int> offsets,indices;std::vector<double> weights;
            for(int ci=0;ci<n;++ci){auto& p=species[0].particle[ci];p.w=1.+ci%5;p.rngState=123+ci;
                offsets.push_back(ci);indices.push_back(ci);weights.push_back(1.);}
            offsets.push_back(n);caches[0].assignCSR(n,offsets,indices,weights);caches[1].assignCSR(0,{0},{},{});};fill();
        BoundaryLogs logs;
        for(int step=0;step<3;++step) {
            if(step==2) {
                DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
                graph.build(mesh);check(mesh.cellIndexNum()==64,"field advance real AMR fixture");
                actual.resetAfterMeshChange();reference.resetAfterMeshChange();
                pol.resetDepositedFields();pol.resetCoefficients(graph.segments().size());
                polRef.resetDepositedFields();polRef.resetCoefficients(graph.segments().size());fill();
                PicStepParticleBalance untouched(species,defaults);PicImplicitStepTiming timing;bool rejected=false;
                try{solver.advance(fieldInput,{species,defaults,caches,norm,dt,step,1,periodic,periodic,false},pol,eval,history,untouched,timing);}
                catch(const std::logic_error&){rejected=true;}
                check(rejected && pol.totalCharge()==0. && untouched.current().stageMacroC==0.
                    && timing.depositOrbit_qp_size_total==0,"stale history rejected before deposition after AMR");
                history.resetAfterMeshChange();historyRef.resetAfterMeshChange();
                eval.resetAfterMeshChange();evalRef.resetAfterMeshChange();solver.resetAfterMeshChange(step);solverRef.resetAfterMeshChange(step);seed();
            }
            species[0].particle[0].w=step+2.;if(!closed)bc.right.value=step+1.;
            const auto current=actual.electricField();
            PicStepParticleBalance balance(species,defaults),balanceRef(species,defaults);
            PicImplicitStepTiming timing;timing.depositOrbit_qp_size_total=7;timing.depositOrbit_nnz_total=11;
            std::vector<double> normal{999.},normalRef;const bool wantNormal=step!=1;
            auto result=solver.advance(fieldInput,{species,defaults,caches,norm,dt,step,1,periodic,periodic,true,wantNormal?&normal:nullptr},
                pol,eval,history,balance,timing);
            // Frozen former driver order, independent of the new composition.
            const auto deposited=polRef.deposit(species,defaults,caches,1,balanceRef);
            polRef.assembleCells(species,norm,dt);if(closed)polRef.removeMeanCharge(norm);
            balanceRef.finishStage(polRef.totalCharge());polRef.assembleFaces(graph);
            const auto expected=solverRef.solve(fieldInput,polRef,step+1);
            solverRef.enforcePeriodicSymmetry(fieldInput,periodic,periodic);
            evalRef.evaluateNext(fieldInput,historyRef,wantNormal?&normalRef:nullptr);historyRef.formHalfStep();
            check(result.report.success && result.report.N==expected.report.N
                && result.report.cgAttempted==expected.report.cgAttempted
                && result.report.cgFallbackToDirect==expected.report.cgFallbackToDirect,"composed stage retains solve policy and fallback");
            if(closure==3)check(result.report.cgAttempted==(step!=1),"auto fallback is suppressed until AMR retry");
            for(int ci=0;ci<mesh.cellIndexNum();++ci) {
                check(actual.potential().value[ci]==reference.potential().value[ci]
                    && pol.chargeAmount().value[ci]==polRef.chargeAmount().value[ci]
                    && pol.chargeDensity().value[ci]==polRef.chargeDensity().value[ci],"field stage potential and deposits match original order exactly");
                const auto a=history.nextElectricField().value[ci],b=historyRef.nextElectricField().value[ci];
                const auto h=history.halfElectricField().value[ci],hr=historyRef.halfElectricField().value[ci];
                check(a.x==b.x && a.y==b.y && h.x==hr.x && h.y==hr.y,"next/half layers match original order exactly");
                check(actual.electricField().value[ci].x==current.value[ci].x && actual.electricField().value[ci].y==current.value[ci].y,
                    "field substep does not accept the next electric field");
            }
            check(pol.cellCoefficients()==polRef.cellCoefficients() && pol.faceCoefficients()==polRef.faceCoefficients(),"same polarization coefficients");
            check(wantNormal?normal==normalRef:normal==std::vector<double>{999.},"normal fields generated only for requested destination");
            check(balance.current().stageMacroC==balanceRef.current().stageMacroC && balance.current().stageGridC==balanceRef.current().stageGridC
                && balance.current().stageRelative==balanceRef.current().stageRelative,"stage balance observes corrected deposition before solve");
            near(balance.current().stageRelative,closed?-1.:0.,"domain-only neutrality correction");
            check(timing.depositOrbit_qp_size_total==7+deposited.particles && timing.depositOrbit_nnz_total==11+deposited.entries
                && timing.polPoisson_total_s==result.solveSeconds && timing.pol_cg_iterations==result.report.cgIterations,
                "stage timings/counts retain accumulation and solver report semantics");
            check(species[0].particle[0].w==step+2. && species[0].particle[0].rngState==123 && caches[0].numParticles()==mesh.cellIndexNum(),
                "field stage does not mutate particle weights, RNG or input cache");
            history.acceptNext();historyRef.acceptNext();
        }
        check((logs.output.str().find("neutrality correction applied")!=std::string::npos)==closed,"original first-step debug notice retained");
        PicStepParticleBalance untouched(species,defaults);PicImplicitStepTiming timing;
        auto reject=[&](ImplicitPicPolarization& p,PicPoissonFieldSolve& e,ImplicitPicFieldHistory& h){bool caught=false;
            try{solver.advance(fieldInput,{species,defaults,caches,norm,dt,3,1,periodic,periodic,false},p,e,h,untouched,timing);}
            catch(const std::invalid_argument& ex){caught=std::string(ex.what()).find("one field state and mesh")!=std::string::npos;}
            check(caught && untouched.current().stageMacroC==0. && timing.depositOrbit_qp_size_total==0,"mixed ownership rejected before work");};
        const double before=pol.totalCharge();reject(pol,evalForeign,history);reject(pol,eval,historyForeign);
        DictionaryMesh otherInput=input;otherInput.Nx=otherInput.Ny=8;QuadtreeMesh otherMesh(otherInput);
        ImplicitPicPolarization otherPol(otherMesh,2,graph.segments().size());reject(otherPol,eval,history);
        check(pol.totalCharge()==before && otherPol.totalCharge()==0.,"ownership rejection preserves polarization state");
        FaceSegmentGraph empty;empty.build(mesh,[](const QuadtreeMesh&,const Cell&){return false;});
        const auto accepted=actual.electricField(),half=history.halfElectricField();std::vector<double> normals{999.};bool failed=false;
        try{solver.advance({empty,bc},{species,defaults,caches,norm,dt,3,1,periodic,periodic,false,&normals},
                pol,eval,history,untouched,timing);}
        catch(const std::runtime_error& ex){failed=std::string(ex.what()).find("[ImplicitPIC] polarization Poisson failed:")==0;}
        check(failed && untouched.current().stageMacroC!=0. && timing.depositOrbit_qp_size_total>0 && normals==std::vector<double>{999.},
            "solve failure propagates after deposition/accounting and before E evaluation, without claiming rollback");
        for(int ci=0;ci<mesh.cellIndexNum();++ci)check(actual.electricField().value[ci].x==accepted.value[ci].x
            && actual.electricField().value[ci].y==accepted.value[ci].y && history.halfElectricField().value[ci].x==half.value[ci].x
            && history.halfElectricField().value[ci].y==half.value[ci].y,"failed field substep does not accept or form a new electric layer");
    }
}
void implicitFieldSolve()
{
    using Stage=ImplicitPicFieldSolve;
    static_assert(!std::is_copy_constructible_v<Stage> && !std::is_move_constructible_v<Stage>);
    static_assert(std::is_same_v<decltype(std::declval<Stage&>().options()),const Stage::Options&>);
    ScopedEnvironment modeEnv("LAMP_IMPLICIT_POL_MODE",nullptr),guessEnv("LAMP_IMPLICIT_POL_CG_TEMPORAL_GUESS",nullptr),
        omegaEnv("LAMP_IMPLICIT_POL_CG_TEMPORAL_OMEGA",nullptr),patternEnv("LAMP_IMPLICIT_POL_CG_REUSE_PATTERN",nullptr);
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);PicFieldState fields(mesh);
    ImplicitPicPolarization polarization(mesh,1,graph.segments().size());
    PicConfig::PolarizationPoissonConfig config;config.mode="DiReCt";config.freezeDirectMatrixUntilPatternChange=true;
    PicPoissonFieldSolve::Options options;Stage direct(fields,options,config,false,false);
    check(direct.mode()=="direct" && !direct.options().useCG && !direct.options().useIC && direct.options().freezeDirectMatrixUntilPatternChange,"direct mode and fixed-geometry freeze policy");
    Stage moving(fields,options,config,true,false);check(!moving.options().freezeDirectMatrixUntilPatternChange,"scheduled geometry prevents matrix freezing");
    PicDomainBC bc;bc.right.value=1.;bc.bottom.type=bc.top.type=PicBCType::Neumann;
    auto result=direct.solve({graph,bc},polarization,1);
    check(result.report.success && result.report.usedSuperLU && !result.report.cgAttempted && result.solveSeconds>=0.,"direct result and isolated solve timing");
    for(int ci=0;ci<16;++ci)near(fields.potential().value[ci],mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x,"analytic implicit potential");
    near(fields.electricField().value[0].x,0.,"polarization solve leaves electric stage untouched");
    auto repeated=direct.solve({graph,bc},polarization,2);check(repeated.report.directMatrixFrozen,"unchanged direct matrix reused");
    // Projection has its own stage: clean opposite cells are averaged once,
    // interior cells and the electric layer remain untouched.
    direct.enforcePeriodicSymmetry({graph,bc},true,false);
    for(int ci=0;ci<16;++ci) {
        const double x=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x;
        near(fields.potential().value[ci],(x<.25 || x>.75)?.5:x,"clean periodic pair projection");
    }
    direct.solve({graph,bc},polarization,3);
    // Nonuniform deposited charge and susceptibility exercise the actual
    // coefficient path, independently compared with the legacy solver port.
    const auto norm=NormalizationContext::MakeLegacy();const double q=norm.scales().rho0_Cpm3;
    std::vector<ParticleType> species;species.emplace_back("positive",q*q,q);species[0].particle.resize(16);
    std::vector<int> offsets,indices;std::vector<double> weights;
    for(int ci=0;ci<16;++ci){species[0].particle[ci].w=1.+ci%5;offsets.push_back(ci);indices.push_back(ci);weights.push_back(1.);}
    offsets.push_back(16);std::vector<CouplingWeightsCache> caches(1);caches[0].assignCSR(16,offsets,indices,weights);
    PicStepParticleBalance balance(0.);polarization.deposit(species,{},caches,10000,balance);
    polarization.assembleCells(species,norm,std::sqrt(PhysicalConstants::epsilon0_Fpm));polarization.assembleFaces(graph);
    result=direct.solve({graph,bc},polarization,4);check(result.report.success && !result.report.directMatrixFrozen,"coefficient change refactors despite fixed sparsity");
    PolarizationPoissonFVSolverCellCentered legacy;CellField<double> expected(mesh,0.);Stage::Report report;
    check(legacy.solve(mesh,graph,bc,polarization.chargeDensity(),expected,&polarization.faceCoefficients(),&polarization.cellCoefficients(),&report,direct.options()),"legacy reference solve");
    for(int ci=0;ci<16;++ci)near(fields.potential().value[ci],expected.value[ci],"module consumes current deposited coefficients");
    // A one-iteration CG budget causes a real fallback. Only auto mode
    // suppresses the next iterative attempt, and only AMR permits a retry.
    PicFieldState autoFields(mesh);config.mode="auto";options.useIC=false;options.maxIter=1;options.relTol=1e-14;
    BoundaryLogs logs;Stage automatic(autoFields,options,config,false,true);
    result=automatic.solve({graph,bc},polarization,1);
    check(result.report.cgAttempted && result.report.cgFallbackToDirect && automatic.iterativeSuppressed(),"real failed iteration suppresses auto mode");
    result=automatic.solve({graph,bc},polarization,2);check(!result.report.cgAttempted,"suppressed next step uses direct solve");
    automatic.clearDirectCache();check(automatic.iterativeSuppressed(),"geometry invalidation preserves suppression");
    result=automatic.solve({graph,bc},polarization,3);check(!result.report.cgAttempted,"geometry rebuild does not retry iterative path");
    config.mode="cg";PicFieldState cgFields(mesh);Stage forced(cgFields,options,config,false,false);
    result=forced.solve({graph,bc},polarization,1);check(result.report.cgAttempted && !forced.iterativeSuppressed(),"explicit CG never latches auto suppression");
    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
    check(mesh.cellIndexNum()==64,"implicit field real AMR");graph.build(mesh);autoFields.resetAfterMeshChange();polarization.resetDepositedFields();polarization.resetCoefficients(graph.segments().size());
    automatic.resetAfterMeshChange(4);check(!automatic.iterativeSuppressed(),"AMR re-enables iterative path");
    result=automatic.solve({graph,bc},polarization,5);check(result.report.success && result.report.N==64 && result.report.cgAttempted,"post-AMR solve retries on current topology");
    check(logs.output.str().find("suppress iterative path after fallback at step=1")!=std::string::npos && logs.output.str().find("re-enable iterative path after AMR at step=4")!=std::string::npos,"retry events keep original step labels");
    FaceSegmentGraph empty;empty.build(mesh,[](const QuadtreeMesh&,const Cell&){return false;});bool rejected=false;
    try{automatic.solve({empty,bc},polarization,6);}catch(const std::runtime_error& error){rejected=std::string(error.what()).find("[ImplicitPIC] polarization Poisson failed:")==0;}
    check(rejected,"implicit failure prefix preserved");
    // IO captures the legacy overrides; physical construction must use only
    // that owned input, even if the process environment changes afterwards.
    ScopedEnvironment::set("LAMP_IMPLICIT_POL_MODE","IC");ScopedEnvironment::set("LAMP_IMPLICIT_POL_CG_TEMPORAL_GUESS","OFF");
    ScopedEnvironment::set("LAMP_IMPLICIT_POL_CG_TEMPORAL_OMEGA","bad");ScopedEnvironment::set("LAMP_IMPLICIT_POL_CG_REUSE_PATTERN","no");
    const auto captured=readImplicitPicOverrides();
    ScopedEnvironment::set("LAMP_IMPLICIT_POL_MODE","direct");
    Stage independent(autoFields,options,config,false,false);
    check(independent.mode()=="cg","programmatic solver input is independent of process overrides");
    Stage overridden(autoFields,options,config,false,false,captured.polarization);
    check(overridden.mode()=="ic" && overridden.options().useIC && !overridden.options().useCG && !overridden.options().useTemporalExtrapolatedGuess && !overridden.options().reuseIterativePattern,"environment policy overrides");
    near(overridden.options().temporalExtrapolationOmega,.8,"malformed omega retains fallback");
    ScopedEnvironment::set("LAMP_IMPLICIT_POL_MODE","invalid");ScopedEnvironment::set("LAMP_IMPLICIT_POL_CG_TEMPORAL_OMEGA","-2");options.useCG=options.useIC=false;
    Stage fallback(autoFields,options,config,false,false,readImplicitPicOverrides().polarization);
    check(fallback.mode()=="auto" && fallback.options().useIC && !fallback.options().useCG,"unknown mode and disabled inherited iterations select IC auto");
    near(fallback.options().temporalExtrapolationOmega,0.,"negative omega clamped");check(overridden.mode()=="ic","environment changes do not mutate prior service");
}
void explicitFieldAdvance()
{
    using Solver=PicPoissonFieldSolve;
    using Policy=PicCurrentWeights::Policy;
    for(bool cached:{false,true}) for(int variant=0;variant<3;++variant) {
        DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input),foreignMesh(input);FaceSegmentGraph graph;graph.build(mesh);
        PicFieldState actual(mesh),reference(mesh),foreign(mesh);
        ExplicitPicChargeDeposition deposit(actual),depositRef(reference),depositForeign(foreign);
        Solver::Options options;options.cacheEnabled=cached;options.cacheLogOnMiss=false;
        Solver solver(actual,options),solverRef(reference,options);
        PicCurrentWeights weights(mesh),weightsRef(mesh),wrongWeights(foreignMesh);
        DepositionConfig config;Policy coupling(mesh,config,Policy::Mode::Explicit,variant==1,false);
        Policy wrongPolicy(foreignMesh,config,Policy::Mode::Explicit,false,false);
        Policy implicitPolicy(mesh,config,Policy::Mode::Implicit,false,false);
        PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);
        DictionaryGeometry dict;std::unique_ptr<ImmersedGeometry<QuadtreeMesh>> geometry;
        if(variant==2) {
            dict.levelSetName={"solid"};dict.innerPoints={{.1,.5}};dict.moving={false};
            dict.linesOnInterface={{{0,false}}};dict.phaseName={"solid"};dict.levelSetOperation={{LevelSetOperation::non}};
            dict.levelSetIndexA={{{0,false}}};dict.levelSetIndexB={{}};dict.boundaryLines={{{0,false}}};
            Line line;line.points[0]={.625,0.};line.points[1]={.625,1.};dict.lines={line};
            DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.kind=GeometricBoundary::fixed;
            boundary.levelSetIndex=0;boundary.phase=0;boundary.lineIds={0};dict.boundaries={boundary};dict.boundaryById[1]=boundary;
            geometry=std::make_unique<ImmersedGeometry<QuadtreeMesh>>(dict,mesh);
            geometry->levelSets[0].initialize(mesh,0.);geometry->LS_dx(0).initialize(mesh,1.);geometry->LS_dy(0).initialize(mesh,0.);
            for(int vi=0;vi<mesh.vertIndexNum();++vi)
                geometry->levelSets[0].value[vi]=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x-.625;
            embedded.resolveSelection(dict,1);embedded.rebuildComposite(mesh,geometry.get());
            embedded.rebuildMetrics(mesh,geometry.get(),graph,{1});
            check(embedded.metricsReady(),"field-stage cut metrics ready");
        }
        const auto norm=NormalizationContext::MakeLegacy();
        std::vector<ParticleType> species;species.emplace_back("electron",1.,-.02*norm.scales().rho0_Cpm3);
        species.emplace_back("ion",2.,.03*norm.scales().rho0_Cpm3);species.emplace_back("empty",1.,0.);
        for(int si=0;si<2;++si) for(int pi=0;pi<3;++pi) {
            Particle p;p.position={.72+.08*pi,.2+.2*pi};p.w=pi==0?0.:double(pi);
            p.inCellL=mesh.findCellLabel(p.position);p.rngState=100+10*si+pi;p.trackId=200+10*si+pi;
            species[si].particle.push_back(p);
        }
        const std::vector<double> macroWeight{2.,3.};
        PicDomainBC domain;domain.bottom.type=domain.top.type=PicBCType::Neumann;
        domain.right.timeDirichlet.enabled=true;domain.right.timeDirichlet.type="sine";
        domain.right.timeDirichlet.dc=2.;domain.right.timeDirichlet.ac=1.;domain.right.timeDirichlet.freq_hz=1.;
        if(variant==1) {
            domain.left.type=domain.right.type=PicBCType::Periodic;
            domain.bottom.type=domain.top.type=PicBCType::Dirichlet;domain.top.value=1.;
        }
        PicEmbeddedBC embeddedBC;PicEmbeddedBCMap byId;
        if(variant==2) {embeddedBC.enabled=true;embeddedBC.type=PicBCType::Dirichlet;
            embeddedBC.timeDirichlet=domain.right.timeDirichlet;byId[1]=embeddedBC;byId[1].timeDirichlet.dc=.5;}
        PicFieldBoundaryState boundaries(domain,embeddedBC,byId),boundariesRef(domain,embeddedBC,byId);
        std::vector<double> normal{99.},normalRef{99.};Solver::Report report,reportRef;
        PicExplicitStepTiming timing;
        for(int step=0;step<3;++step) {
            if(step==1) for(auto& sp:species) for(auto& p:sp.particle) {
                p.position.y+=.03;p.inCellL=mesh.findCellLabel(p.position);
            }
            if(step==2 && variant==0) {
                DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
                check(mesh.cellIndexNum()==64,"explicit field stage real AMR");graph.build(mesh);
                const auto oldWeights=weights.caches()[0].weights();bool rejected=false;
                try {solver.advanceExplicit(graph,{species,macroWeight,norm,coupling,embedded,nullptr,.25,step,true,0},
                    weights,deposit,boundaries,report,&timing);}catch(const std::logic_error&){rejected=true;}
                check(rejected && weights.caches()[0].weights()==oldWeights,"stale AMR fields rejected before cache rebuild");
                actual.resetAfterMeshChange();reference.resetAfterMeshChange();deposit.resetAfterMeshChange();depositRef.resetAfterMeshChange();
                solver.resetAfterMeshChange();solverRef.resetAfterMeshChange();coupling.refreshAfterMeshChange();
                for(auto& sp:species) sp.particle.clear();
            }
            const double time=step*.125;
            const bool collect=step!=1;
            timing={-1.,-1.,-1.,-1.,-1.};
            const Solver::ExplicitStageInput stage{species,macroWeight,norm,coupling,embedded,geometry.get(),time,step,collect,0,
                step==1?nullptr:&normal};
            solver.advanceExplicit(graph,stage,weights,deposit,boundaries,report,&timing);
            // Original driver order, with separately owned fields and cache.
            weightsRef.rebuild(species,coupling,embedded);
            depositRef.deposit(species,macroWeight,weightsRef.caches(),norm,embedded,geometry.get(),collect,0);
            reference.copyChargeDensityToRhs();boundariesRef.update(time);
            const auto* cut=embedded.metricsReady()?&embedded.metrics():nullptr;
            const Solver::Input solveInput{graph,boundariesRef.domain(),cut,
                (cut && boundariesRef.embedded().enabled)?&boundariesRef.embedded():nullptr,
                (cut && !boundariesRef.embeddedById().empty())?&boundariesRef.embeddedById():nullptr};
            solverRef.solveExplicit(solveInput,step,reportRef);
            solverRef.evaluateCurrent(solveInput,Solver::BoundaryPolicy::Configured,step==1?nullptr:&normalRef);
            check(report.success && report.N==reportRef.N && report.nnz==reportRef.nnz
                && report.usedSuperLU==reportRef.usedSuperLU && report.message==reportRef.message,"explicit stage retains solver/cache report");
            check(normal==normalRef,"same-stage normal fields preserve exact values and disabled output");
            for(std::size_t si=0;si<species.size();++si) {
                const auto& a=weights.caches()[si];const auto& b=weightsRef.caches()[si];
                check(a.offsets()==b.offsets() && a.cellIndices()==b.cellIndices() && a.weights()==b.weights(),"explicit stage exact current-position CSR");
                for(std::size_t pi=0;pi<species[si].particle.size();++pi)
                    check(species[si].particle[pi].rngState==100+10*si+pi && species[si].particle[pi].trackId==200+10*si+pi,"field stage preserves particle identity and RNG");
            }
            for(int ci=0;ci<mesh.cellIndexNum();++ci) {
                check(actual.potential().value[ci]==reference.potential().value[ci]
                    && actual.poissonRhs().value[ci]==reference.poissonRhs().value[ci]
                    && actual.chargeDensity().value[ci]==reference.chargeDensity().value[ci]
                    && actual.chargeAmount().value[ci]==reference.chargeAmount().value[ci]
                    && deposit.electronAmount().value[ci]==depositRef.electronAmount().value[ci]
                    && actual.electricField().value[ci].x==reference.electricField().value[ci].x
                    && actual.electricField().value[ci].y==reference.electricField().value[ci].y,"explicit stage exact deposited and solved fields");
                if(step==2 && variant==0) {
                    const double x=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x;
                    near(actual.potential().value[ci],3.*x,"driven AMR analytic potential");
                    near(actual.electricField().value[ci].x,-3.,"driven AMR analytic electric field");
                }
            }
            check(deposit.totals().grid==depositRef.totals().grid && deposit.totals().particles==depositRef.totals().particles,"explicit stage charge totals");
            check(timing.buildWeights_s>=0. && timing.deposit_s>=0. && timing.rhs_s>=0. && timing.evalE_s>=0.
                && timing.amr_indicator_s==0.,"field timing resets without claiming indicator work");
        }
        const auto saved=weights.caches()[0].weights();const auto potential=actual.potential().value[0];
        const auto reject=[&](PicCurrentWeights& cache,ExplicitPicChargeDeposition& d,const Policy& policy) {
            bool caught=false;timing.rhs_s=-7.;const auto voltage=boundaries.domain().right.value;
            try {solver.advanceExplicit(graph,{species,macroWeight,norm,policy,embedded,geometry.get(),.75,5,true,0},
                cache,d,boundaries,report,&timing);}catch(const std::invalid_argument&){caught=true;}
            check(caught && timing.rhs_s==-7. && boundaries.domain().right.value==voltage
                && weights.caches()[0].weights()==saved && actual.potential().value[0]==potential,"mixed field-stage owners rejected before mutation");
        };
        reject(weights,depositForeign,coupling);reject(wrongWeights,deposit,coupling);
        reject(weights,deposit,wrongPolicy);reject(weights,deposit,implicitPolicy);
        // Numerical failure retains already completed deposition/boundary work,
        // but never evaluates E or pretends to roll the field stage back.
        if(variant==0) {
            FaceSegmentGraph empty;empty.build(mesh,[](const QuadtreeMesh&,const Cell&){return false;});
            solver.resetAfterMeshChange();BoundaryLogs logs;bool rejected=false;const auto electric=actual.electricField().value[0];
            try {solver.advanceExplicit(empty,{species,macroWeight,norm,coupling,embedded,nullptr,.75,6,true,0},
                weights,deposit,boundaries,report,&timing);}catch(const std::runtime_error&){rejected=true;}
            check(rejected && !report.success && timing.evalE_s==0.,"failed explicit solve stops before E evaluation");
            near(boundaries.domain().right.value,1.,"failed solve retains evaluated boundary time");
            check(actual.electricField().value[0].x==electric.x,"failed solve preserves preceding electric field");
        }
    }
}

void poissonStage()
{
    using Stage=PicPoissonFieldSolve;
    static_assert(!std::is_copy_constructible_v<Stage> && !std::is_move_constructible_v<Stage>);
    static_assert(std::is_same_v<decltype(std::declval<Stage&>().options()),const Stage::Options&>);
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    PicFieldState fields(mesh),uncachedFields(mesh);ImplicitPicFieldHistory history(fields),otherHistory(uncachedFields);
    Stage::Options options;options.cacheLogOnMiss=false;
    Stage solver(fields,options);options.cacheEnabled=false;Stage uncached(uncachedFields,options);
    check(solver.options().cacheEnabled && !uncached.options().cacheEnabled,"solver owns its prepared options");
    PicDomainBC bc;bc.left.value=0.;bc.right.value=1.;bc.bottom.type=bc.top.type=PicBCType::Neumann;
    Stage::Report report;std::vector<double> normals;
    const auto checkLinear=[&](double slope) {
        for(int ci=0;ci<mesh.cellIndexNum();++ci) {
            const double x=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x;
            near(fields.potential().value[ci],slope*x,"analytic linear potential");
            near(fields.electricField().value[ci].x,-slope,"analytic electric field");
            near(fields.electricField().value[ci].y,0.,"transverse electric field");
        }
    };
    solver.solveExplicit({graph,bc},0,report);
    check(report.success && report.N==16 && report.usedSuperLU && report.message.find("setup() completed")!=std::string::npos,"initial setup report preserves solver metadata");
    near(fields.electricField().value[0].x,0.,"potential solve does not advance electric stage");
    solver.evaluateCurrent({graph,bc},Stage::BoundaryPolicy::Configured,&normals);checkLinear(1.);
    check(normals.size()==graph.segments().size(),"normal field output follows current graph");
    bc.right.value=2.;solver.solveExplicit({graph,bc},1,report);
    check(report.success && report.message.find("cache hit")!=std::string::npos && report.message.find("setup() completed")==std::string::npos,"Dirichlet value change reuses factorization");
    solver.evaluateNext({graph,bc},history,&normals);
    near(fields.electricField().value[0].x,-1.,"next evaluation does not overwrite current layer");
    near(history.nextElectricField().value[0].x,-2.,"next electric layer from new potential");
    history.formHalfStep();near(history.halfElectricField().value[0].x,-1.5,"half layer remains separately timed");
    history.acceptNext();checkLinear(2.);
    bool rejected=false;try{solver.evaluateNext({graph,bc},otherHistory);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"foreign electric history rejected");near(otherHistory.nextElectricField().value[0].x,0.,"foreign history unchanged");
    solver.clearCache();solver.solveExplicit({graph,bc},2,report);
    check(report.success && report.message.find("cache miss")!=std::string::npos && report.message.find("setup() completed")==std::string::npos,"geometry invalidation retains original cached-path rebuild report");
    uncached.solveForRebuild({graph,bc});
    for(int ci=0;ci<16;++ci) {
        near(uncachedFields.potential().value[ci],fields.potential().value[ci],"uncached rebuild potential");
        near(uncachedFields.electricField().value[ci].x,-2.,"rebuild evaluates current field");
    }
    // Exercise nonzero RHS through the production deposition owner as well.
    const auto norm=NormalizationContext::MakeLegacy();
    std::vector<ParticleType> species;species.emplace_back("charge",1.,norm.scales().rho0_Cpm3);species[0].particle.resize(16);
    std::vector<int> offsets,indices;std::vector<double> weights;
    for(int ci=0;ci<16;++ci){species[0].particle[ci].w=1.;offsets.push_back(ci);indices.push_back(ci);weights.push_back(1.);}
    offsets.push_back(16);std::vector<CouplingWeightsCache> caches(1);caches[0].assignCSR(16,offsets,indices,weights);
    for(auto* target:{&fields,&uncachedFields}) {target->depositFullCellCharge(species,{},caches,norm);target->copyChargeDensityToRhs();}
    solver.solveExplicit({graph,bc},3,report);uncached.solveForRebuild({graph,bc});
    for(int ci=0;ci<16;++ci)near(fields.potential().value[ci],uncachedFields.potential().value[ci],"cached and uncached nonzero charge solve");
    check(fields.potential().value[5]>2.*mesh.cells(mesh.cellIndexToLabel[5]).centerCoord().x,"nonzero charge changes field solution");
    // Real AMR forces new fields, graph and setup; zero RHS recovers the same
    // linear solution on a uniformly refined mesh without a new discretization.
    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={1.,1.};box.level=1;mesh.applyRefineBoxes({box});
    check(mesh.cellIndexNum()==64,"real Poisson mesh refinement");graph.build(mesh);fields.resetAfterMeshChange();history.resetAfterMeshChange();solver.resetAfterMeshChange();
    solver.solveExplicit({graph,bc},4,report);check(report.success && report.N==64 && report.message.find("setup() completed")!=std::string::npos,"AMR restarts explicit setup bookkeeping");
    solver.evaluateCurrent({graph,bc},Stage::BoundaryPolicy::Configured);checkLinear(2.);
    FaceSegmentGraph empty;empty.build(mesh,[](const QuadtreeMesh&,const Cell&){return false;});BoundaryLogs logs;
    solver.resetAfterMeshChange();rejected=false;
    try{solver.solveExplicit({empty,bc},5,report);}catch(const std::runtime_error& error){rejected=std::string(error.what())=="Poisson setup failed";}
    check(rejected && !report.success,"failed setup preserves original exception and report");
    rejected=false;try{solver.solveForRebuild({empty,bc});}catch(const std::runtime_error& error){rejected=std::string(error.what()).find("[ImplicitPIC] Poisson solve failed:")==0;}
    check(rejected,"rebuild retains original failure prefix");
}

void chargeDeposition()
{
    static_assert(!std::is_copy_constructible_v<ExplicitPicChargeDeposition>);
    static_assert(!std::is_move_constructible_v<ExplicitPicChargeDeposition>);
    static_assert(std::is_same_v<decltype(std::declval<ExplicitPicChargeDeposition&>().electronAmount()),const CellField<double>&>);
    UnitsScales scales;scales.L0_m=2.;scales.T0_s=.5;scales.Phi0_V=3.;
    auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(scales));
    NormalizationConfig config;config.slabThickness_m=3.;
    const auto norm=NormalizationContext::MakeFromUnits(units,config);
    const double chargeUnit=norm.scales().rho0_Cpm3;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);PicFieldState fields(mesh),full(mesh);
    ExplicitPicChargeDeposition deposition(fields);
    PicEmbeddedGeometry noGeometry(PicEmbeddedGeometry::Mode::Explicit);
    std::vector<ParticleType> species;species.emplace_back("electron",1.,-chargeUnit);
    species.emplace_back("ion",2.,2.*chargeUnit);species.emplace_back("neutral",3.,0.);
    species[0].particle.resize(4);species[1].particle.resize(2);species[2].particle.resize(1);
    const double weights[]={2.,0.,-3.,std::numeric_limits<double>::quiet_NaN()};
    for(int p=0;p<4;++p){species[0].particle[p].w=weights[p];species[0].particle[p].rngState=20+p;}
    species[1].particle[0].w=4.;species[1].particle[1].w=0.;species[2].particle[0].w=7.;
    std::vector<CouplingWeightsCache> caches(3);
    caches[0].assignCSR(4,{0,1,2,3,4},{0,1,2,3},{1.,1.,1.,1.});
    caches[1].assignCSR(2,{0,2,3},{0,1,3},{.25,.75,1.});
    caches[2].assignCSR(1,{0,1},{2},{1.});
    const auto deposit=[&](bool number=true,int electron=0){deposition.deposit(species,{3.},caches,norm,noGeometry,nullptr,number,electron);};
    deposit();full.depositFullCellCharge(species,{3.},caches,norm);
    const double expectedCharge[]={0.,3.,-3.,-1.};const double expectedNumber[]={2.,3.,3.,3.};
    for(int ci=0;ci<4;++ci) {
        near(fields.chargeAmount().value[ci]/chargeUnit,expectedCharge[ci],"CSR species charge and default weights");
        near(fields.chargeDensity().value[ci],expectedCharge[ci]/3.,"2D physical volume and charge normalization");
        near(full.chargeDensity().value[ci],fields.chargeDensity().value[ci],"uncut explicit and rebuild closure");
        near(deposition.electronAmount().value[ci],expectedNumber[ci],"dual electron number deposition");
        near(fields.poissonRhs().value[ci],0.,"deposition does not advance RHS stage");
    }
    fields.copyChargeDensityToRhs();
    for(int ci=0;ci<4;++ci)near(fields.poissonRhs().value[ci],expectedCharge[ci]/3.,"RHS copies accepted deposited density");
    near(deposition.totals().particles/chargeUnit,-1.,"particle-order charge sum");
    near(deposition.totals().grid/chargeUnit,-1.,"grid charge sum");
    check(species[0].particle[0].rngState==20 && std::isnan(species[0].particle[3].w),"deposition does not alter particles or RNG");
    bool rejected=false;try{deposition.deposit(species,{3.},{},norm,noGeometry,nullptr,true,0);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"missing per-species cache rejected before field mutation");near(fields.chargeDensity().value[1],1.,"rejected input preserves charge");
    // A cache's own particle-count mismatch retains DepositQK's skip policy.
    caches[0].assignCSR(1,{0,1},{0},{1.});deposit();
    near(deposition.totals().particles/chargeUnit,-1.,"mismatched cache still counts macro charge");
    near(deposition.totals().grid/chargeUnit,10.,"mismatched cache skips only its deposited contribution");
    for(int ci=0;ci<4;++ci)near(deposition.electronAmount().value[ci],0.,"mismatched dual cache skips number");
    caches[0].assignCSR(4,{0,1,2,3,4},{0,1,2,3},{1.,1.,1.,1.});deposit();deposit(false);
    near(deposition.electronAmount().value[0],2.,"disabled AMR preserves unused prior number field");
    deposit(true,-1);for(int ci=0;ci<4;++ci)near(deposition.electronAmount().value[ci],0.,"no selected electron clears number field");
    for(auto& sp:species)sp.particle.clear();deposit();
    for(int ci=0;ci<4;++ci)near(fields.chargeAmount().value[ci],0.,"empty species clears old charge");
    near(deposition.totals().particles,0.,"empty particle totals");
    DictionaryMesh::RefineBox box;box.min={.05,.05};box.max={.3,.3};box.level=1;mesh.applyRefineBoxes({box});
    check(mesh.cellIndexNum()>4,"deposition contract changes real mesh topology");
    rejected=false;try{deposit();}catch(const std::logic_error&){rejected=true;}check(rejected,"old charge field rejected after AMR");
    fields.resetAfterMeshChange();rejected=false;try{deposit();}catch(const std::logic_error&){rejected=true;}
    check(rejected,"electron field must reset independently at its AMR stage");deposition.resetAfterMeshChange();
    species[0].particle.resize(1);auto& particle=species[0].particle[0];particle.w=5.;particle.position={.18,.18};particle.inCellL=mesh.findCellLabel(particle.position);
    PicCurrentWeights current(mesh);CouplingWeightsCache::BuildOptions options;options.levelSetI=-1;
    current.rebuild(species,{.1},options);
    deposition.deposit(species,{3.},current.caches(),norm,noGeometry,nullptr,true,0);
    double numberSum=0.,chargeSum=0.;
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        numberSum+=deposition.electronAmount().value[ci];chargeSum+=fields.chargeAmount().value[ci]/chargeUnit;
        near(fields.chargeDensity().value[ci],fields.chargeAmount().value[ci]/chargeUnit/(mesh.cells(mesh.cellIndexToLabel[ci]).area()*12.),"AMR density uses new cell volumes");
    }
    near(numberSum,5.,"number conserved on rebuilt coarse-fine cache");near(chargeSum,-5.,"charge conserved after real AMR");

    // Physical cut metrics override center classification; implicit rebuild
    // deliberately keeps the complete cell volume and deposits in solid cells.
    QuadtreeMesh cutMesh(input);PicFieldState cutFields(cutMesh),rebuild(cutMesh);ExplicitPicChargeDeposition cut(cutFields);
    DictionaryGeometry dict;dict.levelSetName={"solid"};dict.innerPoints={{.1,.5}};dict.moving={false};
    dict.linesOnInterface={{{0,false}}};dict.phaseName={"solid"};dict.levelSetOperation={{LevelSetOperation::non}};
    dict.levelSetIndexA={{{0,false}}};dict.levelSetIndexB={{}};dict.boundaryLines={{{0,false}}};
    Line line;line.points[0]={.625,0.};line.points[1]={.625,1.};dict.lines={line};
    DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.kind=GeometricBoundary::fixed;
    boundary.levelSetIndex=0;boundary.phase=0;boundary.lineIds={0};dict.boundaries={boundary};dict.boundaryById[1]=boundary;
    ImmersedGeometry<QuadtreeMesh> geometry(dict,cutMesh);
    geometry.levelSets[0].initialize(cutMesh,0.);
    geometry.LS_dx(0).initialize(cutMesh,0.);geometry.LS_dy(0).initialize(cutMesh,0.);
    for(int vi=0;vi<cutMesh.vertIndexNum();++vi) {
        geometry.levelSets[0].value[vi]=cutMesh.verts(cutMesh.vertIndexToLabel[vi]).coord().x-.625;
        geometry.LS_dx(0).value[vi]=1.;geometry.LS_dy(0).value[vi]=0.;
    }
    PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);embedded.resolveSelection(dict,1);embedded.rebuildComposite(cutMesh,&geometry);
    std::vector<ParticleType> cutSpecies;cutSpecies.emplace_back("positive",1.,chargeUnit);cutSpecies[0].particle.resize(4);
    for(auto& p:cutSpecies[0].particle)p.w=1.;std::vector<CouplingWeightsCache> cutCaches(1);
    cutCaches[0].assignCSR(4,{0,1,2,3,4},{0,1,2,3},{1.,1.,1.,1.});
    cut.deposit(cutSpecies,{},cutCaches,norm,embedded,&geometry,true,0);
    for(int ci=0;ci<4;++ci) {
        const bool solid=cutMesh.cells(cutMesh.cellIndexToLabel[ci]).centerCoord().x<.625;
        near(cutFields.chargeDensity().value[ci],solid?0.:1./3.,"center fallback before cut metrics");
    }
    FaceSegmentGraph graph;graph.build(cutMesh);embedded.rebuildMetrics(cutMesh,&geometry,graph,{1});
    check(embedded.metricsReady(),"real cut metrics constructed");
    cut.deposit(cutSpecies,{},cutCaches,norm,embedded,&geometry,true,0);rebuild.depositFullCellCharge(cutSpecies,{},cutCaches,norm);
    int partial=0,solid=0;
    for(int ci=0;ci<4;++ci) {
        const double alpha=embedded.alphaByCell()[ci];
        near(rebuild.chargeDensity().value[ci],1./3.,"implicit rebuild preserves full-cell-volume policy");
        if(alpha<=0.) {++solid;near(cutFields.chargeAmount().value[ci],0.,"solid charge discarded");near(cut.electronAmount().value[ci],0.,"solid number discarded");}
        else {if(alpha<1.)++partial;near(cutFields.chargeDensity().value[ci],1./(3.*alpha),"explicit density divides by physical cut volume");}
    }
    check(partial>0 && solid>0,"cut-volume and solid policies both exercised");
    near(cut.totals().particles/chargeUnit,4.,"pre-mask particle charge retained");
    near(cut.totals().grid/chargeUnit,4.-solid,"post-mask grid charge retained");
    PicPoissonFieldSolve::Options skip;skip.embeddedBCMissingPolicy=PicPoissonFieldSolve::Options::EmbeddedBCMissingPolicy::Skip;
    PicPoissonFieldSolve policySolver(cutFields,skip);PicDomainBC bc;
    // An entirely absent embedded closure is intentionally ignored by the
    // evaluator. Activate the map with an unrelated ID to exercise a missing
    // entry on the real cut boundary (ID 1), with no fallback.
    PicEmbeddedBCMap incomplete;incomplete[2].enabled=true;
    const PicPoissonFieldSolve::Input boundaryInput{graph,bc,&embedded.metrics(),nullptr,&incomplete};BoundaryLogs logs;
    policySolver.evaluateCurrent(boundaryInput,PicPoissonFieldSolve::BoundaryPolicy::Configured);
    bool missingRejected=false;
    try{policySolver.evaluateCurrent(boundaryInput,PicPoissonFieldSolve::BoundaryPolicy::LegacyDefault);}
    catch(const std::runtime_error& error){missingRejected=std::string(error.what()).find("Missing embedded BC")!=std::string::npos;}
    check(missingRejected,"implicit rebuild keeps default missing-BC policy while explicit uses configured skip");
}

void fieldBoundaryState()
{
    using State = PicFieldBoundaryState;
    static_assert(!std::is_copy_constructible_v<State> && !std::is_move_constructible_v<State>);
    static_assert(std::is_same_v<decltype(std::declval<State&>().domain()),const PicDomainBC&>);
    static_assert(std::is_same_v<decltype(std::declval<State&>().embeddedById()),const PicEmbeddedBCMap&>);
    constexpr double pi=3.141592653589793238462643383279502884;
    PicBoundarySide::TimeFunction drive;
    drive.enabled=true;drive.type="SiNe";drive.dc=3.;drive.ac=2.;drive.freq_hz=2.;drive.t0_s=.125;
    double radians=99.,degrees=99.;
    near(State::evaluate(drive,0.,&radians,&degrees),1.,"shifted sine voltage");
    near(radians,1.5*pi,"negative phase wraps to positive radians");near(degrees,270.,"phase degrees");
    near(State::evaluate(drive,.25),5.,"sine quarter cycle");
    drive.type="CoS";near(State::evaluate(drive,.125),5.,"cosine alias");
    drive.phase_rad=pi;near(State::evaluate(drive,.125),1.,"configured phase offset");drive.phase_rad=0.;
    drive.type="unknown";near(State::evaluate(drive,.25),5.,"legacy unknown waveform uses sine");
    for (const char* type : {"constant","CoNsT"}) {
        drive.type=type;near(State::evaluate(drive,.25,&radians,&degrees),3.,"constant ignores AC");
        near(radians,0.,"constant phase reset");near(degrees,0.,"constant degrees reset");
    }
    drive.type="sine";
    for(double frequency : {0.,-2.,std::numeric_limits<double>::quiet_NaN()}) {
        drive.freq_hz=frequency;near(State::evaluate(drive,.25,&radians,&degrees),3.,"invalid frequency fallback");
        near(radians,0.,"invalid frequency phase reset");
    }
    drive.freq_hz=2.;drive.ac=0.;near(State::evaluate(drive,.25),3.,"zero AC fallback");
    drive.ac=std::numeric_limits<double>::quiet_NaN();check(std::isnan(State::evaluate(drive,.25)),"NaN AC retains original propagation");
    drive.ac=2.;drive.enabled=false;near(State::evaluate(drive,.25),3.,"disabled raw waveform returns dc");
    PicDomainBC domain;domain.left.value=17.;domain.left.timeDirichlet=drive;
    near(State::voltage(domain.left,.25),17.,"disabled side uses stored value instead of waveform dc");
    drive.enabled=true;domain.left.timeDirichlet=drive;
    domain.right=domain.left;domain.right.type=PicBCType::Neumann;
    domain.bottom.value=12.;domain.top.type=PicBCType::Periodic;
    PicEmbeddedBC embedded;embedded.enabled=true;embedded.type=PicBCType::Dirichlet;embedded.g=19.;
    embedded.timeDirichlet=drive;embedded.timeDirichlet.type="cosine";
    PicEmbeddedBCMap byId{{7,embedded},{8,embedded},{9,embedded}};
    byId[7].timeDirichlet.type="sine";byId[8].enabled=false;byId[8].g=44.;
    byId[9].type=PicBCType::Robin;byId[9].g=99.;byId[9].a=2.;byId[9].b=3.;
    State state(domain,embedded,byId),independent(domain,embedded,byId);
    const auto* retainedDomain=&state.domain();const auto* retainedEmbedded=&state.embedded();
    const auto* retainedEntry=&state.embeddedById().at(7);
    near(state.domain().left.value,17.,"construction does not advance initial boundary time");
    domain.left.value=-100.;embedded.g=-100.;byId.clear();
    state.update(0.);
    near(retainedDomain->left.value,1.,"driven domain value updated in place");
    near(retainedDomain->right.value,17.,"Neumann ignores Dirichlet waveform");
    near(retainedDomain->bottom.value,12.,"static Dirichlet preserved");
    near(retainedEmbedded->g,3.,"default embedded cosine updated");near(retainedEntry->g,1.,"ID-specific sine updated");
    near(state.embeddedById().at(8).g,44.,"disabled embedded remains unchanged");
    near(state.embeddedById().at(9).g,99.,"Robin data preserved");near(state.embeddedById().at(9).b,3.,"Robin coefficient preserved");
    near(independent.domain().left.value,17.,"two boundary states do not share mutation");
    state.update(.125);near(retainedDomain->left.value,3.,"next physical time");near(retainedEmbedded->g,5.,"stable embedded address");
    check(retainedEntry==&state.embeddedById().at(7),"update does not replace map entries");
    double value=81.;check(state.embeddedVoltage(7,value),"current ID voltage available");near(value,3.,"current voltage is not reevaluated");
    check(state.embeddedVoltage(8,value),"disabled entry falls back");near(value,5.,"disabled entry uses default");
    check(state.embeddedVoltage(123,value),"missing entry falls back");near(value,5.,"missing entry uses default");
    value=81.;check(!state.embeddedVoltage(9,value),"enabled Robin overrides default without a voltage");near(value,81.,"missing voltage leaves output untouched");
    check(State::embeddedVoltage(state.embedded(),state.embeddedById(),7,.25,value),"explicit sampled voltage");near(value,5.,"explicit sampling uses requested time");
    near(retainedEntry->g,3.,"explicit sampling does not mutate stored implicit state");
    near(State::voltage(state.domain().right,.25),0.,"non-Dirichlet has no driven voltage");
    PicEmbeddedBC disabled;PicEmbeddedBCMap none;value=81.;
    check(!State::embeddedVoltage(disabled,none,1,.25,value),"no enabled electrode");near(value,81.,"no electrode preserves output");
}

void implicitAmrStage(const std::filesystem::path& directory)
{
    using Mode=PicDomain::Mode;using Policy=PicCurrentWeights::Policy;using Solver=PicPoissonFieldSolve;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;input.maxLevel_local={1};
    input.physicsRefinement.enabled=true;input.physicsRefinement.hysteresis=.2;
    input.physicsRefinement.indicators={"DebyeLength"};input.physicsRefinement.maxLevelCaps={{"debyelength",1}};
    const auto units=UnitsContext::MakeSI(UnitsScales{});const auto norm=NormalizationContext::MakeLegacy();
    const double density=PhysicalConstants::epsilon0_Fpm*PhysicalConstants::eV_J/(PhysicalConstants::qe_C*PhysicalConstants::qe_C*.01);
    for(int variant=0;variant<3;++variant)for(bool reset:{false,true}) {
        const bool fast=variant!=0,failOutput=variant==2;
        const auto run=[&](bool legacy) {
            const auto dir=directory/("implicit_amr_"+std::to_string(variant)+"_"+std::to_string(reset)+(legacy?"_legacy":"_stage"));
            std::filesystem::create_directories(dir);
            PicDomain domain(input,input,Mode::Implicit);auto& mesh=domain.mesh();domain.rebuildFaceGraph();
            PicFieldState fields(mesh);ImplicitPicFieldHistory history(fields);Solver poisson(fields,{});
            Policy coupling(mesh,DepositionConfig{},Policy::Mode::Implicit,false,false);
            PicSpeciesState particles(1,42);SpeciesConfig species;species.name="ion";species.m=1.;species.q=.02*norm.scales().rho0_Cpm3;species.numParticles=0;
            particles.initialize({species},norm.scales(),mesh,false);
            for(int i=0;i<3;++i){Particle p;p.position={.2+.2*i,.3};p.velocity={.01,0.};p.vz=2.;p.w=i+1.;p.rngState=700+i;p.trackId=800+i;
                p.inCellL=mesh.findCellLabel(p.position);particles.species()[0].particle.push_back(p);}
            particles.species()[0].number=3;const auto rng=particles.rng();
            PicCurrentWeights current(mesh);ImplicitPicOrbit orbit(mesh,1);ImplicitPicPolarization polarization(mesh,1,domain.faceGraph().segments().size());
            PicDomainBC bc;bc.right.value=2.;bc.top.type=bc.bottom.type=PicBCType::Neumann;
            const Solver::Input field{domain.faceGraph(),bc};Solver::Report solved;poisson.solveExplicit(field,0,solved);
            poisson.reseedImplicitFromPotential(field,history);check(solved.success && poisson.hasCache(),"implicit AMR starts from a solved field");
            const Solver::ImplicitRebuildInput rebuild{particles.species(),particles.weights(),norm,coupling,domain.embedded(),false};
            PicAmrGeometry amr(Mode::Implicit,{});PicMovingGeometry moving(Mode::Implicit);const std::vector<int> ids;
            ImplicitPicRefinement stage(mesh);ImplicitPicRefinement::Options options;
            options.enabled=options.wantDebye=true;options.electronIndex=0;options.temperatureEV=1.;options.indicatorEvery=1;
            std::vector<CellField<double>> amount;amount.emplace_back(mesh,density*.25);
            stage.update(options,{amount,fields.electricField(),nullptr,&units,norm},2,domain.refiner());
            check(stage.pending(),"implicit AMR fixture has refinement targets");
            std::vector<double> x,y;int callbacks=0;
            const auto coordinates=[&](const QuadtreeMesh& adapted) {
                ++callbacks;check(stage.pending() && &adapted==&mesh && fields.potential().size()==mesh.cellIndexNum()
                    && history.nextElectricField().size()==mesh.cellIndexNum() && !current.caches().empty()
                    && orbit.averagedCaches()[0].nnz()>0,"output hook follows field rebuild and precedes completion/transient reset");
                for(int ci=0;ci<mesh.cellIndexNum();++ci)check(history.nextElectricField().value[ci].x==fields.electricField().value[ci].x
                    && history.halfElectricField().value[ci].y==fields.electricField().value[ci].y,"AMR field history seeded before output");
                if(failOutput)throw std::runtime_error("injected coordinate refresh failure");
                x.clear();y.clear();for(int vi=0;vi<mesh.vertIndexNum();++vi){const auto pos=mesh.verts(mesh.vertIndexToLabel[vi]).coord();x.push_back(units.LengthCodeToSI(pos.x));y.push_back(units.LengthCodeToSI(pos.y));}
            };
            const auto bindings=[&]() {return ImplicitPicRefinement::Adaptation{true,fast,false,false,reset,2,domain,amr,moving,ids,
                fields,history,polarization,particles,poisson,coupling,current,orbit,field,rebuild,coordinates};};
            const auto prepareScratch=[&]() {
                current.rebuild(particles.species(),coupling,domain.embedded());orbit.gather(particles.species(),current.caches(),fields.electricField(),0);
                ImplicitPicOrbit::Options orbitOptions;orbitOptions.weights=coupling.physicalOptions(domain.embedded());
                orbit.build(particles.species(),coupling.halfWidths(),.001,orbitOptions,{0,false,false});
                check(orbit.averagedCaches()[0].nnz()>0,"AMR fixture has populated orbit scratch");
            };prepareScratch();
            if(!legacy && variant==0 && reset) {
                auto disabled=bindings();disabled.enabled=false;check(!stage.adaptIfNeeded(2,disabled) && !stage.adaptIfNeeded(1,bindings()),"implicit AMR disabled/cadence gates");
                disabled.enabled=true;disabled.every=0;check(!stage.adaptIfNeeded(2,disabled),"implicit AMR nonpositive cadence gate");
                PicDomain other(input,input,Mode::Implicit);other.rebuildFaceGraph();PicFieldState otherFields(other.mesh());ImplicitPicFieldHistory otherHistory(otherFields);
                Solver otherSolver(otherFields,{});Policy wrongPolicy(mesh,DepositionConfig{},Policy::Mode::Explicit,false,false);
                PicCurrentWeights otherCurrent(other.mesh());ImplicitPicOrbit otherOrbit(other.mesh(),1);ImplicitPicPolarization otherPolarization(other.mesh(),1,other.faceGraph().segments().size());
                PicAmrGeometry wrongGeometry(Mode::Explicit,{});PicMovingGeometry wrongMoving(Mode::Explicit);
                auto copiedSpecies=particles.species();Solver::ImplicitRebuildInput wrongRebuild{copiedSpecies,particles.weights(),norm,coupling,domain.embedded(),false};
                for(int fault=1;fault<=13;++fault) {
                    ImplicitPicRefinement::Adaptation bad{true,fast,false,false,reset,2,fault==1?other:domain,
                        fault==11?wrongGeometry:amr,fault==12?wrongMoving:moving,ids,fault==2?otherFields:fields,
                        fault==3?otherHistory:history,fault==4?otherPolarization:polarization,particles,fault==5?otherSolver:poisson,
                        fault==6?wrongPolicy:coupling,fault==7?otherCurrent:current,fault==8?otherOrbit:orbit,
                        {fault==9?other.faceGraph():domain.faceGraph(),bc},fault==10?wrongRebuild:rebuild,coordinates};
                    if(fault==13)bad.refreshOutputCoordinates={};
                    bool rejected=false;try{stage.adaptIfNeeded(2,bad);}catch(const std::logic_error&){rejected=true;}
                    check(rejected && mesh.cellIndexNum()==4 && stage.pending() && callbacks==0 && poisson.hasCache() && particles.rng()==rng,
                        "invalid implicit AMR binding rejects before mutation and output");
                }
            }
            std::map<std::string,std::string> caches;
            for(int pass=0;pass<2;++pass) {
                const int step=2+2*pass;
                if(pass) {
                    auto coarseInput=input;coarseInput.physicsRefinement.maxLevelCaps={{"debyelength",0}};QuadtreeMeshRefiner coarse(coarseInput);
                    amount[0].initialize(mesh,0.);for(int ci=0;ci<mesh.cellIndexNum();++ci)amount[0].value[ci]=density*mesh.cells(mesh.cellIndexToLabel[ci]).area();
                    stage.update(options,{amount,fields.electricField(),nullptr,&units,norm},step,coarse);stage.requestAdaptation();prepareScratch();
                }
                bool failed=false;
                try {
                    if(!legacy)check(stage.adaptIfNeeded(step,bindings()),"implicit AMR stage runs");
                    else {
                        const int before=mesh.cellIndexNum();CellField<double> mapped;if(fast)mapped=fields.potential();amr.captureBeforeAdapt(nullptr);
                        std::vector<VertField<double>*> vertices;stage.adaptMesh(vertices,fast?&mapped:nullptr);
                        const double change=static_cast<double>(std::abs(mesh.cellIndexNum()-before))/static_cast<double>(std::max(1,before));
                        amr.refresh(mesh,nullptr,domain.embedded(),nullptr,{step,change,false,false,domain.refiner().phaseActive(),false});
                        domain.rebuildFaceGraph();fields.resetAfterMeshChange(fast?&mapped:nullptr);history.resetAfterMeshChange();
                        polarization.resetDepositedFields();stage.resetFieldsAfterMeshChange();polarization.resetCoefficients(domain.faceGraph().segments().size());
                        coupling.refreshAfterMeshChange();particles.relabelAfterMeshChange(mesh,{true,true});
                        if(fast)poisson.reseedImplicitFromPotential(field,history);else{PicCurrentWeights rebuilt(mesh);poisson.rebuildImplicit(field,rebuild,rebuilt,&history);}
                        coordinates(mesh);stage.adaptationCompleted();if(reset){current.clear();orbit.resetScratch();polarization.resetScratch();}
                    }
                }catch(const std::runtime_error& error){check(failOutput && std::string(error.what())=="injected coordinate refresh failure","unexpected AMR failure");failed=true;}
                check(failed==failOutput && mesh.cellIndexNum()==(pass?4:16) && stage.pending()==failOutput && particles.rng()==rng,"implicit AMR real topology transitions, completion and RNG");
                if(!legacy && failed)check(stage.adaptationReport().cellsBefore==0 && stage.adaptationReport().seconds==0.,"failed adaptation does not publish a completed report");
                check(polarization.chargeAmount().size()==mesh.cellIndexNum() && polarization.faceCoefficients().size()==domain.faceGraph().segments().size(),"polarization storage follows adapted mesh");
                if(reset && !failed)check(current.caches().empty() && orbit.averagedCaches()[0].nnz()==0,"optional transient reset completes after output");
                else check(!current.caches().empty() && orbit.averagedCaches()[0].nnz()>0,"disabled or failed completion retains transient caches");
                if(fast)for(int ci=0;ci<mesh.cellIndexNum();++ci){near(fields.electricField().value[ci].x,-2.,"analytic fast-mapped E");check(fields.chargeAmount().value[ci]==0.,"fast remap skips charge rebuild");}
                else {double q=0.;for(int ci=0;ci<mesh.cellIndexNum();++ci)q+=fields.chargeAmount().value[ci];check(q>0.,"ordinary AMR rebuild deposits charge");}
                if(!legacy && !failed){const auto saved=stage.adaptationReport();check(saved.cellsBefore==(pass?16:4) && saved.cellsAfter==(pass?4:16)
                    && saved.useFastRemap==fast && saved.seconds>=0.,"implicit adaptation report");check(!stage.adaptIfNeeded(step,bindings()) && stage.adaptationReport().seconds==saved.seconds,"implicit skipped report is retained");}
                const auto saved=dir/std::to_string(step);std::filesystem::create_directories(saved);const auto prefix=PicCheckpointStore::pathPrefix(saved);
                mesh.writeCache(prefix);particles.species()[0].writeCache(prefix);history.writeCheckpoint(saved);
                for(const auto& entry:std::filesystem::directory_iterator(saved))if(entry.path().extension()==".bin") {std::ifstream f(entry.path(),std::ios::binary);caches[std::to_string(step)+entry.path().filename().string()]={std::istreambuf_iterator<char>(f),{}};}
                std::ostringstream coordinatesText;coordinatesText.precision(17);for(std::size_t i=0;i<x.size();++i)coordinatesText<<x[i]<<','<<y[i]<<';';caches[std::to_string(step)+"coordinates"]=coordinatesText.str();
                if(failed)break;
            }
            return caches;
        };
        check(run(true)==run(false),"implicit AMR exact original mesh, particle, seven field layers and output coordinates");
    }
}

void explicitAmrStage(const std::filesystem::path& directory)
{
    using Mode=PicDomain::Mode;
    const auto units=UnitsContext::MakeSI(UnitsScales{});
    const auto norm=NormalizationContext::MakeLegacy();
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;
    input.maxLevel_local={1};input.physicsRefinement.enabled=true;input.physicsRefinement.hysteresis=.2;
    input.physicsRefinement.indicators={"DebyeLength"};input.physicsRefinement.maxLevelCaps={{"debyelength",1}};
    const double density=PhysicalConstants::epsilon0_Fpm*PhysicalConstants::eV_J/
        (PhysicalConstants::qe_C*PhysicalConstants::qe_C*.01);
    for(bool embeddedEnabled:{false,true}) {
        const auto run=[&](bool legacy) {
            const auto dir=directory/(std::string(embeddedEnabled?"amr_wall_":"amr_domain_")+(legacy?"legacy":"stage"));
            std::filesystem::create_directories(dir);
            PicDomain domain(input,input,Mode::Explicit);auto& mesh=domain.mesh();auto& embedded=domain.embedded();
            DictionaryGeometry dict;dict.levelSetName={"solid"};dict.innerPoints={{.5,.5}};dict.moving={true};
            dict.linesOnInterface={{{0,false},{1,false},{2,false},{3,false}}};
            dict.phaseName={"solid","gas"};dict.levelSetOperation={{LevelSetOperation::non},{LevelSetOperation::non}};
            dict.levelSetIndexA={{{0,false}},{{0,true}}};dict.levelSetIndexB=dict.levelSetIndexA;
            dict.boundaryLines={{{0,false},{1,false},{2,false},{3,false}},{{0,true},{1,true},{2,true},{3,true}}};
            const std::array<Vector2D,4> corners{{{.3,.3},{.7,.3},{.7,.7},{.3,.7}}};
            for(int i=0;i<4;++i){Line line;line.points[0]=corners[i];line.points[1]=corners[(i+1)%4];line.levelSetIndex=0;dict.lines.push_back(line);}
            DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.kind=GeometricBoundary::fixed;
            boundary.levelSetIndex=0;boundary.phase=0;boundary.lineIds={0,1,2,3};dict.boundaries={boundary};dict.boundaryById[1]=boundary;
            bool eb=embeddedEnabled,refine=false;
            check(domain.initializeGeometry(embeddedEnabled?&dict:nullptr,eb,refine,{}),"AMR geometry setup");
            PicMovingGeometry moving(Mode::Explicit);MovingEBConfig motion;motion.enabled=embeddedEnabled;
            motion.metricsPath=(dir/"metrics.csv").string();motion.erosionMetricsPath=(dir/"regions.csv").string();
            motion.impactStatsPath=(dir/"impacts.csv").string();
            MovingEBRegion region;region.name="wall";region.boundaryId=1;region.levelSetIndex=0;
            region.selector.wholeBoundary=true;region.motionMode="rigidShift";region.velocityX_mps=.01;motion.regions={region};
            check(moving.initialize(mesh,domain.geometry(),embeddedEnabled?&dict:nullptr,motion,&units,&norm,std::nullopt),"AMR moving owner setup");
            check(!embeddedEnabled || moving.controller()!=nullptr,"AMR wall fixture has a live moving controller");
            const std::vector<int> ids(4,1);domain.rebuildFaceGraph();
            if(eb)embedded.rebuildMetrics(mesh,domain.geometry(),domain.faceGraph(),ids,moving.controller());
            PicFieldState fields(mesh);ExplicitPicChargeDeposition deposition(fields);PicPoissonFieldSolve poisson(fields,{});
            PicCurrentWeights::Policy coupling(mesh,DepositionConfig{},PicCurrentWeights::Policy::Mode::Explicit,false,false);
            PicSpeciesState particles(1,24);SpeciesConfig species;species.name="electron";species.m=1.;species.q=-1e-12;species.numParticles=0;
            particles.initialize({species},norm.scales(),mesh,false);
            for(auto position:{Vector2D(.1,.2),Vector2D(.6,.8),Vector2D(1.2,.5)}) {
                Particle p;p.position=position;p.velocity={1.,2.};p.vz=3.;p.w=2.;p.rngState=23;p.trackId=45;
                p.inCellL=mesh.findCellLabel(position);particles.species()[0].particle.push_back(p);
            }
            particles.species()[0].number=3;
            const auto initialRng=particles.rng();
            PicDomainBC bc;bc.left.value=1.;PicPoissonFieldSolve::Report solved;
            poisson.solveExplicit({domain.faceGraph(),bc},0,solved);poisson.evaluateCurrent({domain.faceGraph(),bc},PicPoissonFieldSolve::BoundaryPolicy::Configured);
            check(solved.success && poisson.hasCache(),"AMR fixture starts with a real cached solve");
            ExplicitPicRefinement stage(mesh);CellField<double> amount(mesh,density*.25);
            AmrConfig config;config.debye.neFloor=1.;config.debye.avgWindowSteps=1;config.indicatorEverySteps=1;config.persistenceSteps=1;config.adaptEverySteps=2;
            check(stage.update(config,{amount,fields.electricField(),nullptr,units,norm},1.,2,8,domain.refiner()) && stage.pending(),"AMR fixture requests real refinement");
            PicAmrGeometry amr(Mode::Explicit,{});
            const auto state=[&]() {return ExplicitPicRefinement::Adaptation{true,eb,false,true,domain,amr,moving,ids,fields,deposition,particles,poisson,coupling};};
            if(!legacy) {
                // Rebind one participant at a time; all rejections precede mesh mutation.
                PicDomain other(input,input,Mode::Explicit);PicFieldState otherFields(other.mesh());
                ExplicitPicChargeDeposition otherDeposition(otherFields);PicPoissonFieldSolve otherPoisson(otherFields,{});
                PicCurrentWeights::Policy otherPolicy(other.mesh(),DepositionConfig{},PicCurrentWeights::Policy::Mode::Explicit,false,false);
                PicCurrentWeights::Policy implicitPolicy(mesh,DepositionConfig{},PicCurrentWeights::Policy::Mode::Implicit,false,false);
                PicAmrGeometry implicitGeometry(Mode::Implicit,{});PicMovingGeometry implicitMoving(Mode::Implicit);PicSpeciesState empty(0,1);
                auto disabled=state();disabled.enabled=false;
                check(!stage.adaptIfNeeded(config,2,disabled) && !stage.adaptIfNeeded(config,1,state()),"disabled and off-cadence AMR do nothing");
                auto disabledConfig=config;disabledConfig.adaptEverySteps=0;
                check(!stage.adaptIfNeeded(disabledConfig,2,state()),"nonpositive AMR cadence does nothing");
                for(int fault=1;fault<=9;++fault) {
                    ExplicitPicRefinement::Adaptation invalid{true,eb,false,true,fault==1?other:domain,
                        fault==7?implicitGeometry:amr,fault==8?implicitMoving:moving,ids,
                        fault==2?otherFields:fields,fault==3?otherDeposition:deposition,fault==9?empty:particles,
                        fault==4?otherPoisson:poisson,fault==5?otherPolicy:(fault==6?implicitPolicy:coupling)};
                    bool rejected=false;try{stage.adaptIfNeeded(config,2,invalid);}catch(const std::logic_error&){rejected=true;}
                    check(rejected && mesh.cellIndexNum()==4 && stage.pending() && poisson.hasCache()
                        && particles.rng()==initialRng && fields.potential().size()==4,"invalid AMR binding rejects before any transition");
                }
            }
            const auto oldSequence=[&](int step) {
                std::vector<VertField<double>*> vertices;amr.captureBeforeAdapt(domain.geometry());
                if(domain.geometry())for(auto& field:domain.geometry()->levelSets)vertices.push_back(&field);
                stage.adaptMesh(vertices);if(moving.controller())moving.onMeshChanged(mesh);
                amr.refresh(mesh,domain.geometry(),embedded,moving.controller(),{step,0.,eb,false,domain.refiner().phaseActive(),false});
                domain.rebuildFaceGraph();if(embedded.ready())embedded.rebuildMetrics(mesh,domain.geometry(),domain.faceGraph(),ids,moving.controller());
                fields.resetAfterMeshChange();deposition.resetAfterMeshChange();
                particles.relabelAfterMeshChange(mesh,PicSpeciesState::RelabelOptions::fromExplicitConfig(config.particleRelabel.mode,config.particleRelabel.fallback));
                poisson.resetAfterMeshChange();coupling.refreshAfterMeshChange();stage.adaptationCompleted();
            };
            if(legacy)oldSequence(2);else check(stage.adaptIfNeeded(config,2,state()),"complete explicit AMR stage runs");
            check(mesh.cellIndexNum()==16 && !stage.pending() && !poisson.hasCache()
                && fields.potential().size()==16 && deposition.electronAmount().size()==16 && coupling.minimumLeafWidth()==.25,
                "real 4 to 16 AMR resets fields and cache and refreshes coupling width");
            check(!eb || (embedded.metricsReady() && embedded.alphaByCell().size()==16
                && domain.geometry()->levelSets[0].size()==mesh.vertIndexNum()),"AMR transfers level sets and rebuilds cut metrics");
            for(int ci=0;ci<16;++ci)check(fields.potential().value[ci]==0. && fields.chargeAmount().value[ci]==0.
                && fields.electricField().value[ci].x==0. && deposition.electronAmount().value[ci]==0.,"explicit AMR discards the old fields without solving");
            check(particles.species()[0].number==3 && particles.species()[0].particle[2].inCellL==-1 && particles.rng()==initialRng,"AMR relabel does not cull out-of-domain particles or draw RNG");
            if(!legacy) {
                const auto saved=stage.adaptationReport();
                check(saved.seconds>=0. && saved.levelsBefore==std::vector<int>{4,0} && saved.levelsAfter==std::vector<int>{0,16},"AMR before/after level report");
                check(!stage.adaptIfNeeded(config,4,state()) && stage.adaptationReport().seconds==saved.seconds,"no pending adaptation retains report");
            }
            // A second real transition coarsens the transferred target/history.
            auto coarseInput=input;coarseInput.physicsRefinement.maxLevelCaps={{"debyelength",0}};
            QuadtreeMeshRefiner coarse(coarseInput);amount.initialize(mesh,0.);
            for(int ci=0;ci<16;++ci)amount.value[ci]=density*mesh.cells(mesh.cellIndexToLabel[ci]).area();
            stage.update(config,{amount,fields.electricField(),nullptr,units,norm},1.,4,8,coarse);stage.requestAdaptation();
            if(legacy)oldSequence(4);else {auto noReport=state();noReport.reportLevels=false;check(stage.adaptIfNeeded(config,4,noReport),"second AMR stage runs");}
            check(mesh.cellIndexNum()==4 && coupling.minimumLeafWidth()==.5 && particles.rng()==initialRng,"real 16 to 4 coarsening retains RNG and refreshes width");
            if(!legacy)check(stage.adaptationReport().levelsAfter.empty(),"disabled console skips after-level collection");
            const std::string prefix=PicCheckpointStore::pathPrefix(dir);
            mesh.writeCache(prefix);domain.writeGeometryCache(prefix);moving.writeCache(dir);
            particles.species()[0].writeCache(prefix);ImplicitPicFieldHistory history(fields);history.writeCheckpoint(dir);
            std::map<std::string,std::string> result;
            for(const auto& entry:std::filesystem::directory_iterator(dir))if(entry.path().extension()==".bin") {
                std::ifstream file(entry.path(),std::ios::binary);result[entry.path().filename().string()]={std::istreambuf_iterator<char>(file),{}};
            }
            return result;
        };
        check(run(true)==run(false),"AMR mesh, geometry, particle and field caches exactly match the original sequence");
    }
}

void refinementState()
{
    static_assert(!std::is_copy_constructible_v<ExplicitPicRefinement>);
    static_assert(!std::is_move_constructible_v<ImplicitPicRefinement>);
    static_assert(std::is_same_v<decltype(std::declval<ExplicitPicRefinement&>().targetLevels().value[0]),const double&>);
    const auto units=UnitsContext::MakeSI(UnitsScales{});
    const auto norm=NormalizationContext::MakeLegacy();
    const double density=PhysicalConstants::epsilon0_Fpm*PhysicalConstants::eV_J/
        (PhysicalConstants::qe_C*PhysicalConstants::qe_C*.01); // lambda_D = 0.1 m at 1 eV
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;
    input.maxLevel_local={1};input.physicsRefinement.enabled=true;input.physicsRefinement.hysteresis=.2;
    input.physicsRefinement.indicators={"DebyeLength","ElectricField"};
    input.physicsRefinement.maxLevelCaps={{"debyelength",1},{"electricfield",1}};
    QuadtreeMesh mesh(input);QuadtreeMeshRefiner refiner(input);ExplicitPicRefinement stage(mesh);
    CellField<double> number(mesh,density*.25);CellField<Vector2D> field(mesh,{2.,0.});
    AmrConfig config;config.debye.neFloor=1.;config.debye.avgWindowSteps=3;
    config.indicator.eFieldAvgWindow=3;config.indicator.eFieldRefAvgWindow=3;
    config.indicatorEverySteps=2;config.persistenceSteps=2;
    const auto update=[&](int step) {return stage.update(config,{number,field,nullptr,units,norm},1.,step,20,refiner);};
    check(!update(1) && !stage.pending() && stage.indicators()[0].name.empty(),"EMA-only step does not publish targets");
    for(int ci=0;ci<4;++ci){number.value[ci]*=4.;field.value[ci]={6.,0.};}
    check(update(2) && !stage.pending(),"first heavy update starts persistence");
    near(stage.report().lambdaD_min_m,.1/std::sqrt(2.5),"density EMA before heavy update");
    near(stage.indicators()[1].cellValue[0],4.,"electric EMA");
    check(!update(3) && !stage.pending(),"skipped heavy update preserves persistence");
    check(update(4) && stage.pending(),"second heavy update reaches persistence");
    near(stage.report().lambdaD_min_m,.1/std::sqrt(3.625),"density history advances on skipped heavy steps");
    near(stage.indicators()[1].reference,4.75,"reference EMA uses only heavy updates");
    std::vector<VertField<double>*> vertices;stage.adaptMesh(vertices);
    check(mesh.cellIndexNum()==16 && stage.pending(),"real explicit refinement retains pending until completion");
    stage.adaptationCompleted();check(!stage.pending(),"explicit completion resets pending");
    number.initialize(mesh,0.);field.initialize(mesh,{0.,0.});
    for(int ci=0;ci<mesh.cellIndexNum();++ci)number.value[ci]=density*mesh.cells(mesh.cellIndexToLabel[ci]).area();
    config.indicatorEverySteps=1;update(5);
    near(stage.report().lambdaD_min_m,.1/std::sqrt(2.3125),"AMR transfers density EMA instead of restarting it");
    near(stage.indicators()[1].cellValue[0],2.75,"AMR transfers electric EMA");
    near(stage.indicators()[1].reference,3.75,"AMR retains reference history");
    stage.requestAdaptation();check(stage.pending(),"moving geometry can request explicit adaptation");

    // Uniform 2x2 support sees all four cells. Use effective-volume weights;
    // the zero-alpha cell must not contribute to either smoothed field.
    QuadtreeMesh smoothMesh(input);QuadtreeMeshRefiner smoothRefiner(input);ExplicitPicRefinement smooth(smoothMesh);
    CellField<double> amounts(smoothMesh,0.);CellField<Vector2D> electric(smoothMesh,{0.,0.});
    std::vector<double> alpha{1.,.5,0.,std::numeric_limits<double>::quiet_NaN()};
    for(int ci=0;ci<4;++ci) {amounts.value[ci]=density*(ci+1)*(ci+1)*.25*(ci==1?.5:1.);electric.value[ci]={2.*(ci+1),0.};}
    config.debye.smoothWindowCells=1;config.indicator.eFieldSmoothWindowCells=1;config.indicator.eFieldAbsFloor=10.;
    smooth.update(config,{amounts,electric,&alpha,units,norm},1.,1,10,smoothRefiner);
    for(int ci=0;ci<4;++ci) {
        near(smooth.indicators()[0].cellValue[ci],ci==2?0.:.06,"effective-volume smoothing of Debye length");
        near(smooth.indicators()[1].cellValue[ci],ci==2?0.:4.8,"effective-volume smoothing of electric magnitude");
    }
    near(smooth.indicators()[1].reference,10.,"electric reference absolute floor");

    QuadtreeMesh implicitMesh(input);QuadtreeMeshRefiner implicitRefiner(input);ImplicitPicRefinement implicit(implicitMesh);
    std::vector<CellField<double>> numbers{CellField<double>(implicitMesh,density*.25)};
    CellField<Vector2D> implicitField(implicitMesh,{0.,0.});CellField<double> potential(implicitMesh,7.);
    for(int ci=0;ci<4;++ci)implicitField.value[ci]={static_cast<double>(1<<ci),0.};
    ImplicitPicRefinement::Options options;options.enabled=options.wantDebye=options.wantEField=true;
    options.electronIndex=0;options.temperatureEV=1.;options.electricPercentile=50.;options.maxChanges=1;options.persistenceSteps=2;
    const auto updateImplicit=[&](int step) {implicit.update(options,{numbers,implicitField,nullptr,&units,norm},step,implicitRefiner);};
    updateImplicit(1);check(!implicit.pending(),"implicit persistence first update");
    for(int ci=0;ci<4;++ci)check(implicit.targetLevels().value[ci]==(ci==3?1.:0.),"implicit budget selects highest combined indicator");
    updateImplicit(2);check(implicit.pending(),"implicit persistence second update");
    options.gateElectricRelative=options.gateDebyeRelative=1e-6;updateImplicit(3);
    check(implicit.skippedByGate() && !implicit.pending(),"unchanged indicators gate and clear pending persistence");
    implicit.requestAdaptation();check(implicit.pending(),"geometry request overrides pending decision");
    implicit.adaptMesh(vertices,&potential);check(implicitMesh.cellIndexNum()==7,"implicit budget refines only one leaf");
    for(int ci=0;ci<potential.size();++ci)near(potential.value[ci],7.,"optional fast-remap potential transfer");
    implicit.resetFieldsAfterMeshChange();implicit.adaptationCompleted();
    numbers[0].initialize(implicitMesh,0.);implicitField.initialize(implicitMesh,{2.,0.});
    for(int ci=0;ci<implicitMesh.cellIndexNum();++ci)numbers[0].value[ci]=density*implicitMesh.cells(implicitMesh.cellIndexToLabel[ci]).area();
    updateImplicit(4);check(!implicit.skippedByGate() && !implicit.pending(),"changed mesh rebuilds indicator history and restarts persistence");
    updateImplicit(5);check(implicit.skippedByGate(),"rebuilt indicators can gate again");
    options.enabled=false;updateImplicit(6);check(!implicit.skippedByGate(),"inactive update clears per-step gate marker");
    options.enabled=true;options.wantDebye=options.wantEField=false;BoundaryLogs logs;updateImplicit(7);
    check(implicit.indicators()[0].name.empty() && implicit.indicators()[1].name.empty(),"disabled physics clears owned indicator histories");
}
void refinementInputs()
{
    static_assert(std::is_same_v<decltype(std::declval<PicRefinementIndicators&>().fields()[0]),const RefinementIndicatorField&>);
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=64;input.maxInitialIterations=0;
    input.maxLevel_local={2};input.physicsRefinement.enabled=true;
    input.physicsRefinement.indicators={"dEbYeLeNgTh","ElectricFIELD"};input.physicsRefinement.hysteresis=.2;
    input.physicsRefinement.maxLevelCaps={{"debyelength",2},{"electricfield",2}};
    QuadtreeMesh mesh(input),otherMesh(input);QuadtreeMeshRefiner refiner(input),otherRefiner(input);
    CellField<double> target(mesh,0.),score(mesh,0.),otherTarget(otherMesh,0.),otherScore(otherMesh,0.);
    PicRefinementIndicators dense,dilute;
    std::vector<double> values(4,.01);dense.setDebye(values,.5);values.assign(4,100.);
    dilute.setDebye(values,.5);BoundaryLogs logs;
    const auto build=[&](const PicRefinementIndicators* supplied) {
        return refiner.buildTargetLevels(mesh,target,&score,nullptr,false,nullptr,supplied?&supplied->fields():nullptr);
    };
    check(build(&dense),"dense Debye input must produce actual refinement targets");
    for(int ci=0;ci<4;++ci) {near(score.value[ci],100.,"Debye ratio");check(target.value[ci]==2.,"Debye cap");}
    otherRefiner.buildTargetLevels(otherMesh,otherTarget,&otherScore,nullptr,false,nullptr,&dilute.fields());
    for(int ci=0;ci<4;++ci)check(otherTarget.value[ci]==0. && target.value[ci]==2.,"live refiners must not use another session's indicators");
    check(!build(nullptr),"omitting indicators must not reuse the preceding call's inputs");
    for(int ci=0;ci<4;++ci)check(target.value[ci]==0. && score.value[ci]==0.,"missing inputs preserve original geometry-only fallback");
    build(&dense);check(target.value[0]==2.,"another refiner and an omitted call must not mutate owned indicators");
    dense.setDebye(std::vector<double>(4,1.),.5);dense.setEField({1.,2.,4.,8.},2.);
    build(&dense);
    for(int ci=0;ci<4;++ci) {
        const double expectedScore=ci<2?1.:static_cast<double>(1<<(ci-1));
        near(score.value[ci],expectedScore,"combined Debye/E ratio");
        check(target.value[ci]==(ci<2?0.:ci-1.),"combined indicator target levels");
    }
    dense.clearDebye();dense.setEField({1.,2.,4.,8.},-1.,50.);build(&dense);
    for(int ci=0;ci<4;++ci)near(score.value[ci],static_cast<double>(1<<ci)/4.,"electric percentile reference");
    dense.clearEField();check(!build(&dense),"clearing both indicators makes them unavailable");
    dense.setDebye({.01},.5);check(!build(&dense),"wrong-sized input retains warning/ignore behavior");
    check(logs.errors.str().find("indicator size mismatch")!=std::string::npos,"wrong-sized input warning");
    QuadtreeMeshRefiner::CurrentFields flow;target.value[0]=17.;bool rejected=false;
    try {refiner.buildTargetLevels(mesh,target,&score,nullptr,true,&flow,&dense.fields());}
    catch(const std::invalid_argument& e) {rejected=std::string(e.what()).find("must not be combined")!=std::string::npos;}
    check(rejected && target.value[0]==17.,"mixed flow/PIC input rejects before changing target storage");
    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={.4,.4};box.level=1;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()>4,"refinement input contract changes real topology");
    dense.setDebye(std::vector<double>(mesh.cellIndexNum(),.01),.5);build(&dense);
    check(target.size()==mesh.cellIndexNum() && score.size()==mesh.cellIndexNum(),"target storage follows rebuilt mesh inputs");
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        const auto& c=mesh.cells(mesh.cellIndexToLabel[ci]);
        near(score.value[ci],std::min(c.size().x,c.size().y)/.005,"Debye score after AMR");
        check(target.value[ci]==2.,"new mesh target cap");
    }
    check(dilute.debyeValues()==values && dilute.electricFieldValues().empty(),"AMR in another session leaves owned inputs intact");
}
void currentWeights()
{
    static_assert(!std::is_copy_constructible_v<PicCurrentWeights>);
    static_assert(!std::is_move_constructible_v<PicCurrentWeights>);
    static_assert(std::is_same_v<decltype(std::declval<PicCurrentWeights&>().caches()[0]),const CouplingWeightsCache&>);
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);PicCurrentWeights weights(mesh);
    const auto& caches=weights.caches();
    std::vector<ParticleType> species;
    species.emplace_back("large",1.,-1.);species.emplace_back("small",2.,1.);species.emplace_back("empty",1.,0.);
    for(int si=0;si<2;++si) {
        species[si].particle.resize(si==0?300:7);
        for(std::size_t j=0;j<species[si].particle.size();++j) {
            auto& p=species[si].particle[j];p.position={j%2?.97:.03,.35+.1*(j%3)};
            p.inCellL=mesh.findCellLabel(p.position);p.rngState=123+j;p.w=2.;
        }
    }
    std::vector<int> ignoredHosts(300,-1);std::vector<double> alpha(mesh.cellIndexNum(),1.);
    CouplingWeightsCache::BuildOptions options;options.periodicX=options.periodicY=true;
    options.levelSetI=-1;options.alphaByCellIndex=&alpha;options.cutByAlpha=true;
    options.enableMomentExactCut=false;options.hostCellLabels=&ignoredHosts;
    const std::vector<double> hp{.2};
    const auto verify=[&](bool parallel) {
        const auto counts=weights.rebuild(species,hp,options,parallel);
        std::uint64_t fallback=0,particles=0;
        check(caches.size()==species.size(),"current caches must follow current species count");
        for(std::size_t si=0;si<species.size();++si) {
            std::vector<Vector2D> positions;std::vector<int> hosts;
            for(const auto& p:species[si].particle) {positions.push_back(p.position);hosts.push_back(p.inCellL);}
            auto referenceOptions=options;referenceOptions.hostCellLabels=&hosts;
            CouplingWeightsCache expected;expected.build(mesh,positions,hp,referenceOptions);
            const auto& actual=caches[si];
            check(actual.offsets()==expected.offsets() && actual.cellIndices()==expected.cellIndices()
                && actual.weights()==expected.weights(),"current coupling must preserve order, hosts and periodic/cut policy");
            fallback+=expected.lastBuildFallbackCount();particles+=expected.lastBuildParticleCount();
            for(std::size_t p=0;p<actual.numParticles();++p) {
                double sum=0.;for(int j=actual.offsets()[p];j<actual.offsets()[p+1];++j)sum+=actual.weights()[j];
                near(sum,1.,"current-position partition of unity");
            }
        }
        check(counts.fallback==fallback && counts.particles==particles,"aggregate coupling counts");
        check(options.hostCellLabels==&ignoredHosts,"rebuild must not publish its scratch pointer in caller options");
        return counts;
    };
    verify(false);const auto first=caches[0].weights();verify(true);
    check(caches[0].weights()==first && species[0].particle[0].rngState==123
        && species[0].particle[0].position.x==.03,"parallel position copy preserves weights, particles and RNG");
    CouplingWeightsCache diagnostic;PicCurrentWeights::buildDiagnostic(mesh,species[1],hp,options,diagnostic);
    check(caches[0].weights()==first && caches.size()==3,"diagnostics must not replace physical caches");
    std::reverse(species[0].particle.begin(),species[0].particle.end());species[1].particle.clear();verify(true);
    check(caches[1].offsets()==std::vector<int>({0}),"empty species clears its previous entries");
    weights.clear();check(caches.empty(),"explicit reset invalidates borrowed caches");
    DictionaryMesh::RefineBox box;box.min={0.,.2};box.max={.3,.6};box.level=1;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()>16,"current coupling contract refines real topology");
    alpha.assign(mesh.cellIndexNum(),1.);
    for(auto& sp:species)for(auto& p:sp.particle)p.inCellL=mesh.findCellLabel(p.position);
    verify(true);
    alpha.assign(mesh.cellIndexNum(),0.);BoundaryLogs logs;
    check(verify(true).fallback==300,"zero-fluid support uses original host fallback and reports every particle");
    species.clear();const auto empty=verify(false);
    check(caches.empty() && empty.fallback==0 && empty.particles==0,"empty population resets caches and counts");
}
void couplingPolicy()
{
    using Policy=PicCurrentWeights::Policy;using Mode=Policy::Mode;
    static_assert(!std::is_copy_constructible_v<Policy> && !std::is_move_constructible_v<Policy>);
    static_assert(std::is_same_v<decltype(std::declval<Policy&>().halfWidths()),const std::vector<double>&>);
    for(auto mode:{Mode::Explicit,Mode::Implicit})for(int surfaces:{1,2}) {
        DictionaryMesh input;input.Nx=8;input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);PicCurrentWeights current(mesh);
        const auto embeddedMode=mode==Mode::Explicit?PicEmbeddedGeometry::Mode::Explicit:PicEmbeddedGeometry::Mode::Implicit;
        PicEmbeddedGeometry embedded(embeddedMode);
        DepositionConfig config;config.compositeBand=false;config.momentExactCut=true;
        config.debugParticleId=7;config.dumpWeightsOnce=false;config.dumpWeightsPath="coupling_snapshot.csv";
        Policy snapshot(mesh,config,mode,true,false);
        config.debugParticleId=-1;config.dumpWeightsOnce=true;config.dumpWeightsPath.clear();
        config.compositeBand=true;config.momentExactCut=false;
        const auto captured=snapshot.physicalOptions(embedded),diagnostic=snapshot.diagnosticOptions(embedded);
        check(captured.periodicX && !captured.periodicY && !captured.enableCompositeBand && captured.enableMomentExactCut,
            "coupling policy owns configured values after caller mutation");
        check(captured.debugParticleId==(mode==Mode::Explicit?7:std::numeric_limits<std::size_t>::max())
            && captured.dumpWeightsPath==(mode==Mode::Explicit?"coupling_snapshot.csv":""),"only explicit physical coupling carries debug selection");
        check(diagnostic.debugParticleId==std::numeric_limits<std::size_t>::max() && diagnostic.dumpWeightsPath.empty()
            && diagnostic.dumpWeightsOnce==(mode==Mode::Implicit),"diagnostic debug defaults and isolation");
        check(captured.levelSetI==(mode==Mode::Explicit?-1:0) && !captured.levelSetVert && !captured.alphaByCellIndex,
            "uninitialized geometry retains distinct physical defaults");
        config.momentExactCut=true;Policy policy(mesh,config,mode,false,true);
        const auto& widths=policy.halfWidths();const auto* widthAddress=widths.data();
        check(widths==std::vector<double>{.125},"initial minimum leaf support");
        DictionaryGeometry dict;dict.levelSetName={"left"};dict.innerPoints={{.1,.5}};dict.moving={false};
        dict.linesOnInterface={{{0,false}}};dict.phaseName={"left"};dict.levelSetOperation={{LevelSetOperation::non}};
        dict.levelSetIndexA={{{0,false}}};dict.levelSetIndexB={{}};dict.boundaryLines={{{0,false}}};
        if(surfaces==2) {
            dict.levelSetName.push_back("right");dict.innerPoints.push_back({.9,.5});dict.moving.push_back(false);
            dict.linesOnInterface.push_back({{1,false}});dict.phaseName.push_back("right");
            dict.levelSetOperation.push_back({LevelSetOperation::non});dict.levelSetIndexA.push_back({{1,false}});
            dict.levelSetIndexB.push_back({});dict.boundaryLines.push_back({{1,false}});
        }
        for(int k=0;k<surfaces;++k) {
            Line line;line.points[0]={k?.7:.3,0.};line.points[1]={k?.7:.3,1.};dict.lines.push_back(line);
            DictionaryGeometry::BoundaryEntity b;b.id=k+1;b.kind=GeometricBoundary::fixed;b.levelSetIndex=k;b.phase=k;b.lineIds={k};
            dict.boundaries.push_back(b);dict.boundaryById[k+1]=b;
        }
        ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh);FaceSegmentGraph graph;
        embedded.resolveSelection(dict,surfaces);
        auto rebuildGeometry=[&] {
            for(int k=0;k<surfaces;++k) {
                geometry.levelSets[k].initialize(mesh,0.);geometry.LS_dx(k).initialize(mesh,0.);geometry.LS_dy(k).initialize(mesh,0.);
                for(int vi=0;vi<mesh.vertIndexNum();++vi) {
                    const double x=mesh.verts(mesh.vertIndexToLabel[vi]).coord().x;
                    geometry.levelSets[k].value[vi]=k?.7-x:x-.3;geometry.LS_dx(k).value[vi]=k?-1.:1.;
                }
            }
            graph.build(mesh);embedded.rebuildComposite(mesh,&geometry);embedded.rebuildMetrics(mesh,&geometry,graph,{1,2});
        };
        rebuildGeometry();check(embedded.ready() && embedded.metricsReady() && embedded.usesComposite()==(surfaces==2),"real selected cut geometry");
        std::vector<ParticleType> species;species.emplace_back("charged",1.,-1.);species.emplace_back("empty",1.,0.);
        for(int j=0;j<300;++j) {
            Particle p;p.position={.32+.035*(j%11),j%2?.97:.03};p.w=2.;p.rngState=123+j;
            species[0].particle.push_back(p);
        }
        const auto equal=[](const CouplingWeightsCache& a,const CouplingWeightsCache& b) {
            check(a.offsets()==b.offsets() && a.cellIndices()==b.cellIndices() && a.weights()==b.weights(),
                "prepared coupling matches independent kernel inputs in CSR order and values");
        };
        auto verify=[&](double halfWidth) {
            std::vector<Vector2D> positions;std::vector<int> hosts;
            for(auto& p:species[0].particle){p.inCellL=mesh.findCellLabel(p.position);positions.push_back(p.position);hosts.push_back(p.inCellL);}
            // Independent inputs for the old physical paths, including the
            // implicit composite level-set index and disabled hp bookkeeping.
            CouplingWeightsCache::BuildOptions expected;
            expected.periodicY=true;expected.enableMomentExactCut=true;
            expected.levelSetI=mode==Mode::Explicit && surfaces==2?-1:0;
            expected.levelSetVert=embedded.levelSet();expected.alphaByCellIndex=&embedded.alphaByCell();
            expected.cutByAlpha=mode==Mode::Implicit || surfaces==2;expected.hostCellLabels=&hosts;
            if(mode==Mode::Implicit)expected.recordHpStats=expected.recordPerParticleHp=expected.recordPerParticleCompositeBand=false;
            CouplingWeightsCache reference;reference.build(mesh,positions,{halfWidth},expected);
            const auto counts=current.rebuild(species,policy,embedded,true);
            equal(current.caches()[0],reference);
            check(counts.particles==300 && counts.fallback==reference.lastBuildFallbackCount() && current.caches()[1].numParticles()==0,
                "policy rebuild keeps species and fallback bookkeeping");
            check(current.caches()[0].hasHpStats()==(mode==Mode::Explicit)
                && current.caches()[0].effectiveHp().size()==(mode==Mode::Explicit?300:0),"physical-only implicit bookkeeping policy");
            auto diagExpected=expected;diagExpected.levelSetI=surfaces==2?-1:0;diagExpected.cutByAlpha=surfaces==2;
            diagExpected.recordHpStats=diagExpected.recordPerParticleHp=diagExpected.recordPerParticleCompositeBand=true;
            CouplingWeightsCache diagReference,actual;diagReference.build(mesh,positions,{halfWidth},diagExpected);
            const auto saved=current.caches()[0].weights();PicCurrentWeights::buildDiagnostic(policy,embedded,species[0],actual,true);
            equal(actual,diagReference);check(actual.hasHpStats() && current.caches()[0].weights()==saved,"diagnostic build is independent of physical cache");
            for(std::size_t j=0;j<actual.numParticles();++j) {
                double sum=0.;for(int k=actual.offsets()[j];k<actual.offsets()[j+1];++k)sum+=actual.weights()[k];
                near(sum,1.,"cut periodic partition of unity");
            }
            check(species[0].particle[0].rngState==123 && species[0].particle[0].position.x==.32,"coupling does not mutate particles or RNG");
        };
        verify(.125);
        const auto gated=policy.physicalOptions(embedded,false);
        check(mode==Mode::Explicit?gated.alphaByCellIndex!=nullptr:!gated.alphaByCellIndex && !gated.levelSetVert,
            "implicit cut gate retains its original physical semantics");
        DictionaryMesh::RefineBox box;box.min={.25,0.};box.max={.5,.5};box.level=1;mesh.applyRefineBoxes({box});
        check(mesh.cellIndexNum()>16 && widths[0]==.125,"topology does not silently refresh support before its caller phase");
        check(policy.diagnosticOptions(embedded,true).alphaByCellIndex==nullptr,"diagnostic excludes stale alpha after topology change");
        if(mode==Mode::Explicit)check(policy.physicalOptions(embedded).alphaByCellIndex==nullptr,"explicit physical alpha size gate");
        policy.refreshAfterMeshChange();check(widths.data()==widthAddress && widths[0]==.0625,"AMR refresh preserves borrowed support vector");
        rebuildGeometry();verify(.0625);
        QuadtreeMesh other(input);Policy wrong(other,config,mode,false,true);const auto saved=current.caches()[0].weights();
        bool rejected=false;try{current.rebuild(species,wrong,embedded,true);}catch(const std::invalid_argument&){rejected=true;}
        check(rejected && current.caches()[0].weights()==saved,"mismatched mesh policy rejected before replacing physical cache");
    }
}
void nestedCouplingChunks()
{
    DictionaryMesh input;input.Nx=8;input.Ny=4;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.5,.75};box.level=1;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()>32,"nested coupling uses adaptive topology");
    std::vector<Vector2D> positions;for(int j=0;j<300;++j)positions.push_back({.01+.049*(j%20),.03+.067*(j%15)});
    std::vector<double> alpha(mesh.cellIndexNum(),0.);
    for(bool fallback:{false,true}) {
        CouplingWeightsCache::BuildOptions options;options.periodicX=true;options.levelSetI=-1;
        options.cutByAlpha=fallback;options.alphaByCellIndex=fallback?&alpha:nullptr;
        CouplingWeightsCache reference,nested;reference.enableRawWeights(true);nested.enableRawWeights(true);
        reference.build(mesh,positions,{.125},options);std::exception_ptr error;
#ifdef _OPENMP
        const int oldLevels=omp_get_max_active_levels();omp_set_max_active_levels(1);
        #pragma omp parallel num_threads(2)
        {
            #pragma omp single
            {
                try { nested.build(mesh,positions,{.125},options); }
                catch(...) { error=std::current_exception(); }
            }
        }
        omp_set_max_active_levels(oldLevels);
#else
        nested.build(mesh,positions,{.125},options);
#endif
        if(error)std::rethrow_exception(error);
        check(nested.numParticles()==300 && nested.offsets()==reference.offsets()
            && nested.cellIndices()==reference.cellIndices() && nested.weights()==reference.weights()
            && nested.rawWeights()==reference.rawWeights(),"serialized nested coupling processes all logical chunks with exact CSR and raw weights");
        check(nested.effectiveHp()==reference.effectiveHp() && nested.inCompositeBand()==reference.inCompositeBand()
            && nested.lastHpStats().count==300 && nested.lastHpStats().min==reference.lastHpStats().min
            && nested.lastHpStats().max==reference.lastHpStats().max && nested.lastHpStats().mean==reference.lastHpStats().mean,
            "nested chunk scheduling preserves adaptive support bookkeeping");
        check(nested.lastBuildFallbackCount()==(fallback?300:0) && nested.lastBuildFallbackCount()==reference.lastBuildFallbackCount(),
            "nested chunk scheduling retains every zero-weight fallback");
        for(int j=0;j<300;++j) {
            double sum=0.;for(int k=nested.offsets()[j];k<nested.offsets()[j+1];++k)sum+=nested.weights()[k];
            near(sum,1.,"nested coupling retains every particle's normalized stencil");
        }
    }
}
void orbitConfiguration()
{
    using Stage=ImplicitPicOrbit;
    DictionaryMesh input;input.Nx=8;input.Ny=2;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);PicEmbeddedGeometry empty(PicEmbeddedGeometry::Mode::Implicit);
    UnitsScales scales;scales.L0_m=2.;const auto units=UnitsContext::MakeSI(scales);
    PicConfig probe;ImplicitPicOverrides overrides;
    probe.segmenterMaxSegmentsPerParticle=9;probe.segmenterMaxLevel=3;probe.orbitAvgFused=true;
    probe.orbitAvgPatch.enabled=true;probe.orbitAvgAdaptive.sheathPad_m=.08;
    overrides.orbitAvgFused=false;overrides.orbitAvgSpeciesPar=true;
    overrides.segmenterMaxSegments=-5;overrides.segmenterMaxLevel=-2;
    overrides.segmenterTrace=true;overrides.segmenterIgnoreHostHints=true;
    Stage snapshot(mesh,0,probe,overrides,&units,true,false);
    probe.orbitAvgPatch.enabled=false;probe.orbitAvgAdaptive.sheathPad_m=-1.;overrides={};
    const auto stored=snapshot.optionsFor(empty,false,true);
    check(stored.averaging.useFusedPath && stored.averaging.segmentPatchEnabled && stored.parallelSpecies,
        "orbit configuration owns snapshot and patch forces fused despite explicit false override");
    near(stored.averaging.adaptiveCrossSplitSheathPad,.04,"orbit sheath width converts SI to code once");
    check(stored.segmenter.maxSegmentsPerParticle==16 && stored.segmenter.segmentCellMaxLevel==-2
        && stored.segmenter.trace && stored.segmenter.ignoreHostHints,"orbit override priority and original segment floor");
    check(stored.segmenter.periodicX && stored.weights.periodicX && !stored.segmenter.periodicY && !stored.weights.periodicY,
        "predictor and quadrature share configured periodicity");
    check(!stored.segmenter.enableEmbeddedIntercept && !stored.segmenter.embeddedLevelSetVert
        && !stored.weights.alphaByCellIndex && stored.weights.levelSetI==0 && !stored.weights.cutByAlpha,
        "disabled cut gate leaves original geometry defaults despite requested interception");
    Stage fallback(mesh,0,probe,overrides,nullptr,false,true);
    const auto defaulted=fallback.optionsFor(empty,false,false);
    check(defaulted.segmenter.maxSegmentsPerParticle==16 && defaulted.segmenter.segmentCellMaxLevel==3
        && defaulted.averaging.useFusedPath && !defaulted.parallelSpecies,"absent overrides preserve configured defaults");
    near(defaulted.averaging.adaptiveCrossSplitSheathPad,0.,"negative sheath pad retains original clamp");
    probe.orbitAvgAdaptive.sheathPad_m=.08;
    Stage noUnits(mesh,0,probe,overrides,nullptr,false,false);
    near(noUnits.optionsFor(empty,false,false).averaging.adaptiveCrossSplitSheathPad,.08,"absent unit context preserves raw pad");

    for(int surfaces:{1,2})for(int mode=0;mode<4;++mode) {
        DictionaryGeometry dict;dict.levelSetName={"left"};dict.innerPoints={{.1,.5}};dict.moving={false};
        dict.linesOnInterface={{{0,false}}};dict.phaseName={"left"};dict.levelSetOperation={{LevelSetOperation::non}};
        dict.levelSetIndexA={{{0,false}}};dict.levelSetIndexB={{}};dict.boundaryLines={{{0,false}}};
        if(surfaces==2) {
            dict.levelSetName.push_back("right");dict.innerPoints.push_back({.9,.5});dict.moving.push_back(false);
            dict.linesOnInterface.push_back({{1,false}});dict.phaseName.push_back("right");
            dict.levelSetOperation.push_back({LevelSetOperation::non});dict.levelSetIndexA.push_back({{1,false}});
            dict.levelSetIndexB.push_back({});dict.boundaryLines.push_back({{1,false}});
        }
        for(int k=0;k<surfaces;++k) {
            Line line;line.points[0]={k?.7:.3,0.};line.points[1]={k?.7:.3,1.};dict.lines.push_back(line);
            DictionaryGeometry::BoundaryEntity b;b.id=k+1;b.kind=GeometricBoundary::fixed;b.levelSetIndex=k;b.phase=k;b.lineIds={k};
            dict.boundaries.push_back(b);dict.boundaryById[k+1]=b;
        }
        QuadtreeMesh local(input);ImmersedGeometry<QuadtreeMesh> geometry(dict,local);FaceSegmentGraph graph;
        PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Implicit);embedded.resolveSelection(dict,surfaces);
        auto refresh=[&] {
            for(int k=0;k<surfaces;++k) {
                geometry.levelSets[k].initialize(local,0.);geometry.LS_dx(k).initialize(local,0.);geometry.LS_dy(k).initialize(local,0.);
                for(int vi=0;vi<local.vertIndexNum();++vi) {
                    const double x=local.verts(local.vertIndexToLabel[vi]).coord().x;
                    geometry.levelSets[k].value[vi]=k?.7-x:x-.3;geometry.LS_dx(k).value[vi]=k?-1.:1.;
                }
            }
            graph.build(local);embedded.rebuildComposite(local,&geometry);embedded.rebuildMetrics(local,&geometry,graph,{1,2});
        };
        refresh();check(embedded.ready() && embedded.metricsReady(),"orbit uses real cut geometry");
        std::vector<ParticleType> species;species.emplace_back("negative",1.,-1.);species.emplace_back("positive",1.,1.);species.emplace_back("empty",1.,0.);
        for(int si=0;si<2;++si)for(int j=0;j<96;++j) {
            Particle p;p.position={.34+.03*(j%11),j%2?.97:.03};p.velocity={si?-.7:.7,.5};p.rngState=123+j;
            species[si].set_kappa_s(.5);species[si].particle.push_back(p);
        }
        PicConfig config;config.orbitAvgFused=mode==1;config.orbitAvgPatch.enabled=mode==2;
        config.orbitAvgAdaptive.enabled=mode==3;config.orbitAvgAdaptive.tol=.002;config.orbitAvgAdaptive.maxRefineDepth=4;
        config.orbitAvgAdaptive.cheapRefine=true;config.orbitAvgAdaptive.crossSplitAtMaxDepth=true;
        config.orbitAvgAdaptive.crossSplitMaxIters=7;config.orbitAvgAdaptive.crossSplitMaxSplits=2;config.orbitAvgAdaptive.sheathPad_m=.08;
        config.orbitAvgPatch.tol=.004;config.orbitAvgPatch.maxRefineDepth=3;
        config.segmenterMaxSegmentsPerParticle=512;config.segmenterMaxLevel=-1;config.deposition.momentExactCut=true;
        ImplicitPicOverrides controls;controls.orbitAvgSpeciesPar=true;controls.segmenterIgnoreHostHints=true;
        Stage actual(local,3,config,controls,&units,false,true),reference(local,3);
        Stage::Options expected;
        expected.averaging.useFusedPath=mode==1 || mode==2;expected.averaging.segmentPatchEnabled=mode==2;
        expected.averaging.adaptiveEnabled=mode==3;expected.averaging.adaptiveTol=.002;expected.averaging.adaptiveMaxRefineDepth=4;
        expected.averaging.adaptiveCheapRefine=true;expected.averaging.adaptiveCrossSplitAtMaxDepth=true;
        expected.averaging.adaptiveCrossSplitMaxIters=7;expected.averaging.adaptiveCrossSplitMaxSplits=2;expected.averaging.adaptiveCrossSplitSheathPad=.04;
        expected.averaging.segmentPatchTol=.004;expected.averaging.segmentPatchMaxRefineDepth=3;
        expected.parallelSpecies=true;expected.segmenter.ignoreHostHints=true;
        expected.segmenter.periodicY=expected.weights.periodicY=true;expected.weights.enableMomentExactCut=true;
        config={};controls={}; // Subsequent caller changes cannot alter this run.
        CellField<Vector2D> field(local,{.1,.2});PicCurrentWeights current(local);
        auto compare=[&](double width,bool cut,bool intercept) {
            auto raw=expected;
            if(cut) {
                raw.segmenter.enableEmbeddedIntercept=intercept;raw.segmenter.embeddedLevelSetVert=embedded.levelSet();
                raw.segmenter.alphaByCellIndex=&embedded.alphaByCell();raw.segmenter.embeddedInsideIsNegative=embedded.insideIsNegative();
                raw.weights.levelSetI=embedded.levelSetIndex();raw.weights.levelSetVert=embedded.levelSet();
                raw.weights.alphaByCellIndex=&embedded.alphaByCell();raw.weights.insideIsNegative=embedded.insideIsNegative();raw.weights.cutByAlpha=true;
            }
            const auto fresh=actual.optionsFor(embedded,cut,intercept);
            check(fresh.segmenter.embeddedLevelSetVert==raw.segmenter.embeddedLevelSetVert
                && fresh.weights.levelSetI==raw.weights.levelSetI && fresh.weights.alphaByCellIndex==raw.weights.alphaByCellIndex
                && fresh.weights.recordHpStats && fresh.weights.recordPerParticleHp,
                "orbit keeps its own cut indices, fresh geometry and bookkeeping rather than current-weight policy");
            for(auto& sp:species)for(auto& p:sp.particle)p.inCellL=local.findCellLabel(p.position);
            current.rebuild(species,{width},raw.weights);
            actual.gather(species,current.caches(),field,1);reference.gather(species,current.caches(),field,1);
            std::vector<Stage::SpeciesPerf> right;
            try { right=reference.build(species,{width},.15,raw,{1,false,false}); }
            catch(const std::exception& error) { throw std::runtime_error(std::string("original raw orbit path: ")+error.what()); }
            const auto left=actual.build(species,{width},.15,embedded,cut,intercept,{1,false,false});
            for(int si=0;si<3;++si) {
                const auto& a=actual.averagedCaches()[si];const auto& b=reference.averagedCaches()[si];
                check(a.numParticles()==species[si].particle.size() && a.offsets()==b.offsets() && a.cellIndices()==b.cellIndices() && a.weights()==b.weights(),"configured orbit exactly matches original option path");
                check(left[si].predictor_segments_total==right[si].predictor_segments_total
                    && left[si].orbitAvg_sampleCount==right[si].orbitAvg_sampleCount
                    && left[si].orbitAvg_adaptAcceptCount==right[si].orbitAvg_adaptAcceptCount
                    && left[si].orbitAvg_adaptRefineCount==right[si].orbitAvg_adaptRefineCount
                    && left[si].orbitAvg_fused_mode==right[si].orbitAvg_fused_mode,"configured orbit preserves integration path and event counts");
                for(std::size_t j=0;j<a.numParticles();++j) {
                    double sum=0.;for(int k=a.offsets()[j];k<a.offsets()[j+1];++k)sum+=a.weights()[k];
                    near(sum,1.,"configured moving orbit conserves per-particle weight");
                }
            }
            check(species[0].particle[0].position.x==.34 && species[0].particle[0].velocity.x==.7
                && species[0].particle[0].rngState==123,"configured prediction leaves physical particles and RNG unchanged");
        };
        compare(.125,true,true);compare(.125,false,true);compare(.125,true,false);
        DictionaryMesh::RefineBox box;box.min={.25,0.};box.max={.5,.5};box.level=1;local.applyRefineBoxes({box});
        check(local.cellIndexNum()>16,"configured orbit exercises real mesh change");
        actual.resetScratch();reference.resetScratch();refresh();field.initialize(local,{.1,.2});compare(.0625,true,true);
    }
}
void orbitStage()
{
    static_assert(!std::is_copy_constructible_v<ImplicitPicOrbit>);
    static_assert(!std::is_move_constructible_v<ImplicitPicOrbit>);
    static_assert(std::is_same_v<decltype(std::declval<ImplicitPicOrbit&>().averagedCaches()[0]),const CouplingWeightsCache&>);
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);CellField<Vector2D> field(mesh,{0.,0.});
    std::vector<ParticleType> species;
    species.emplace_back("one",1.,1.);species.emplace_back("two",1.,-1.);species.emplace_back("empty",1.,0.);
    for(int si=0;si<2;++si) {
        species[si].set_kappa_s(.5);species[si].particle.resize(1);
        auto& p=species[si].particle[0];p.position={.25+.5*si,.375};p.velocity={0.,0.};p.rngState=123+si;
        p.inCellL=mesh.findCellLabel(p.position);
    }
    ImplicitPicOrbit stage(mesh,species.size());ImplicitPicOrbit::Options options;
    const std::vector<double> hp{.125};std::vector<CouplingWeightsCache> current(3);
    const auto prepare=[&] {
        for(std::size_t si=0;si<species.size();++si) {
            std::vector<Vector2D> positions;for(const auto& p:species[si].particle)positions.push_back(p.position);
            current[si].build(mesh,positions,hp,options.weights);
        }
        ImplicitPicOrbit::validateCurrentCaches(species,current);
        stage.gather(species,current,field,0);
    };
    const auto normalized=[&] {
        for(int si=0;si<2;++si) {
            const auto& cache=stage.averagedCaches()[si];
            check(cache.numParticles()==1 && cache.nnz()>0,"orbit stage publishes owned nonempty CSR");
            double total=0.;for(double w:cache.weights()){check(w>=0.,"negative orbit weight");total+=w;}
            near(total,1.,"orbit weights preserve particle number");
        }
        check(stage.averagedCaches()[2].numParticles()==0 && stage.averagedCaches()[2].nnz()==0,
              "empty species has no stale orbital entries");
    };
    for(int mode=0;mode<4;++mode) {
        options.averaging.adaptiveEnabled=mode==3;
        options.averaging.useFusedPath=mode==1 || mode==2;
        options.averaging.segmentPatchEnabled=mode==2;
        options.parallelSpecies=mode==3;
        prepare();const auto performance=stage.build(species,hp,.1,options,{});
        normalized();check(performance.size()==3 && performance[0].predictor_segments_total>0,"per-species orbit report");
        for(int si=0;si<2;++si) {
            const auto& averaged=stage.averagedCaches()[si];const auto& atRest=current[si];
            check(averaged.cellIndices()==atRest.cellIndices(),"stationary orbit preserves instantaneous cell support");
            for(std::size_t j=0;j<atRest.nnz();++j)near(averaged.weights()[j],atRest.weights()[j],"stationary orbit quadrature weights");
            check(species[si].particle[0].velocity.x==0. && species[si].particle[0].rngState==123+si,
                  "predictor must not push particles or consume their RNG");
        }
    }
    options.averaging.segmentPatchEnabled=false;options.parallelSpecies=false;
    options.weights.periodicX=options.segmenter.periodicX=true;
    auto& particle=species[0].particle[0];particle.position.x=.95;particle.velocity.x=.9;
    particle.inCellL=mesh.findCellLabel(particle.position);
    field.initialize(mesh,{2.,0.});prepare();const auto oldCurrent=current[0].weights();
    stage.build(species,hp,.2,options,{});normalized();
    check(particle.position.x==.95 && particle.velocity.x==.9 && particle.rngState==123,
          "periodic kick/drift prediction leaves current particle layer unchanged");
    check(current[0].weights()==oldCurrent,"orbit averaging must not overwrite current-position cache");
    bool left=false,right=false;
    for(int ci:stage.averagedCaches()[0].cellIndices()) {
        const double x=mesh.cells(mesh.cellIndexToLabel[ci]).centerCoord().x;
        left=left || x<.25;right=right || x>.75;
    }
    check(left && right,"periodic orbit crosses both sides of the box");

    // A path across 40 fine cells forces the original segment-cap recovery.
    {
        DictionaryMesh retryInput;retryInput.Nx=40;retryInput.Ny=1;retryInput.maxCellNum=160;retryInput.maxInitialIterations=0;
        QuadtreeMesh retryMesh(retryInput);CellField<Vector2D> zero(retryMesh,{0.,0.});
        std::vector<ParticleType> retrySpecies;retrySpecies.emplace_back("retry",1.,1.);
        retrySpecies[0].particle.resize(1);auto& p=retrySpecies[0].particle[0];
        p.position={.0125,.5};p.velocity={.975,0.};p.inCellL=retryMesh.findCellLabel(p.position);
        ImplicitPicOrbit retry(retryMesh,1);ImplicitPicOrbit::Options retryOptions;
        retryOptions.segmenter.maxSegmentsPerParticle=16;
        std::vector<CouplingWeightsCache> retryCurrent(1);
        retryCurrent[0].build(retryMesh,{p.position},{.0125},retryOptions.weights);
        retry.gather(retrySpecies,retryCurrent,zero,0);
        BoundaryLogs logs;const auto perf=retry.build(retrySpecies,{.0125},1.,retryOptions,{0,true,false});
        check(logs.output.str().find("retry baseline")!=std::string::npos && perf[0].predictor_segments_total>16,
              "segment cap retries with the original enlarged budget");
    }
    normalized();particle.velocity.x=0.;
    stage.resetScratch();
    check(stage.averagedCaches()[0].nnz()==0,"AMR scratch reset invalidates averaged cache");
    DictionaryMesh::RefineBox box;box.min={.8,.2};box.max={1.,.5};box.level=1;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()>16,"orbit contract refines real topology");
    field.initialize(mesh,{0.,0.});
    for(auto& sp:species)for(auto& p:sp.particle)p.inCellL=mesh.findCellLabel(p.position);
    prepare();stage.build(species,hp,.1,options,{});normalized();
    const auto rejects=[&](auto action,const char* why) {
        bool caught=false;try{action();}catch(const std::exception&){caught=true;}check(caught,why);
    };
    const auto before=stage.averagedCaches()[0].weights();auto bad=options;bad.weights.periodicX=false;
    rejects([&]{stage.build(species,hp,.1,bad,{});},"mismatched predictor/weight periodicity rejected");
    check(stage.averagedCaches()[0].weights()==before,"option rejection precedes cache mutation");
    rejects([&]{ImplicitPicOrbit::validateCurrentCaches(species,{});},"missing current species cache rejected");
    auto badCache=current;badCache[0].assignCSR(0,{0},{},{});
    rejects([&]{ImplicitPicOrbit::validateCurrentCaches(species,badCache);},"current cache particle identity mismatch rejected");
    ImplicitPicOrbit ungathered(mesh,species.size());bad=options;bad.parallelSpecies=true;
    rejects([&]{ungathered.build(species,hp,.1,bad,{});},"species worker errors propagate after parallel join");

    // Large populations activate inner worker chunks. A serialized nested
    // team must still process all chunks when species use the outer team.
#ifdef _OPENMP
    const int oldLevels=omp_get_max_active_levels();omp_set_max_active_levels(1);
#endif
    auto many=species;
    for(int si=0;si<2;++si)many[si].particle.assign(96,species[si].particle[0]);
    std::vector<CouplingWeightsCache> manyCurrent(3);
    for(std::size_t si=0;si<many.size();++si) {
        std::vector<Vector2D> positions;for(const auto& p:many[si].particle)positions.push_back(p.position);
        manyCurrent[si].build(mesh,positions,hp,options.weights);
    }
    ImplicitPicOrbit serial(mesh,3),parallel(mesh,3);
    for(int mode=0;mode<4;++mode) {
        options.averaging.adaptiveEnabled=mode==3;
        options.averaging.useFusedPath=mode==1 || mode==2;
        options.averaging.segmentPatchEnabled=mode==2;
        options.parallelSpecies=false;serial.gather(many,manyCurrent,field,0);
        const auto reference=serial.build(many,hp,.1,options,{});
        options.parallelSpecies=true;parallel.gather(many,manyCurrent,field,0);
        const auto actual=parallel.build(many,hp,.1,options,{});
        for(int si=0;si<3;++si) {
            const auto& a=serial.averagedCaches()[si];const auto& b=parallel.averagedCaches()[si];
            check(a.numParticles()==many[si].particle.size() && a.offsets()==b.offsets()
                && a.cellIndices()==b.cellIndices() && a.weights()==b.weights(),
                "nested-team orbit processing retains all particles and exact CSR weights");
            check(reference[si].predictor_segments_total==actual[si].predictor_segments_total
                && reference[si].orbitAvg_sampleCount==actual[si].orbitAvg_sampleCount,
                "nested-team predictor and quadrature counts match full inner team");
        }
    }
#ifdef _OPENMP
    omp_set_max_active_levels(oldLevels);
#endif
}
void polarizationStage()
{
    static_assert(!std::is_copy_constructible_v<ImplicitPicPolarization>);
    static_assert(!std::is_move_constructible_v<ImplicitPicPolarization>);
    static_assert(std::is_same_v<decltype(std::declval<ImplicitPicPolarization&>().chargeAmount().value[0]),const double&>);
    static_assert(std::is_same_v<decltype(std::declval<ImplicitPicPolarization&>().cellCoefficients()[0]),const double&>);
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    std::vector<ParticleType> species;
    species.emplace_back("negative",2.,-2.);species.emplace_back("positive",4.,1.);
    species.emplace_back("neutral",1.,0.);species.emplace_back("zero_mass",0.,3.);
    species.emplace_back("empty",2.,5.);
    species[0].particle.resize(2);species[0].particle[0].w=2.;species[0].particle[1].w=0.;
    species[1].particle.resize(1);species[1].particle[0].w=3.;
    species[2].particle.resize(1);species[2].particle[0].w=-1.;
    species[3].particle.resize(1);species[3].particle[0].w=1.;
    std::vector<CouplingWeightsCache> caches(5);
    caches[0].assignCSR(2,{0,2,3},{0,1,1},{.25,.75,1.});
    caches[1].assignCSR(1,{0,2},{1,2},{.5,.5});
    caches[2].assignCSR(1,{0,1},{0},{1.});
    caches[3].assignCSR(1,{0,1},{2},{1.});caches[4].assignCSR(0,{0},{},{});
    const auto norm=NormalizationContext::MakeLegacy();
    ImplicitPicPolarization stage(mesh,species.size(),graph.segments().size());
    const std::vector<double> expectedQ{-1.,-9.5,4.5,0.};
    const std::vector<double> expectedMu{.04,.455,.015,0.};
    std::vector<double> coefficients;
    for(int threshold:{1000000,1}) {
        PicStepParticleBalance balance(0.);
        const auto stats=stage.deposit(species,{4.},caches,threshold,balance);
        check(stats.particles==5 && stats.entries==7,"deposition counts only nonempty species and their CSR entries");
        near(balance.current().stageMacroC,-6.,"macro charge uses positive per-particle weights or species/default fallback");
        near(stage.totalCharge(),-6.,"orbit deposition conserves uncorrected charge");
        for(int ci=0;ci<4;++ci)near(stage.chargeAmount().value[ci],expectedQ[ci],"analytical orbit-weight deposition");
        near(stage.numberAmounts()[0].value[0],.5,"fractional electron number amount");
        near(stage.numberAmounts()[0].value[1],5.5,"fallback electron weight");
        near(stage.numberAmounts()[2].value[0],1.,"neutral number deposit with omitted default weight");
        stage.assembleCells(species,norm,.2);
        for(int ci=0;ci<4;++ci) {
            near(stage.chargeDensity().value[ci],norm.RhoSIToTilde(expectedQ[ci]/.25),"density uses full SI cell volume");
            near(stage.susceptibilitySI()[ci],expectedMu[ci],"charge squared/mass susceptibility excludes neutral and zero-mass species");
        }
        stage.assembleFaces(graph);
        if(coefficients.empty())coefficients=stage.faceCoefficients();
        else check(coefficients==stage.faceCoefficients(),"OpenMP particle fill preserves deterministic cell/face coefficients");
        bool interior=false,boundary=false;
        for(std::size_t sid=0;sid<graph.segments().size();++sid) {
            const auto& face=graph.segments()[sid];
            const int k=mesh.cells(face.K).index();
            const double a=1.+expectedMu[k]/PhysicalConstants::epsilon0_Fpm;
            double expected=a;
            if(face.L>=0) {
                interior=true;const int l=mesh.cells(face.L).index();
                const double b=1.+expectedMu[l]/PhysicalConstants::epsilon0_Fpm;
                expected=2.*a*b/(a+b);
            } else boundary=true;
            check(std::abs(stage.faceCoefficients()[sid]-expected)<=5e-13*std::max(1.,std::abs(expected)),
                  "face coefficient uses harmonic interior and one-sided boundary values");
        }
        check(interior && boundary,"polarization fixture covers both face closures");
        const auto correction=stage.removeMeanCharge(norm);
        check(correction.applied,"Neumann compatibility correction applied");
        near(correction.totalQ,-6.,"correction reports original total charge");
        near(correction.meanRhoSI,-6.,"correction uses total full-cell volume");
        near(stage.totalCharge(),0.,"mean correction enforces charge compatibility");
        balance.finishStage(stage.totalCharge());
        near(balance.current().stageResidualC,6.,"stage diagnostic sees correction without rewriting macro charge");
        for(int ci=0;ci<4;++ci)near(stage.chargeAmount().value[ci],expectedQ[ci]+1.5,"uniform-volume neutrality correction");
        check(stage.faceCoefficients()==coefficients,"neutrality correction does not change susceptibility");
        stage.resetScratch();
    }
    const auto rejects=[&](auto action,const char* reason) {
        bool caught=false;try{action();}catch(const std::exception&){caught=true;}check(caught,reason);
    };
    PicStepParticleBalance untouched(0.);
    rejects([&]{stage.deposit(species,{4.},{},1,untouched);},"short cache array rejected before reset");
    near(stage.chargeAmount().value[0],.5,"invalid input preserves deposited values");
    near(untouched.current().stageMacroC,0.,"invalid input preserves charge diagnostic");
    rejects([&]{stage.assembleCells({},norm,.2);},"changed species cardinality rejected");
    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={.49,.49};box.level=1;
    mesh.applyRefineBoxes({box});graph.build(mesh);
    check(mesh.cellIndexNum()>4,"polarization contract changes real mesh topology");
    rejects([&]{stage.deposit(species,{4.},caches,1,untouched);},"stale polarization cannot index refined mesh");
    stage.resetDepositedFields();stage.resetCoefficients(graph.segments().size());stage.resetScratch();
    for(int ci=0;ci<mesh.cellIndexNum();++ci) {
        check(stage.chargeAmount().label(ci)==mesh.cellIndexToLabel[ci]
            && stage.numberAmounts()[0].label(ci)==mesh.cellIndexToLabel[ci],"AMR rebuild rebinds deposit labels");
        near(stage.chargeAmount().value[ci],0.,"AMR rebuild resets charge");
        near(stage.cellCoefficients()[ci],1.,"AMR rebuild starts with vacuum coefficient");
    }
    stage.deposit(species,{4.},caches,1,untouched);stage.assembleCells(species,norm,.2);stage.assembleFaces(graph);
    near(stage.totalCharge(),-6.,"deposit workspaces work after mesh reconstruction");
    check(stage.faceCoefficients().size()==graph.segments().size(),"refined face coefficient size");
}
void boundarySetup() {
    using namespace PicBoundarySetup;
    DictionaryGeometry geometry;geometry.domainMin={0.,0.};geometry.domainMax={1.,1.};
    const auto line=[&](Vector2D a,Vector2D b) {Line l;l.points[0]=a;l.points[1]=b;geometry.lines.push_back(l);};
    line({0.,0.},{0.,1.});line({1.,0.},{1.,1.});line({0.,0.},{1.,0.});line({0.,1.},{1.,1.});line({.5,0.},{.5,1.});
    const auto entity=[&](int id,const char* name,GeometricBoundary kind,std::vector<int> lines,bool ignored=false) {
        DictionaryGeometry::BoundaryEntity e;e.id=id;e.name=name;e.kind=kind;e.lineIds=std::move(lines);e.ignored=ignored;
        geometry.boundaries.push_back(e);geometry.boundaryById[id]=e;
    };
    entity(10,"not_a_side",GeometricBoundary::domain,{0});entity(11,"split_left",GeometricBoundary::domain,{0});
    entity(20,"right",GeometricBoundary::domain,{1});entity(30,"bottom",GeometricBoundary::domain,{2});
    entity(40,"top",GeometricBoundary::domain,{3});entity(50,"wall",GeometricBoundary::fixed,{-1,4,99});
    entity(51,"same_line",GeometricBoundary::moving,{4});entity(60,"ignored",GeometricBoundary::fixed,{0},true);
    entity(70,"ToP",GeometricBoundary::domain,{});entity(71,"left",GeometricBoundary::domain,{0,1});
    PicConfig cfg;cfg.domainBC.left.value=9.;cfg.embeddedBC.enabled=true;cfg.embeddedBC.g=12.;
    PicBoundarySide left;left.type=PicBCType::Dirichlet;left.value=2.;left.timeDirichlet.enabled=true;
    left.timeDirichlet.type="SINE";left.timeDirichlet.ac=3.;left.timeDirichlet.freq_hz=4.;
    PicBoundarySide same=left;same.value+=1e-13;same.timeDirichlet.type="sine";
    PicBoundarySide neumann;neumann.type=PicBCType::Neumann;neumann.value=5.;
    cfg.fieldBCById={{10,left},{11,same},{50,neumann},{51,neumann},{60,left},{999,left}};
    for(auto mode:{Mode::Explicit,Mode::Implicit}) {
        BoundaryLogs logs;const auto mapped=mapFields(cfg,geometry,true,mode);
        check(mapped.domain.left.value==same.value && mapped.domainIds[0]==10,"geometry side inference, last equal BC and first assigned ID");
        check(mapped.domain.left.timeDirichlet.type=="sine" && mapped.embedded.g==5.
            && mapped.embedded.a==0. && mapped.embedded.b==1.,"waveform case equality and common Neumann fallback");
        check(mapped.lineToBoundaryId==std::vector<int>({0,0,0,0,50}) && mapped.embeddedById.size()==2,"line first-owner, invalid lines and ignored boundaries");
        check((logs.errors.str().find("multiple boundaryIds")!=std::string::npos)==(mode==Mode::Explicit),"mode-specific duplicate-line warning");
        check((logs.errors.str().find("ignoring Field BC")!=std::string::npos)==(mode==Mode::Explicit),"mode-specific missing-ID warning");
        const auto absent=mapFields(cfg,geometry,false,mode);
        check(absent.domain.left.value==9. && absent.embedded.g==12. && absent.embeddedById.empty()
            && absent.lineToBoundaryId.empty() && absent.domainIds[0]==-1,"no geometry retains global defaults without mapping");
        PicConfig distinct=cfg;distinct.fieldBCById[3].bc.value=6.;
        check(mapFields(distinct,geometry,true,mode).embedded.g==12.,"different embedded IDs preserve global fallback");
        PicConfig fallback;fallback.fieldBCById={{70,left}};
        check(mapFields(fallback,geometry,true,mode).domain.top.value==2.,"legacy case-insensitive name fallback without lines");
        const auto reject=[&](PicConfig bad,const std::string& message) {
            bool threw=false;
            try {mapFields(bad,geometry,true,mode);} catch(const std::runtime_error& e) {
                threw=true;const std::string text=e.what();
                check(text.find(mode==Mode::Explicit?"[ExplicitPIC][BC] ERROR:":"[ImplicitPIC][BC] ERROR:")==0,"boundary error prefix");
                check(text.find(message)!=std::string::npos,"specific boundary error reason");
            }
            check(threw,"invalid boundary setup must fail");
        };
        auto conflict=cfg;conflict.fieldBCById[1].bc.timeDirichlet.ac=4.;reject(conflict,"conflict (ids 10 and 11)");
        conflict=cfg;conflict.fieldBCById.push_back({50,left});reject(conflict,"conflicting embedded BC definitions");
        conflict=cfg;conflict.fieldBCById[2].bc.type=PicBCType::Robin;reject(conflict,"BoundarySide has no (a,b)");
        conflict=cfg;conflict.fieldBCById={{71,left}};reject(conflict,"cannot infer domain side");
        // Disabled waveform parameters are ignored when comparing otherwise equal BCs.
        conflict=cfg;conflict.fieldBCById[0].bc.timeDirichlet.enabled=false;
        conflict.fieldBCById[1].bc.timeDirichlet.enabled=false;conflict.fieldBCById[1].bc.timeDirichlet.ac=999.;
        mapFields(conflict,geometry,true,mode);
    }
    BoundaryLogs logs;
    PicConfig periodic;PicBoundarySide wrap;wrap.type=PicBCType::Periodic;
    periodic.fieldBCById={{10,wrap},{20,wrap},{30,wrap},{40,wrap}};
    geometry.periodicMateByBoundary[20]=10; // One direction is sufficient in the original explicit gate.
    auto mapped=prepareExplicit(periodic,geometry,true);
    check(mapped.periodicX && !mapped.periodicY,"explicit requires geometry mates as well as field periodicity");
    check(logs.errors.str().find("Periodic wrap/ghost will remain disabled")!=std::string::npos,"unpaired periodic warning retained");
    geometry.periodicMateByBoundary[30]=40;
    check(prepareExplicit(periodic,geometry,true).periodicY,"paired y wrap enabled");
    check(!prepareExplicit(periodic,geometry,false).periodicX,"no geometry disables explicit periodicity");
    PicConfig particles;particles.see.enabled=false;particles.see.gamma=.1;
    PicConfig::BoundaryIdParticleBC off;off.id=50;off.seeSpecified=true;off.see.gamma=.7;
    particles.particleBCById={off};
    check(prepareExplicit(particles,geometry,true).particles.see.gamma==.7
        && prepareImplicitParticles(particles,geometry).see.gamma==.1,"explicit disabled template copy differs from implicit enabled-template selection");
    auto on=off;on.see.enabled=true;on.see.gamma=.8;particles.particleBCById.push_back(on);
    on.id=10;on.see.gamma=.9;particles.particleBCById.push_back(on);
    on.id=60;particles.particleBCById.push_back(on);on.id=999;on.model="reflect";particles.particleBCById.push_back(on);
    auto exp=prepareExplicit(particles,geometry,true);auto imp=prepareImplicitParticles(particles,geometry);
    check(exp.particles.see.gamma==.8 && imp.see.gamma==.8 && exp.particles.embeddedIds==std::unordered_set<int>{50}
        && imp.embeddedIds==exp.particles.embeddedIds,"first enabled SEE template and ignored/missing boundaries");
    check(exp.particles.see.cathodeSides==std::vector<std::string>{"left"} && imp.see.cathodeSides.empty(),"only explicit adds domain SEE side");
    particles.see.enabled=true;particles.see.gamma=.3;
    check(prepareExplicit(particles,geometry,true).particles.see.gamma==.3
        && prepareImplicitParticles(particles,geometry).see.gamma==.3,"enabled global SEE parameters retained");
    particles.fieldBCById.clear();particles.particleBCById[0].model="reflect";
    prepareExplicit(particles,geometry,true);
    check(logs.errors.str().find("is not implemented (only absorb)")!=std::string::npos,"unsupported explicit particle model warning");
    // Snapshots remain valid when their input containers are changed or released.
    particles.particleBCById.clear();geometry.lines.clear();geometry.boundaryById.clear();
    check(exp.fields.lineToBoundaryId.size()==5 && exp.particles.see.gamma==.8 && imp.embeddedIds.count(50)==1,"setup results own their data");
}
struct Totals { double weight=0.,px=0.,py=0.,pz=0.,energy=0.; };
Totals totals(const ParticleType& species)
{
    Totals out;
    for (const auto& p:species.particle)
    {
        out.weight+=p.w; out.px+=p.w*p.velocity.x; out.py+=p.w*p.velocity.y; out.pz+=p.w*p.vz;
        out.energy+=p.w*(p.velocity.mod2()+p.vz*p.vz);
    }
    return out;
}
void conserved(const Totals& a,const Totals& b)
{
    near(a.weight,b.weight,"macro weight conservation"); near(a.px,b.px,"momentum x");
    near(a.py,b.py,"momentum y"); near(a.pz,b.pz,"momentum z"); near(a.energy,b.energy,"equal-velocity kinetic energy");
}
void speciesOwnership(const QuadtreeMesh& mesh)
{
    static_assert(!std::is_copy_constructible_v<PicSpeciesState> && !std::is_move_constructible_v<PicSpeciesState>);
    static_assert(std::is_same_v<decltype(std::declval<const PicSpeciesState&>().species()),const std::vector<ParticleType>&>);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<PicSpeciesState&>().weights())>>);
    NormalizationScales scales;scales.L0_m=2.;scales.T0_s=.5;scales.Phi0_V=3.;scales.v0_mps=4.;
    std::vector<SpeciesConfig> config(4);
    for(auto& s:config){s.m=2.*PhysicalConstants::eV_J;s.Te_eV=8.;s.numParticles=7;}
    config[0].name="shared";config[0].q=1.;config[0].weight=2.;
    config[1].name="shared";config[1].q=-1e-31;config[1].weight=3.;config[1].numParticles=9;
    config[2].name="PhotoELECTRON";config[2].q=1.;config[2].weight=4.;config[2].numParticles=0;
    config[3].name="neutral";config[3].q=0.;config[3].weight=5.;config[3].Te_eV=0.;
    PicSpeciesState owner(config.size(),0);PicSpeciesState::Engine reference(12345);
    check(owner.rng()==reference && owner.species().empty() && owner.weights().empty(),"species owner uses original zero-seed fallback without sampling");
    auto* speciesIdentity=&owner.species();auto* rngIdentity=&owner.rng();const auto* weightIdentity=&owner.weights();const auto* indexIdentity=&owner.index();
    const auto exclude=[](const Vector2D& point,int){return point.x<.25;};
    std::size_t observed=0;
    owner.initialize(config,scales,mesh,false,exclude,[&](const SpeciesConfig& s,const ParticleType& p){
        check(&s==&config[observed] && owner.species().size()==observed && owner.weights().size()==observed,"initial diagnostic callback precedes species commit");
        check(p.name()==s.name && p.number==p.particle.size(),"initial diagnostic sees sampled species");++observed;
    });
    check(observed==4 && owner.species().size()==4 && owner.weights()==std::vector<double>({2.,3.,4.,5.}),"species order and configured weights owned together");
    for(std::size_t j=0;j<config.size();++j){ParticleType expected(config[j].name,config[j].m,config[j].q);
        initializePicSpecies(expected,config[j],scales,mesh,reference,false,exclude);const auto& actual=owner.species()[j];
        check(actual.number==expected.number && actual.kappa_s()==expected.kappa_s(),"owner preserves original species initialization");
        for(std::size_t k=0;k<expected.particle.size();++k){const auto& a=actual.particle[k];const auto& b=expected.particle[k];
            check(a.position.x==b.position.x && a.position.y==b.position.y && a.velocity.x==b.velocity.x && a.velocity.y==b.velocity.y && a.vz==b.vz && a.w==b.w && a.inCellL==b.inCellL && a.inBackgroundCellL==b.inBackgroundCellL,"owner preserves samples and rejected-candidate RNG across species");}
    }
    check(owner.rng()==reference,"owned master RNG matches original sequential calls");
    check(owner.ions().count("shared") && owner.electrons().count("shared") && owner.neutrals().count("shared") && owner.neutrals().count("neutral"),"configured classifications preserve overlapping sign and near-neutral rules");
    check(owner.index().empty(),"name index is not built before checkpoint restoration");owner.rebuildIndex();
    check(owner.index().size()==3 && owner.index().at("shared")==0 && owner.electronIndex()==2,"first duplicate wins while electron name precedes charge fallback");
    config[0].weight=99.;config[0].name="changed";
    check(owner.weights()[0]==2. && owner.ions().count("shared") && !owner.ions().count("changed"),"configured metadata survives later input edits");
    owner.species()[2]=ParticleType("restored-ion",1.,1.);owner.rebuildIndex();
    check(owner.index().count("PhotoELECTRON")==0 && owner.index().at("restored-ion")==2 && owner.electronIndex()==1,"restored runtime names replace index and electron lookup falls back to charge");
    check(owner.ions().count("PhotoELECTRON") && owner.weights()[2]==4.,"restored metadata does not rewrite original classification or default weights");
    const auto savedRng=owner.rng();bool duplicate=false;
    try{owner.initialize(config,scales,mesh,false);}catch(const std::logic_error&){duplicate=true;}
    check(duplicate && owner.rng()==savedRng && owner.species().size()==4,"duplicate initialization fails before draws or mutation");
    check(&owner.species()==speciesIdentity && &owner.rng()==rngIdentity && &owner.weights()==weightIdentity && &owner.index()==indexIdentity,"borrowed container and RNG identities remain stable");
    PicSpeciesState restored(config.size(),12345);int excluded=0;
    restored.initialize(config,scales,mesh,true,[&](const auto&,int){++excluded;return false;});
    check(excluded==0 && restored.rng()==PicSpeciesState::Engine(12345),"checkpoint initialization never samples or calls exclusion");
    for(std::size_t j=0;j<config.size();++j){ParticleType expected(config[j].name,config[j].m,config[j].q);PicSpeciesState::Engine unused(12345);
        initializePicSpecies(expected,config[j],scales,mesh,unused,true);
        check(restored.species()[j].particle.empty() && restored.species()[j].number==0 && restored.species()[j].kappa_s()==expected.kappa_s(),"checkpoint population retains prepared normalization without sampling");}
    std::ostringstream serialized;serialized<<savedRng;std::istringstream input(serialized.str());input>>restored.rng();
    check(restored.rng()==savedRng && owner.rng()==savedRng,"checkpoint IO can restore only the selected owner's master engine");
    (void)restored.rng()();check(owner.rng()==savedRng && restored.rng()!=owner.rng(),"simultaneous sessions do not share RNG");
    PicSpeciesState empty(0,9);bool premature=false;try{empty.rebuildIndex();}catch(const std::logic_error&){premature=true;}
    check(premature && empty.index().empty(),"index publication requires initialization");empty.initialize({},scales,mesh,true);empty.rebuildIndex();
    check(empty.species().empty() && empty.index().empty() && empty.electronIndex()==-1 && empty.rng()==PicSpeciesState::Engine(9),"empty checkpoint population has no implicit sampling or electron");
}
void initialization(const QuadtreeMesh& mesh)
{
    speciesOwnership(mesh);
    NormalizationScales scales;
    scales.L0_m=2.; scales.T0_s=.5; scales.Phi0_V=3.; scales.v0_mps=4.;
    SpeciesConfig config;
    config.name="initial"; config.m=2.*PhysicalConstants::eV_J; config.q=-4.*PhysicalConstants::eV_J;
    config.weight=7.; config.numParticles=17; config.drift={8.,-4.};
    std::mt19937_64 rng(24680), positionOnly(rng);
    ParticleType cold(config.name,config.m,config.q);
    initializePicSpecies(cold,config,scales,mesh,rng,false);
    check(cold.number==17 && cold.particle.size()==17,"cold initial count");
    near(cold.kappa_s(),-.375,"normalized initial acceleration coefficient");
    std::uniform_real_distribution<double> draw(0.,1.);
    for(const auto& p:cold.particle) {
        draw(positionOnly);draw(positionOnly);
        check(p.velocity.x==2. && p.velocity.y==-1. && p.vz==0. && p.w==7.,"cold drift and macro weight");
        check(p.inCellL==mesh.findCellLabel(p.position)
            && p.inBackgroundCellL==mesh.findBackgroundCellLabel(p.position),"initial mesh identity");
        check(p.position.x>1e-6 && p.position.x<1.-1e-6
            && p.position.y>1e-6 && p.position.y<1.-1e-6,"initial boundary margin");
    }
    check(rng==positionOnly,"cold initialization must draw positions only");

    // Three velocity components should each have code variance 1/4:
    // kT/m = 4 m^2/s^2, v0 = 4 m/s. This checks the physical conversion.
    config.Te_eV=8.; config.numParticles=20001;
    ParticleType thermal(config.name,config.m,config.q);
    initializePicSpecies(thermal,config,scales,mesh,rng,false);
    double sx=0.,sy=0.,sz=0.,sxx=0.,syy=0.,szz=0.;
    for(const auto& p:thermal.particle) {
        const double x=p.velocity.x-2.,y=p.velocity.y+1.,z=p.vz;
        sx+=x;sy+=y;sz+=z;sxx+=x*x;syy+=y*y;szz+=z*z;
    }
    const double n=thermal.number;
    check(n==config.numParticles && std::abs(sx/n)<.02 && std::abs(sy/n)<.02 && std::abs(sz/n)<.02,
        "thermal means around physical drift");
    check(std::abs(sxx/n-.25)<.015 && std::abs(syy/n-.25)<.015 && std::abs(szz/n-.25)<.015,
        "three dimensional thermal variance in code units");

    std::mt19937_64 allRng(43210),filteredRng(allRng);
    int rejected=0;
    for(int index=0;index<4;++index) {
        config.numParticles=7+index*2;
        config.Te_eV=(index==1)?-2.:8.;
        config.m=(index==3)?0.:2.*PhysicalConstants::eV_J;
        ParticleType all(config.name,config.m,config.q),filtered(config.name,config.m,config.q);
        initializePicSpecies(all,config,scales,mesh,allRng,false);
        int visits=0;
        initializePicSpecies(filtered,config,scales,mesh,filteredRng,false,[&](const Vector2D& point,int cell) {
            ++visits;check(cell==mesh.findCellLabel(point),"exclusion sees the sampled cell");
            if(point.x<.5) {++rejected;return true;} return false;
        });
        check(visits==config.numParticles && allRng==filteredRng,"exclusion must not retry or change the RNG stream");
        std::size_t kept=0;
        for(const auto& p:all.particle) if(p.position.x>=.5) {
            const auto& q=filtered.particle.at(kept++);
            check(p.position.x==q.position.x && p.position.y==q.position.y
                && p.velocity.x==q.velocity.x && p.velocity.y==q.velocity.y && p.vz==q.vz
                && p.w==q.w && p.inCellL==q.inCellL && p.inBackgroundCellL==q.inBackgroundCellL,
                "filter preserves surviving samples across odd hot/cold species");
        }
        check(kept==filtered.particle.size() && kept==filtered.number,"filtered initial count");
        if(index==3) {
            check(all.kappa_s()==0.,"nonpositive mass acceleration fallback");
            for(const auto& p:all.particle) check(p.velocity.x==2. && p.velocity.y==-1. && p.vz==0.,"nonpositive mass thermal fallback");
        }
    }
    check(rejected>0,"particle exclusion exercised");
    config.m=2.*PhysicalConstants::eV_J;config.numParticles=31;
    ParticleType restored(config.name,config.m,config.q);
    const auto beforeRestore=rng;int exclusionCalls=0;
    initializePicSpecies(restored,config,scales,mesh,rng,true,[&](const auto&,int) {++exclusionCalls;return false;});
    check(restored.number==0 && restored.particle.empty() && restored.particle.capacity()>=31
        && rng==beforeRestore && exclusionCalls==0,"checkpoint initialization must not sample or exclude particles");
    near(restored.kappa_s(),-.375,"checkpoint species normalization");
    config.numParticles=-1;
    ParticleType empty(config.name,config.m,config.q);
    initializePicSpecies(empty,config,scales,mesh,rng,false);
    check(empty.number==0 && rng==beforeRestore,"nonpositive initial count draws nothing");
}
using TrajectoryRow = std::array<double,14>;
std::vector<TrajectoryRow> trajectoryRows(const std::filesystem::path& path) {
    std::ifstream input(path);check(bool(input),"missing trajectory file");
    std::string line;std::getline(input,line);check(line=="Title = \"Particle Trajectories\"","trajectory title");
    std::getline(input,line);check(line.find("\"initial_index\"")!=std::string::npos,"trajectory columns");
    std::vector<TrajectoryRow> rows;
    while(std::getline(input,line)) {
        if(line.rfind("ZONE ",0)==0)continue;
        std::istringstream values(line);TrajectoryRow row;
        for(auto& value:row)check(bool(values>>value),"truncated trajectory record");
        std::string extra;check(!(values>>extra),"extra trajectory columns");rows.push_back(row);
    }
    return rows;
}
void trajectories(const std::filesystem::path& root) {
    std::vector<ParticleType> species;
    species.emplace_back("Neg - Species",1.,-1.);species.emplace_back("Positive",2.,1.);species.emplace_back("Neutral",3.,0.);
    for(int si=0;si<3;++si)for(int pi=0;pi<3;++pi) {
        Particle p;p.position={.1*(pi+1),.2*(si+1)};p.velocity={double(pi+1),-double(si+1)};
        p.vz=pi+si+.5;p.w=3.+pi;p.inCellL=pi+si*10;species[si].particle.push_back(p);
    }
    species[0].particle[0].trackId=900;
    std::mt19937_64 rng(7251);const auto initialRng=rng;
    PicParticleTrajectoryConfig off;off.baseDir="disabled";
    std::filesystem::create_directories(root/off.baseDir);std::ofstream(root/off.baseDir/"keep")<<"keep";
    PicParticleTrajectories disabled;disabled.assignMissingIds(species);
    check(species[0].particle[0].trackId==900 && species[0].particle[1].trackId==1
        && species[2].particle[2].trackId==8,"legacy track ID assignment order/preservation");
    for(const auto& sp:species)check(sp.number==3,"ID assignment synchronizes species count");
    disabled.configure(off,1,root,"contract");disabled.resetOutput();
    disabled.sample(species,rng,1,.5,2.,4.);disabled.finish();
    check(rng==initialRng && std::filesystem::exists(root/off.baseDir/"keep"),"disabled trajectories changed RNG or output");
    std::map<std::uint64_t,Particle> particles;
    for(const auto& sp:species)for(const auto& p:sp.particle)particles.emplace(p.trackId,p);

    PicParticleTrajectoryConfig config;config.enabled=true;config.baseDir="active";
    config.electronCount=config.ionCount=1;config.electronMaxSamples=config.ionMaxSamples=5;
    config.electronEvery=2;config.ionEvery=3;config.electronSampleEvery=1;config.ionSampleEvery=2;
    const auto dir=root/config.baseDir;
    std::filesystem::create_directories(dir/"old");std::ofstream(dir/"old/stale")<<"old";
    PicParticleTrajectories recorder;recorder.configure(config,7,root,"contract");recorder.resetOutput();
    check(std::filesystem::is_empty(dir),"enabled startup must clear prior trajectory files");
    recorder.sample(species,rng,1,.5,2.,4.);check(rng==initialRng,"unscheduled batch consumed RNG");
    recorder.sample(species,rng,2,1.,2.,4.);check(rng!=initialRng,"random selection must use the supplied RNG");
    check(std::filesystem::is_empty(dir),"active trajectory written before rollover/completion");
    std::reverse(species[0].particle.begin(),species[0].particle.end());
    recorder.sample(species,rng,3,1.5,2.,4.);
    recorder.sample(species,rng,4,2.,2.,4.);
    recorder.sample(species,rng,5,2.5,2.,4.);recorder.finish();
    const auto e2=trajectoryRows(dir/"particle_trajectories_electron_0000002.dat");
    const auto e4=trajectoryRows(dir/"particle_trajectories_electron_0000004.dat");
    const auto i3=trajectoryRows(dir/"particle_trajectories_ion_0000003.dat");
    check(e2.size()==2 && e4.size()==2 && i3.size()==1,"batch rollover/final flush/sample cadence");
    check(e2[0][0]==2 && e2[1][0]==3 && e4[0][0]==4 && e4[1][0]==5 && i3[0][0]==4,"sampled step numbers");
    check(e2[0][3]==e2[1][3],"particle reordering changed tracked identity");
    for(const auto* rows:{&e2,&e4,&i3})for(const auto& row:*rows) {
        const auto& p=particles.at(static_cast<std::uint64_t>(row[3]));
        check(row[1]==row[0]*.5 && row[2]==1 && row[5]==1,"sample time/slot/alive");
        near(row[6],p.position.x*2.,"trajectory length scale");near(row[7],p.position.y*2.,"trajectory y scale");
        near(row[8],p.velocity.x*4.,"trajectory velocity scale");near(row[9],p.velocity.y*4.,"trajectory vy scale");
        near(row[10],p.vz*4.,"trajectory vz scale");
        check(row[11]==p.inCellL && row[12]==p.w,"trajectory mesh/weight");
        check(row[4]==(rows==&i3?1.:0.),"charge role or electron/ion ordering");
    }
    bool rejected=false;try {recorder.configure(config,1,root,"contract");}catch(const std::logic_error&){rejected=true;}
    check(rejected,"run-owned trajectory state must not be reconfigured");

    config.ionCount=0;config.electronEvery=10;config.baseDir="lost";
    PicParticleTrajectories lost;lost.configure(config,1,root,"contract");lost.resetOutput();
    lost.sample(species,rng,10,5.,2.,4.);const auto saved=species[0].particle;species[0].particle.clear();
    lost.sample(species,rng,11,5.5,2.,4.);species[0].particle=saved;
    lost.sample(species,rng,12,6.,2.,4.);lost.finish();
    check(trajectoryRows(root/config.baseDir/"particle_trajectories_electron_0000010.dat").size()==1,
        "lost particle should finish without a fabricated dead sample or later resumption");
    config.baseDir="limited";config.electronEvery=1;config.electronSampleEvery=2;config.electronMaxSamples=1;
    {
        PicParticleTrajectories limited;limited.configure(config,1,root,"contract");limited.resetOutput();
        limited.sample(species,rng,1,.5,2.,4.);limited.sample(species,rng,2,1.,2.,4.);
        check(trajectoryRows(root/config.baseDir/"particle_trajectories_electron_0000001.dat").empty(),"unsampled rollover retains empty header");
        check(trajectoryRows(root/config.baseDir/"particle_trajectories_electron_0000002.dat").size()==1,"maximum samples must finish the batch");
        limited.sample(species,rng,3,1.5,2.,4.);
    }
    check(!std::filesystem::exists(root/config.baseDir/"particle_trajectories_electron_0000003.dat"),"destructor must not flush an unfinished batch");
    config.baseDir="no_candidates";config.electronEvery=1;config.electronSampleEvery=1;
    std::vector<ParticleType> empty;const auto beforeEmpty=rng;
    PicParticleTrajectories noCandidates;noCandidates.configure(config,1,root,"contract");noCandidates.resetOutput();
    noCandidates.sample(empty,rng,1,.5,2.,4.);
    check(rng==beforeEmpty && trajectoryRows(root/config.baseDir/"particle_trajectories_electron_0000001.dat").empty(),"empty selection consumes no RNG and writes only header");
}
}
namespace {
void boundaryRates()
{
    using Losses=PicBoundaryLosses;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;
    input.maxCoord={2.,3.};QuadtreeMesh mesh(input);
    UnitsScales scales;scales.L0_m=2.;scales.T0_s=.5;
    auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(scales));
    NormalizationConfig config;config.slabThickness_m=3.;
    const auto normalization=NormalizationContext::MakeFromUnits(units,config);
    std::vector<ParticleType> species;
    species.emplace_back("electron",2.,-3.);species.emplace_back("ion",2.,3.);
    species.emplace_back("neutral",2.,0.);
    Losses losses(3,2.,false,false,Losses::Mode::Implicit);
    Particle p;p.velocity={1.,2.};p.vz=2.;p.rngState=123;
    for(int side=0;side<4;++side)for(std::size_t si=0;si<3;++si) {
        p.w=double(si+2);losses.recordDomain(si,species[si],p,1.,side);
    }
    const auto wall=losses.summarizeImplicitWall(species,mesh,units.get(),&normalization,.5);
    check(wall.areaM2==std::array<double,4>{18.,18.,12.,12.},"physical side area uses length scale and slab thickness");
    for(int side=0;side<4;++side) {
        check(wall.electronCharge[side]==-6. && wall.ionCharge[side]==9.,"signed charge excludes neutrals");
        near(wall.electronCurrent[side]*wall.areaM2[side]*.5,-6.,"electron current integrates to lost charge");
        near(wall.ionCurrent[side]*wall.areaM2[side]*.5,9.,"ion current integrates to lost charge");
        const auto& t=losses.totals();
        const auto flux=PicBoundaryFlux<std::uint64_t>::sample(t.domainMacro[0][side],t.domainPhysical[0][side],-3.,wall.areaM2[side],.5);
        near(flux.current,wall.electronCurrent[side],"wide/long species current agrees with charge accounting");
        const auto energy=PicBoundaryEnergy::sample({1,2.},t.domainEnergy[0][side],wall.areaM2[side],.5);
        check(energy.out.macro==1 && energy.absorbedEnergyJ==72. && energy.absorbedPowerW==144.,"absorption power from physical three-component energy");
        near(energy.energyFluxWpm2*energy.surfaceAreaM2,144.,"energy flux integrates to absorbed power");
        near(energy.meanEnergyEV*PhysicalConstants::eV_J,36.,"mean absorbed energy in electron volts");
    }
    const auto legacy=losses.summarizeImplicitWall(species,mesh,nullptr,nullptr,2.);
    check(legacy.areaM2==std::array<double,4>{3.,3.,2.,2.},"missing contexts retain unit slab and code lengths");
    const auto clamped=losses.summarizeImplicitWall(species,mesh,nullptr,nullptr,0.);
    near(clamped.electronCurrent[0]*3e-300,-6.,"implicit nonpositive dt retains 1e-300 clamp");
    check(losses.totals().domainChargeE==std::array<double,4>{} && p.rngState==123,"reporting does not mutate losses or RNG");
    bool rejected=false;
    try { Losses(3,1.,false,false,Losses::Mode::Explicit).summarizeImplicitWall(species,mesh,nullptr,nullptr,1.); }
    catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"explicit event charge cannot be silently resampled as implicit species charge");
    rejected=false;
    try { Losses(2,1.,false,false,Losses::Mode::Implicit).summarizeImplicitWall(species,mesh,nullptr,nullptr,1.); }
    catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"species mismatch rejected before reduction");
    std::vector<ParticleType> decimal;decimal.emplace_back("decimal",1.,-.1);
    Losses reduced(1,1.,false,false,Losses::Mode::Implicit);double eventCharge=0.;p.w=1.;
    for(int i=0;i<10;++i)eventCharge+=reduced.recordDomain(0,decimal[0],p,1.,0).chargeC;
    check(reduced.summarizeImplicitWall(decimal,mesh,nullptr,nullptr,1.).electronCharge[0]==-1.
        && eventCharge!=-1.,"implicit species-first reduction must not become event-wise summation");
    const std::uint64_t large=(std::uint64_t{1}<<53)+1;
    const auto integer=PicBoundaryFlux<std::uint64_t>::sample(large,12.,-2.,3.,2.);
    check(integer.macro==large && integer.physicalFlux==2. && integer.current==-4.,"large integer macro count survives normalization");
    const auto interval=PicBoundaryFlux<double>::sample(6.,12.,-2.,3.,2.);
    const auto instant=PicBoundaryFlux<double>::sample(2.,4.,-2.,3.,.5);
    check(interval.current==-4. && instant.current!=interval.current,"report interval is distinct from instantaneous step");
    const double nan=std::numeric_limits<double>::quiet_NaN();
    for(double area:{0.,-1.,nan}) {
        const auto flux=PicBoundaryFlux<double>::sample(1.,2.,-3.,area,1.);
        const auto energy=PicBoundaryEnergy::sample({1,2.},6.,area,1.);
        check(flux.current==0. && flux.physical==2. && std::isnan(energy.energyFluxWpm2),"undefined surface preserves distinct particle/energy flux semantics");
    }
    for(double dt:{0.,-1.,nan})
        check(PicBoundaryFlux<double>::sample(1.,2.,-3.,2.,dt).current==0.,"instantaneous nonpositive/nonfinite dt has zero flux");
    check(PicBoundaryEnergy::sample({0,0.},0.,1.,1.).meanEnergyEV==0.,"empty energy sample has zero mean");
}
void boundaryLossAccounting()
{
    using Losses=PicBoundaryLosses;using Mode=Losses::Mode;
    static_assert(std::is_same_v<decltype(std::declval<Losses&>().totals()),const Losses::Totals&>);
    ParticleType electron("electron",2.,-3.),ion("ion",2.,3.),neutral("neutral",2.,0.);
    Particle p;p.velocity={1.,2.};p.vz=2.;p.w=2.;p.rngState=123;
    for(auto mode:{Mode::Explicit,Mode::Implicit}) {
        Losses losses(3,2.,false,false,mode);const auto& totals=losses.totals();
        for(int side=0;side<4;++side) {
            p.w=2.;const auto e=losses.recordDomain(0,electron,p,7.,side);
            p.w=0.;const auto i=losses.recordDomain(1,ion,p,3.,side);
            p.w=-1.;const auto n=losses.recordDomain(2,neutral,p,5.,side);
            check(e.accepted && i.accepted && n.accepted && e.weight==2. && i.weight==3. && n.weight==5.,
                "accepted domain particles and configured weight fallback");
            near(e.energyPerParticleJ,36.,"three-component SI kinetic energy");
            check(e.chargeC==-6. && i.chargeC==9. && n.chargeC==0.,"signed impact charge");
            check(totals.domainMacroE[side]==1 && totals.domainMacroI[side]==1
                && totals.domainMacro[2][side]==1,"charged/neutral macro counting");
            near(totals.domainEnergy[0][side],72.,"electron absorption energy");
            near(totals.domainEnergy[1][side],108.,"ion absorption energy");
            near(totals.domainEnergy[2][side],180.,"neutral absorption energy");
            check(totals.domainChargeE[side]==(mode==Mode::Explicit?-6.:0.)
                && totals.domainChargeI[side]==(mode==Mode::Explicit?9.:0.),"implicit side charge remains a later species reduction");
        }
        near(totals.domainEnergyJ,1440.,"domain event sum");
        check(!losses.recordDomain(0,electron,p,7.,-1).accepted
            && !losses.recordDomain(1,ion,p,7.,4).accepted
            && !losses.recordDomain(2,neutral,p,7.,99).accepted,"unknown side cannot become physical outflow");
        check(totals.unknownE==1 && totals.unknownI==1 && totals.domainEnergyJ==1440.,"unknown counters without energy double count");
        p.w=2.;const auto e=losses.recordEmbedded(0,electron,p,7.);
        p.w=0.;const auto i=losses.recordEmbedded(1,ion,p,3.);
        p.w=-1.;const auto n=losses.recordEmbedded(2,neutral,p,5.);
        losses.recordEmbeddedBoundary(7,e,true,true);
        losses.recordEmbeddedBoundary(-1,i,false,true);
        losses.recordEmbeddedBoundary(8,n,true,false);
        check(totals.embeddedMacro==std::vector<long long>({1,1,1})
            && totals.embeddedPhysical==std::vector<double>({2.,3.,5.}),"owned per-species embedded loss");
        near(totals.embeddedEnergyJ,360.,"embedded total independent of domain");
        check(totals.embeddedChargeById.at(7)==-6. && totals.embeddedChargeById.at(8)==0.
            && totals.embeddedChargeById.size()==2,"only resolved positive IDs accumulate charge, including zero-charge neutral");
        check(totals.embeddedCounts.byId.at(7)==1,"resolved impact count");
        if(mode==Mode::Explicit)
            check(totals.embeddedCounts.inside==1 && totals.embeddedCounts.crossing==1 && totals.embeddedCounts.unknown==1
                && totals.embeddedCounts.byId.count(8)==0,"explicit classification cadence independent of electrical accounting");
        else check(totals.embeddedCounts.byId.at(8)==1 && totals.embeddedCounts.unknown==0
                && totals.embeddedCounts.inside==0 && totals.embeddedCounts.crossing==0,"implicit unconditional ID count and original classification scope");
        Losses next(3,2.,false,false,mode);
        check(next.totals().domainEnergyJ==0. && next.totals().embeddedChargeById.empty()
            && next.totals().embeddedMacro==std::vector<long long>(3,0),"fresh macro-step owns fresh counters");

        for(bool px:{false,true})for(bool py:{false,true}) {
            Losses periodic(2,2.,px,py,mode);
            for(int side=0;side<4;++side) {
                const bool excluded=side<2?px:py;
                const auto e=periodic.recordDomain(0,electron,p,1.,side);
                const auto i=periodic.recordDomain(1,ion,p,1.,side);
                const auto& t=periodic.totals();
                check(e.accepted!=excluded && i.accepted!=excluded && t.domainMacroE[side]==!excluded
                    && t.domainMacroI[side]==!excluded,"periodic impacts excluded on configured axis only");
                check(t.periodicLeakE[side]==(excluded && mode==Mode::Explicit)
                    && t.periodicLeakI[side]==(excluded && mode==Mode::Explicit),"explicit-only periodic leak accounting");
            }
        }
        p.w=std::numeric_limits<double>::infinity();
        Losses invalidWeight(1,2.,false,false,mode);
        const auto weighted=invalidWeight.recordDomain(0,electron,p,3.,0);
        if(mode==Mode::Explicit)check(std::isinf(weighted.weight) && std::isinf(invalidWeight.totals().domainEnergyJ),"explicit infinite weights retained");
        else check(weighted.weight==3. && invalidWeight.totals().domainEnergyJ==108.,"implicit finite-weight fallback");
        p.w=-1.;Losses nanWeight(1,2.,false,false,mode);
        nanWeight.recordEmbedded(0,electron,p,std::numeric_limits<double>::quiet_NaN());
        check(std::isnan(nanWeight.totals().embeddedPhysical[0]),"invalid default weight retains physical-count semantics");
        check(mode==Mode::Explicit?std::isnan(nanWeight.totals().embeddedEnergyJ):nanWeight.totals().embeddedEnergyJ==0.,"mode-specific weighted energy filtering");
        p.w=2.;p.vz=std::numeric_limits<double>::infinity();
        Losses invalidVelocity(1,2.,false,false,mode);
        const auto impact=invalidVelocity.recordEmbedded(0,electron,p,3.);
        check(std::isinf(impact.energyPerParticleJ),"SEE input keeps original per-particle energy even when aggregate filters it");
        check(mode==Mode::Explicit?std::isinf(invalidVelocity.totals().embeddedEnergyJ):invalidVelocity.totals().embeddedEnergyJ==0.,"mode-specific finite energy guard");
        p.vz=2.;Losses negativeScale(1,-2.,false,false,mode);
        negativeScale.recordDomain(0,electron,p,3.,0);
        check(negativeScale.totals().domainEnergyJ==(mode==Mode::Explicit?72.:0.),"mode-specific velocity scale guard");
        check(p.rngState==123 && p.velocity.x==1. && p.velocity.y==2.,"accounting does not mutate the particle");
    }
}
}
int main() try
{
    resamplingStage();
    populationStage();
    boundaryRates();
    boundaryLossAccounting();
    spawnManagement();
    particleSources();
    secondaryEmission();
    implicitEnvironmentInput();
    explicitAmrSetup();
    implicitAmrSetup();
    segmenterLookup();
    implicitFieldSolve();
    implicitFieldAdvance();
    speciesPush();
    speciesInitializationAndRelabel();
    poissonStage();
    explicitFieldAdvance();
    chargeDeposition();
    fieldBoundaryState();
    refinementState();
    refinementInputs();
    currentWeights();
    couplingPolicy();
    nestedCouplingChunks();
    orbitConfiguration();
    orbitStage();
    polarizationStage();
    boundarySetup();
    DictionaryMesh input; input.Nx=input.Ny=8; input.maxCellNum=256; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    initialization(mesh);
    Particle p; p.position={.3125,.3125}; p.velocity={-3.,4.}; p.vz=12.; p.w=2.5;
    p.inCellL=mesh.findCellLabel(p.position); p.inBackgroundCellL=p.inCellL;
    p.discreteVertL=3; p.emitChance=.25; p.rngState=0x123456789abcdef0ULL; p.trackId=123;
    p.levelSetValues={-.125,.25,.5};
    p.updateWeight(mesh,mesh.cells(p.inCellL));
    for (double w:p.weight) near(w,.25,"bilinear weights at cell center");
    Particle copy(p); copy.levelSetValues[0]=99.;
    check(copy.rngState==p.rngState && copy.trackId==123 && p.levelSetValues[0]==-.125,"particle copy ownership/identity");

    ParticleType original("contract",2.,3.); original.set_kappa_s(.125);
    original.particle={p}; original.number=1;
    const auto directory=std::filesystem::current_path()/
        ("particle_contract_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    explicitAmrStage(directory);
    implicitAmrStage(directory);
    normalizationInputs(directory);
    speciesRunStartup(directory);
    implicitFieldRebuild(directory);
    trajectories(directory);
    original.writeCache(directory.generic_string()+"/");
    ParticleType restored("contract",9.,-9.); restored.readCache(directory.generic_string()+"/");
    check(restored.number==1 && restored.particle.size()==1 && restored.name()=="contract","cache species identity");
    near(restored.mass(),2.,"mass cache"); near(restored.charge(),3.,"charge cache");
    near(restored.q_m(),1.5,"charge/mass cache"); near(restored.kappa_s(),.125,"normalization cache");
    const auto& r=restored.particle[0];
    near(r.position.x,p.position.x,"position cache x"); near(r.position.y,p.position.y,"position cache y");
    near(r.velocity.x,p.velocity.x,"velocity cache x"); near(r.velocity.y,p.velocity.y,"velocity cache y");
    near(r.vz,p.vz,"velocity cache z"); near(r.w,p.w,"macro weight cache"); near(r.emitChance,p.emitChance,"emission cache");
    for (int i=0;i<4;++i) near(r.weight[i],p.weight[i],"interpolation cache");
    check(r.rngState==p.rngState && r.levelSetValues==p.levelSetValues && r.inCellL==p.inCellL
        && r.inBackgroundCellL==p.inBackgroundCellL && r.discreteVertL==p.discreteVertL,"cache geometry/RNG");
    check(r.trackId==0,"existing v4 cache does not store runtime track IDs");

    ParticleType population("ions",2.,3.); population.particle={p,p,p,p}; population.number=4;
    for (int i=0;i<4;++i) population.particle[i].rngState=10+i;
    const auto initial=totals(population);
    ParticleManagementConfig config; config.enabled=true; config.perCellMax=1; config.splitEnabled=false;
    ParticleManagementStats stats; std::mt19937_64 rng(24680);
    ParticleManager::Apply(config,mesh,population,0,2.,1.,rng,stats);
    check(population.number==1 && stats.mergedPairs==3,"actual particle merge");
    conserved(initial,totals(population));
    config.mergeEnabled=false; config.splitEnabled=true; config.splitWMax=5.; config.splitMaxPerCell=1;
    ParticleManager::Apply(config,mesh,population,0,2.,1.,rng,stats);
    check(population.number==1 && stats.splitCount==0,"legacy single-particle weight-only split guard");
    config.splitByCountEnabled=true; config.perCellMin=2;
    ParticleManager::Apply(config,mesh,population,0,2.,1.,rng,stats);
    check(population.number==2 && stats.splitCount==1 && stats.splitEmergencyCount==1,"actual particle split");
    check(population.particle[0].rngState!=population.particle[1].rngState,"split RNG identity");
    conserved(initial,totals(population));
    config.mergeEnabled=true; config.splitEnabled=false;
    std::vector<unsigned char> protectedCells(mesh.cellIndexNum(),1);
    ParticleManager::Apply(config,mesh,population,0,2.,1.,rng,stats,true,&protectedCells);
    check(population.number==2 && stats.mergedPairs==3,"protected cell must not merge");
    conserved(initial,totals(population));

    DictionaryGeometry dict;
    Line wall; wall.points[0]={.5,0.}; wall.points[1]={.5,1.}; wall.setType(LineType::straight);
    dict.lines={wall}; dict.levelSetName={"wall"}; dict.innerPoints={{.25,.5}};
    dict.linesOnInterface={{{0,false}}}; dict.moving={false};
    ImmersedGeometry<QuadtreeMesh> geometry(dict,mesh); geometry.lines[0].levelSetIndex=0;
    for (int i=0;i<mesh.vertIndexNum();++i)
        geometry.levelSets[0].value[i]=mesh.verts(mesh.vertIndexToLabel[i]).coord().x-.5;
    geometry.refreshLevelSetDerivatives(mesh);
    WallImpactSample impact;
    const Vector2D impactPoint(.5,.4);
    const auto resolve=[&](double charge) {
        return MovingEBController::resolveWallImpact(mesh,geometry,nullptr,{7},nullptr,impactPoint,
            mesh.findCellLabel(impactPoint),p,2,2.*1.602176634e-19,charge,1.,2.,.75,true,impact);
    };
    check(resolve(1.) && impact.resolved && impact.boundaryId==7 && impact.levelSetIndex==0,"ion wall resolution");
    near(impact.normal.x,-1.,"wall normal"); near(impact.vn,6.,"physical impact speed");
    near(impact.energy_eV,676.,"3D velocity-space impact energy");
    near(impact.theta_rad,std::acos(3./13.),"impact angle"); near(impact.macroWeight,2.5,"impact weight");
    near(impact.time_s,.75,"impact time");
    check(!resolve(-1.) && !impact.ionSpecies,"electron excluded from ion erosion");
    std::cout<<"Particle state/cache, population conservation and wall-impact checks passed\n";
    return 0;
}
catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
