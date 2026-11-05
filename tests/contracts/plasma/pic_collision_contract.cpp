#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/Physics/BuiltInCrossSections.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/IO/PicCheckpointWriter.h"
#include "Modules/Plasma/IO/PicCollisionDiagnostics.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include "Modules/Plasma/Physics/PicParticleSources.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/DictionaryMesh.h"
#include <chrono>
#include <array>
#include <limits>
#include <fstream>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <sstream>

namespace {
template<class T,class=void> struct PublicRngRestore : std::false_type {};
template<class T> struct PublicRngRestore<T,std::void_t<decltype(&T::restoreRngState)>> : std::true_type {};
template<class T,class=void> struct PublicCollisionRead : std::false_type {};
template<class T> struct PublicCollisionRead<T,std::void_t<decltype(&T::readCollisions)>> : std::true_type {};
static_assert(!PublicRngRestore<PicCollisions>::value);
static_assert(!PublicCollisionRead<PicCheckpointReader>::value);
void check(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
SpeciesConfig electron()
{
    SpeciesConfig species;
    species.name = "electron";
    species.q = -PhysicalConstants::eV_J;
    species.m = 9.1093837015e-31;
    species.weight = 7.0;
    species.Te_eV = 1.0;
    return species;
}
PlasmaCollisionInputs inputs(CollisionsRngMode mode, bool ionization=false)
{
    BackgroundGas gas("Ar", 1000.0, 300.0);
    gas.SetParticleMassKg(6.6335209e-26);
    auto xs = std::make_shared<CrossSectionDB>();
    std::string error;
    check(BuiltInCrossSections::RegisterArV1(xs.get(), "", false, &error), "cross sections");
    auto reactions = std::make_shared<ReactionDatabase>();
    ReactionSet set;
    set.name = "electron-Ar";
    set.projectile = "electron";
    set.target = "Ar";
    ReactionChannel channel;
    channel.name = "elastic";
    channel.cross_section_key = "e-Ar:ELASTIC";
    channel.scatter_model = "isotropic";
    if (ionization) {
        channel.name="ionization";channel.cross_section_key="e-Ar:IONIZATION";
        channel.threshold_eV=15.76;channel.products={"electron","electron","ion"};
    }
    set.channels.push_back(channel);
    check(reactions->Register(set, false, &error), "reaction set");
    auto config = std::make_shared<CollisionsConfig>();
    config->enabled = true;
    config->parallel = false;
    config->rng_mode = mode;
    config->rng_seed = 24680;
    config->mappings.push_back({"electron", "electron-Ar"});
    if(ionization) {
        ReactionSet ions;ions.name="ion-Ar";ions.projectile="ion";ions.target="Ar";
        ReactionChannel elastic;elastic.name="elastic";elastic.cross_section_key="Ar+-Ar:ELASTIC";elastic.scatter_model="isotropic";
        ions.channels.push_back(elastic);check(reactions->Register(ions,false,&error),"ion reaction set");
        config->mappings.push_back({"ion","ion-Ar"});
    }
    return {std::make_shared<GasState>(gas.Evaluate()), xs, reactions, config};
}
void compareAdvance(CollisionsRngMode mode, double velocityScale)
{
    auto ownedInputs = inputs(mode);
    const std::weak_ptr<const GasState> gasOwner = ownedInputs.gas;
    const std::weak_ptr<const CrossSectionDB> xsOwner = ownedInputs.crossSections;
    const std::weak_ptr<const ReactionDatabase> reactionOwner = ownedInputs.reactions;
    const std::weak_ptr<const CollisionsConfig> configOwner = ownedInputs.collisions;
    std::string error;
    auto service = PicCollisions::create(ownedInputs, {electron()}, velocityScale, error);
    check(service != nullptr && error.empty(), "service initialization");
    check(service->rngState().applyCount==0,"new collision owner starts at zero");
    ownedInputs = {};
    check(!gasOwner.expired() && !xsOwner.expired() && !reactionOwner.expired() && !configOwner.expired(),
        "operator dependencies outlived caller inputs");
    const auto* stableOperator = &service->operation();
    auto movedOwner = std::move(service);
    check(&movedOwner->operation() == stableOperator, "operator address changed with owner transfer");

    // Direct existing operator is the reference. It independently owns neither
    // the new PIC service nor its inputs; preserve its old per-call contract.
    auto referenceInputs = inputs(mode);
    SpeciesDatabase database;
    Species definition;
    definition.name = "electron";
    definition.mass_kg = electron().m;
    definition.charge_C = electron().q;
    definition.macro_weight = 7.0;
    definition.role = SpeciesRole::kElectron;
    check(database.Add(definition, &error), "reference species");
    MCCCollisionOperator direct(database, *referenceInputs.gas, *referenceInputs.crossSections, *referenceInputs.reactions);
    check(direct.Initialize(*referenceInputs.collisions, &error), "reference initialization");
    ParticleType actual("electron", electron().m, electron().q);
    ParticleType expected("electron", electron().m, electron().q);
    for (int i = 0; i < 512; ++i)
    {
        Particle particle;
        particle.position = Vector2D(0.25, 0.75);
        particle.velocity = Vector2D(524288.0 / velocityScale, 131072.0 / velocityScale);
        particle.vz = 65536.0 / velocityScale;
        particle.w = 7.0;
        particle.trackId = i + 1;
        actual.particle.push_back(particle);
    }
    actual.number = 512;
    expected.particle = actual.particle;
    expected.number = 512;
    std::mt19937_64 rng(12345), referenceRng(12345);
    std::unique_ptr<PicCollisions> resumed;
    ParticleType restarted("electron", electron().m, electron().q);
    std::mt19937_64 restoredRng;
    std::uint64_t collisions = 0;
    for (int step = 0; step < 4; ++step)
    {
        if (step == 2)
        {
            const auto root=std::filesystem::current_path()/
                ("collision_restart_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
            PicCheckpointStore store(root,{});std::filesystem::path dir;
            check(store.createDirectory(2,dir) && PicCheckpointStore::writeMetadata(dir,2,2e-11,1e-11),"collision checkpoint fixture");
            PicCheckpointWriter writer(store,rng,true,0,"contract");
            writer.writeCollisions(dir,*movedOwner);
            double checkpointDt=1e-11;
            auto reader=PicCheckpointReader::select(store,2,checkpointDt,"contract");check(bool(reader),"collision checkpoint selection");
            resumed=PicCollisions::create(referenceInputs,{electron()},velocityScale,error,&*reader,true);
            check(resumed && resumed->rngState().applyCount==2 && resumed->rngState().mode==mode,
                "factory publishes restored collision call count");
            const auto path=dir/PicCheckpointStore::collisionName;
            for(const auto& bad:std::vector<std::string>{
                "LAMP_MCC_RNG 2\n2 24680 99", "LAMP_MCC_RNG 1\n2 24680 -1",
                "LAMP_MCC_RNG 1\n2 24680 18446744073709551616", "LAMP_MCC_RNG 1\n2 24680",
                "LAMP_MCC_RNG 1\n2 24680 99 extra", "LAMP_MCC_RNG 1\n3 24680 99",
                "LAMP_MCC_RNG 1\n"+std::to_string(static_cast<int>(mode))+" 123 99",
                "LAMP_MCC_RNG 1\n"+std::to_string((static_cast<int>(mode)+1)%3)+" 24680 99"})
            {
                {std::ofstream file(path);file<<bad;}
                std::unique_ptr<PicCollisions> candidate;bool rejected=false;
                try {candidate=PicCollisions::create(referenceInputs,{electron()},velocityScale,error,&*reader,true);}
                catch(const std::exception&) {rejected=true;}
                check(rejected && !candidate && resumed->rngState().applyCount==2,
                    "failed restore must not publish a fresh collision owner or mutate the live owner");
            }
            auto fresh=PicCollisions::create(referenceInputs,{electron()},velocityScale,error,&*reader,false);
            check(fresh && fresh->rngState().applyCount==0,"restoreRng=false ignores invalid collision checkpoint");
            auto incomplete=referenceInputs;incomplete.gas.reset();
            check(!PicCollisions::create(incomplete,{electron()},velocityScale,error,&*reader,true) && !error.empty(),
                "failed initialization must precede checkpoint restore");
            auto disabled=referenceInputs;disabled.collisions=std::make_shared<CollisionsConfig>();
            check(!PicCollisions::create(disabled,{},velocityScale,error,&*reader,true) && error.empty(),
                "disabled collisions must not read the selected checkpoint");
            std::filesystem::remove(path);
            fresh.reset();bool rejected=false;
            try {fresh=PicCollisions::create(referenceInputs,{electron()},velocityScale,error,&*reader,true);}
            catch(const std::exception&) {rejected=true;}
            check(mode==CollisionsRngMode::kCounterBased ? rejected && !fresh
                : !rejected && fresh && fresh->rngState().applyCount==0,
                "factory preserves legacy missing-component policy");
            PicCheckpointWriter noRng(store,rng,false,0,"contract");noRng.writeCollisions(dir,*resumed);
            fresh=PicCollisions::create(referenceInputs,{electron()},velocityScale,error,&*reader,false);
            check(!std::filesystem::exists(path) && fresh && fresh->rngState().applyCount==0
                && resumed->rngState().applyCount==2,"saveRng gate applies to the operator-owned state");
            writer.writeCollisions(dir,*movedOwner);
            resumed=PicCollisions::create(referenceInputs,{electron()},velocityScale,error,&*reader,true);
            check(resumed && resumed->rngState().applyCount==2,"valid restore after rejection must resume the saved history");
            restarted.particle=actual.particle;restarted.number=actual.number;restoredRng=rng;
        }
        ParticleSpawnBuffer spawned, referenceSpawned;
        MCCCollisionOperator::StepStats stats, referenceStats;
        stats.Clear(); referenceStats.Clear();
        const double dt = (step % 2 == 0) ? 1e-11 : 5e-12;
        check(movedOwner->apply(actual, 7.0, dt, rng, spawned, stats, &error), "PIC apply");
        check(direct.ApplyCodeUnits(expected, 7.0, referenceInputs.gas->number_density_m3, dt,
            referenceRng, velocityScale, 1.0 / velocityScale, referenceSpawned, referenceStats, &error), "reference apply");
        check(&movedOwner->operation() == stableOperator && rng == referenceRng, "operator/RNG continuity");
        check(stats.n_collisions == referenceStats.n_collisions && stats.n_attempts == referenceStats.n_attempts
            && stats.n_null == referenceStats.n_null && stats.per_channel == referenceStats.per_channel,
            "collision events changed");
        if (resumed)
        {
            ParticleSpawnBuffer restoredSpawn;MCCCollisionOperator::StepStats restoredStats;restoredStats.Clear();
            check(resumed->apply(restarted,7.0,dt,restoredRng,restoredSpawn,restoredStats,&error),"restored collision application");
            check(restoredRng==rng && restoredStats.n_collisions==stats.n_collisions && restoredStats.per_channel==stats.per_channel
                && restoredSpawn.empty() && resumed->rngState().applyCount==movedOwner->rngState().applyCount,"resumed events and random state agree");
            for(std::size_t i=0;i<actual.particle.size();++i)
            {
                const auto& a=actual.particle[i];const auto& b=restarted.particle[i];
                check(a.velocity.x==b.velocity.x && a.velocity.y==b.velocity.y && a.vz==b.vz && a.rngState==b.rngState,
                    "collision restoration must reproduce particle results exactly");
            }
        }
        check(spawned.empty() && referenceSpawned.empty() && actual.number == 512 && actual.particle.size() == 512,
            "elastic particle budget changed");
        collisions += stats.n_collisions;
        for (std::size_t i = 0; i < actual.particle.size(); ++i)
        {
            const auto& a = actual.particle[i]; const auto& b = expected.particle[i];
            auto equal = [](double x, double y) { return std::abs(x-y) <= 1e-12 * std::max(1.0, std::abs(y)); };
            check(equal(a.velocity.x,b.velocity.x) && equal(a.velocity.y,b.velocity.y) && equal(a.vz,b.vz),
                "velocity scale or collision result changed");
            check(a.rngState == b.rngState && a.trackId == b.trackId && a.w == b.w
                && a.position.x == b.position.x && a.position.y == b.position.y, "particle state changed");
        }
    }
    check(collisions > 0, "reference did not exercise collisions");
    movedOwner.reset();
    check(gasOwner.expired() && xsOwner.expired() && reactionOwner.expired() && configOwner.expired(),
        "collision owner retained dependencies after destruction");
}
void rejectedInitialization()
{
    std::string error = "stale";
    check(!PicCollisions::create({}, {}, 1.0, error) && error.empty(), "vacuum result");
    auto source = inputs(CollisionsRngMode::kPerParticle);
    auto disabled = std::make_shared<CollisionsConfig>();
    check(!PicCollisions::create({nullptr,nullptr,nullptr,disabled}, {}, 1.0, error) && error.empty(), "disabled result");
    for (int i = 0; i < 3; ++i)
    {
        auto missing = source;
        if (i == 0) missing.gas.reset();
        if (i == 1) missing.crossSections.reset();
        if (i == 2) missing.reactions.reset();
        check(!PicCollisions::create(missing, {electron()}, 1.0, error) && !error.empty(), "partial inputs published");
    }
    auto bad = electron(); bad.m = 0.0;
    check(!PicCollisions::create(source, {bad}, 1.0, error) && !error.empty(), "invalid species accepted");
    check(!PicCollisions::create(source, {electron(),electron()}, 1.0, error) && !error.empty(), "duplicate species accepted");
    // The existing operator permits an empty registry and resolves projectile
    // mass from ParticleType during Apply. Keep that initialization behavior.
    auto emptyRegistry = PicCollisions::create(source, {}, 1.0, error);
    check(emptyRegistry && error.empty() && emptyRegistry->operation().HasProjectile("electron"), "legacy empty registry changed");
    auto badMapping = std::make_shared<CollisionsConfig>(*source.collisions);
    badMapping->mappings.front().reaction_set = "missing";
    source.collisions = badMapping;
    check(!PicCollisions::create(source, {electron()}, 1.0, error) && !error.empty(), "missing reaction accepted");
}
void collisionStage(CollisionsRngMode mode,bool implicit,bool unknownProduct)
{
    using Stage=PicCollisions::Stage;
    static_assert(!std::is_copy_constructible_v<Stage>);
    static_assert(std::is_same_v<decltype(std::declval<Stage&>().pending()),const ParticleSpawnBuffer&>);
    static_assert(std::is_same_v<decltype(std::declval<Stage&>().batchStats()),const MCCCollisionOperator::StepStats&>);
    auto ion=electron();ion.name="ion";ion.m=6.6335209e-26;ion.q=-ion.q;ion.weight=11.;
    auto config=inputs(mode,true);std::string error;
    auto owner=PicCollisions::create(config,{electron(),ion},3.,error);
    auto reference=PicCollisions::create(config,{electron(),ion},3.,error);
    check(owner && reference,"ionization stage initialization");
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=2;dictionary.maxCellNum=16;dictionary.maxInitialIterations=0;
    QuadtreeMesh mesh(dictionary);
    std::unordered_map<std::string,std::size_t> names{{"electron",0}};
    if(!unknownProduct)names.emplace("ion",1);
    const std::vector<double> weights{7.,11.};Stage stage(owner.get(),names,weights);
    ParticleManagementConfig pm;PicSpawnManagement management(implicit?PicSpawnManagement::Mode::Implicit:PicSpawnManagement::Mode::Explicit,pm,names,weights);
    std::vector<ParticleType> actual,expected;
    actual.emplace_back("electron",electron().m,electron().q);actual.emplace_back("ion",ion.m,ion.q);
    for(int j=0;j<256;++j) {
        Particle p;p.position={implicit || j%2?.25:-.25,.25};p.velocity={2000000.,0.};p.vz=100000.;
        p.w=7.;p.inCellL=implicit?mesh.findCellLabel(p.position):-1;p.inBackgroundCellL=99;p.trackId=j+1;
        actual[0].particle.push_back(p);
    }
    actual[0].number=256;
    for(int j=0;j<512;++j) {
        Particle p;p.position={.75,.75};p.velocity={1000.,0.};p.w=11.;
        p.inCellL=p.inBackgroundCellL=mesh.findCellLabel(p.position);actual[1].particle.push_back(p);
    }
    actual[1].number=512;expected=actual;
    std::mt19937_64 rng(86420),referenceRng(rng);std::uniform_real_distribution<double> uniform,referenceUniform;
    MCCCollisionOperator::StepStats totals;
    const auto sameParticle=[](const Particle& a,const Particle& b,bool hosts) {
        check(a.position.x==b.position.x && a.position.y==b.position.y && a.velocity.x==b.velocity.x
            && a.velocity.y==b.velocity.y && a.vz==b.vz && a.w==b.w && a.rngState==b.rngState,
            "collision stage changed particles or per-particle RNG");
        if(hosts)check(a.inCellL==b.inCellL && a.inBackgroundCellL==b.inBackgroundCellL,"spawn host policy changed");
    };
    for(int batch=0;batch<2;++batch) {
        // Start from the same charged population but retain operator/global RNG history.
        const auto counts=std::array<std::size_t,2>{actual[0].particle.size(),actual[1].particle.size()};
        if(implicit)stage.applyImplicitHalf(actual,1e-10,rng);
        else {
            stage.beginBatch();
            for(std::size_t si=0;si<actual.size();++si)
                check(stage.applyExplicit(mesh,actual[si],weights[si],1e-10,rng,batch+1),"explicit collision stage");
        }
        ParticleSpawnBuffer pending;MCCCollisionOperator::StepStats stats;
        for(std::size_t si=0;si<expected.size();++si)if(reference->operation().HasProjectile(expected[si].name())) {
            ParticleSpawnBuffer raw;
            check(reference->apply(expected[si],weights[si],1e-10,referenceRng,raw,stats,&error),"independent raw collision application");
            for(const auto& item:raw.items()) {
                if(names.count(item.species_name)==0)continue;
                auto p=item.particle;
                if(!implicit) {
                    if(p.position.x<0. || p.position.x>1. || p.position.y<0. || p.position.y>1.)continue;
                    p.inCellL=mesh.findCellLabel(p.position);p.inBackgroundCellL=mesh.findBackgroundCellLabel(p.position);
                }
                pending.Push(item.species_name,p);
            }
        }
        check(stats.n_collisions>0 && !pending.empty(),"stage fixture must generate actual secondaries");
        check(stats.per_projectile.size()==2 && stats.per_projectile.at("ion").n_collisions>0,
            "stage fixture must process multiple colliding projectiles");
        check(actual[0].particle.size()==counts[0] && actual[1].particle.size()==counts[1],"collision stage committed spawns before caller phase");
        check(rng==referenceRng && owner->rngState().applyCount==reference->rngState().applyCount,"stage changed random call sequence");
        const auto& observed=stage.batchStats();
        check(observed.n_collisions==stats.n_collisions && observed.per_channel==stats.per_channel
            && observed.spawn_per_species==stats.spawn_per_species && observed.energy_loss_J_total==stats.energy_loss_J_total,
            "collision stage changed physical statistics");
        check(stage.pending().size()==pending.size(),"unknown/outside spawned particle policy changed");
        for(std::size_t k=0;k<pending.size();++k) {
            check(stage.pending().items()[k].species_name==pending.items()[k].species_name,"spawn insertion order changed");
            sameParticle(stage.pending().items()[k].particle,pending.items()[k].particle,true);
        }
        const auto savedRng=rng;const auto savedCount=stage.pending().size();
        for(double invalid:{0.,-1.,std::numeric_limits<double>::infinity()})stage.applyImplicitHalf(actual,invalid,rng);
        check(rng==savedRng && stage.pending().size()==savedCount && stage.batchStats().n_collisions==stats.n_collisions,
            "invalid implicit half-step changed existing batch");
        const auto kept=stage.thinSpawns(management,true,.5,rng,uniform);
        const auto referenceKept=management.thin(pending,true,.5,referenceRng,referenceUniform);
        check(kept.drop>0 && kept.kept>0 && kept.kept==referenceKept.kept && kept.physKept==referenceKept.physKept
            && rng==referenceRng,"owned spawn thinning changed RNG or compensated weights");
        std::vector<unsigned char> touched(2),referenceTouched(2);
        if(implicit)management.commit(stage.pending(),mesh,actual,&touched,true);
        else {
            PicSecondaryEmission emission({},PicSecondaryEmission::Mode::Explicit);
            emission.initialize({},actual,weights,names,1.,3.);
            PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);
            PicParticleSources sources(PicParticleSources::Mode::Explicit,{},names,nullptr,mesh,embedded,nullptr);
            std::vector<Particle> see;ParticleSpawnBuffer source;
            const auto accepted=stage.acceptExplicit(management,{mesh,actual,rng,uniform,3.,{},false},emission,see,sources,source);
            touched=accepted.touchedSpecies;
            check(accepted.collision.in==static_cast<long long>(stage.pending().size()) && accepted.injectionEnergyJ==0.
                && accepted.keepProbability==1. && accepted.sourceChargeC==std::array<double,4>{},"explicit stage joins private collision output");
            check(stage.batchStats().n_collisions==stats.n_collisions && rng==referenceRng,"acceptance must not clear batch or draw disabled thinning");
        }
        management.commit(pending,mesh,expected,&referenceTouched,implicit);
        check(touched==referenceTouched,"deferred commit touched different species");
        for(std::size_t si=0;si<actual.size();++si) {
            check(actual[si].number==expected[si].number && actual[si].particle.size()==expected[si].particle.size(),"deferred commit counts changed");
            for(std::size_t k=0;k<actual[si].particle.size();++k)sameParticle(actual[si].particle[k],expected[si].particle[k],true);
        }
        totals.Accumulate(stats);stage.accumulateBatch();
        check(stage.accumulatedStats().n_collisions==totals.n_collisions,"caller-selected accumulation phase");
    }
    stage.clearAccumulated();check(stage.accumulatedStats().n_collisions==0 && !stage.pending().empty(),"clear aggregate must preserve batch");
    stage.beginBatch();check(stage.pending().empty() && stage.batchStats().n_collisions==0,"begin batch retained stale output");
    Stage disabled(nullptr,names,weights);const auto before=rng;
    disabled.applyImplicitHalf(actual,1e-10,rng);
    check(disabled.applyExplicit(mesh,actual[0],7.,1e-10,rng,1) && disabled.pending().empty()
        && disabled.batchStats().n_collisions==0 && rng==before,"disabled stage consumed RNG or published output");
}
}
namespace {
void sameAcceptedSpecies(const std::vector<ParticleType>& a,const std::vector<ParticleType>& b)
{
    check(a.size()==b.size(),"accepted species size");
    for(std::size_t si=0;si<a.size();++si) {
        check(a[si].number==b[si].number && a[si].particle.size()==b[si].particle.size(),"accepted species metadata and count");
        for(std::size_t j=0;j<a[si].particle.size();++j) {
            const auto& p=a[si].particle[j];const auto& q=b[si].particle[j];
            check(p.position.x==q.position.x && p.position.y==q.position.y && p.velocity.x==q.velocity.x && p.velocity.y==q.velocity.y && p.vz==q.vz && p.w==q.w,"accepted particle physical state");
            check(p.inCellL==q.inCellL && p.inBackgroundCellL==q.inBackgroundCellL && p.discreteVertL==q.discreteVertL && p.trackId==q.trackId && p.rngState==q.rngState && p.emitChance==q.emitChance && p.levelSetValues==q.levelSetValues,"accepted particle host/cache/identity");
            for(int k=0;k<4;++k)check(p.weight[k]==q.weight[k],"accepted interpolation cache");
        }
    }
}
void sameSpawnStatistics(const PicSpawnManagement::Stats& a,const PicSpawnManagement::Stats& b)
{
    check(a.in==b.in && a.kept==b.kept && a.drop==b.drop && a.physIn==b.physIn && a.physKept==b.physKept && a.keepProbWeightedMacro==b.keepProbWeightedMacro && a.keepProbMacroCount==b.keepProbMacroCount,"spawn acceptance statistics and probability accumulation");
}
// Original explicit driver sequence: notably SEE is thinned before source,
// whereas source is committed before SEE, and all use one initial probability.
PicSpawnManagement::ExplicitResult referenceExplicitAcceptance(const PicSpawnManagement& management,
    const ParticleManagementConfig& config,const PicSpawnManagement::ExplicitStage& stage,
    ParticleSpawnBuffer& collision,PicSecondaryEmission& emission,std::vector<Particle>& secondary,
    PicParticleSources& sources,ParticleSpawnBuffer& source,const std::vector<double>& weights,
    const std::unordered_map<std::string,std::size_t>& names)
{
    PicSpawnManagement::ExplicitResult result;
    const long long count=static_cast<long long>(collision.size())+static_cast<long long>(secondary.size())+static_cast<long long>(source.size());
    const bool thin=stage.managementActive && config.thinningEnabled;
    result.keepProbability=management.keepProbability(stage.species,count,thin);
    result.collision=management.thin(collision,thin && config.thinningMccSpawn,result.keepProbability,stage.rng,stage.uniform);
    result.secondary=emission.thinExplicit(management,secondary,thin && config.thinningSeeSpawn,result.keepProbability,stage.rng,stage.uniform);
    result.source=sources.thinExplicit(management,source,thin && config.thinningMccSpawn,result.keepProbability,stage.rng,stage.uniform);
    emission.accountRetained(secondary,stage.species,weights);
    const auto injection=sources.accountRetained(source,stage.species,weights,names,stage.velocityScale,stage.domainBoundaryIds);
    result.injectionEnergyJ=emission.stats().emittedEnergyJ+injection.energyJ;result.sourceChargeC=injection.domainChargeC;
    result.touchedSpecies.assign(stage.species.size(),0u);
    management.commit(collision,stage.mesh,stage.species,&result.touchedSpecies);
    management.commit(source,stage.mesh,stage.species,&result.touchedSpecies);
    if(emission.enabled() && emission.emittedIndex()!=static_cast<std::size_t>(-1) && !secondary.empty())
        management.commitExplicitSee(secondary,emission.emittedIndex(),stage.species,result.touchedSpecies);
    return result;
}
void explicitSpawnAcceptance(int flags,bool active,int groups,int failure=0)
{
    using Manager=PicSpawnManagement;using Emission=PicSecondaryEmission;using Sources=PicParticleSources;
    using Side=ExplicitPIC::Pusher::DomainSide;
    DictionaryMesh input;input.Nx=input.Ny=4;input.maxCellNum=128;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    DictionaryGeometry dictionary;Line line;line.points[0]={0.,1.};line.points[1]={0.,0.};dictionary.lines={line};
    for(int id=1;id<=5;++id){DictionaryGeometry::BoundaryEntity b;b.id=id;b.lineIds={0};dictionary.boundaryById[id]=b;}
    PicEmbeddedGeometry embedded(PicEmbeddedGeometry::Mode::Explicit);
    std::vector<ParticleType> actual{{"electron",PhysicalConstants::eV_J,-2.},{"ion",2.*PhysicalConstants::eV_J,1.},{"untouched",1.,0.}};
    Particle particle;particle.position={.25,.25};particle.velocity={3.,4.};particle.vz=12.;particle.w=0.;
    particle.inCellL=particle.inBackgroundCellL=mesh.findCellLabel(particle.position);
    actual[0].particle.push_back(particle);actual[1].particle.push_back(particle);
    for(auto& sp:actual)sp.number=77;
    auto expected=actual;const auto initial=actual;
    const std::vector<double> weights{2.}; // Missing ion default retains the legacy fallback.
    const std::unordered_map<std::string,std::size_t> names{{"electron",0},{"ion",1},{"untouched",2}};
    ParticleManagementConfig config;config.thinningEnabled=(flags&1)!=0;config.thinningMccSpawn=(flags&2)!=0;config.thinningSeeSpawn=(flags&4)!=0;
    config.softCapTotal=12;config.hardCapTotal=9;config.minKeepProb=.2;
    Manager manager(failure==3?Manager::Mode::Implicit:Manager::Mode::Explicit,config,names,weights);
    SEEConfig seeConfig;seeConfig.enabled=(groups&2)!=0;seeConfig.gamma=2.;seeConfig.emittedSpecies="electron";
    seeConfig.cathodeSides={"left","right","bottom","top"};seeConfig.emittedEnergyModel="mono";seeConfig.emittedEnergy_eV=5.;
    seeConfig.angular="fixed";seeConfig.injectEps_m=.06;seeConfig.maxEmitMacro=3;
    Emission emission(seeConfig,Emission::Mode::Explicit),referenceEmission(seeConfig,Emission::Mode::Explicit);
    emission.initialize({},actual,weights,names,2.,2.);referenceEmission.initialize({},expected,weights,names,2.,2.);
    std::vector<SpeciesConfig> emitters;
    for(int id=1;id<=5;++id){SpeciesConfig s;s.name="electron";s.q=-2.;s.m=PhysicalConstants::eV_J;s.weight=2.;s.drift={2.,0.};
        s.source.enabled=(groups&4)!=0;s.source.boundaryId=id;s.source.current_A=16.;s.source.maxEmitMacro=4;emitters.push_back(s);}
    Sources sources(Sources::Mode::Explicit,emitters,names,&dictionary,mesh,embedded,nullptr);
    Sources referenceSources(Sources::Mode::Explicit,emitters,names,&dictionary,mesh,embedded,nullptr);
    ParticleSpawnBuffer collision,referenceCollision,source,referenceSource;std::vector<Particle> secondary,referenceSecondary;
    if(groups&1)for(int j=0;j<21;++j){auto p=particle;p.trackId=100+j;p.position.x=j%2?-.25:.25;
        p.inCellL=-1;p.w=j%3?0.:3.;collision.Push(j%3==0?"missing":j%3==1?"electron":"ion",p);}
    referenceCollision=collision;
    std::mt19937_64 rng(833),referenceRng(rng);std::uniform_real_distribution<double> uniform,referenceUniform;
    ExplicitPIC::Pusher::DomainImpactBuffer impacts;
    for(auto side:{Side::kLeft,Side::kRight,Side::kBottom,Side::kTop})impacts.Push(particle,side);
    const Emission::BoundaryGeometry geometry{mesh,nullptr,dictionary,embedded};
    const Emission::DomainBox box{{0.,0.},{1.,1.},{1e-6,1e-6}};
    const auto prepare=[&](auto& model,auto& emit,auto& see,auto& src,auto& engine,auto& draw){
        model.beginStep();model.processExplicit(geometry,initial,1,6.,impacts,nullptr,box,{},see,engine,draw);
        emit.emit({0,1.,1.,1.},mesh,embedded,nullptr,engine,draw,src);
        for(std::size_t j=0;j<see.size();++j)see[j].trackId=2000+j;
        for(std::size_t j=0;j<src.size();++j)src.items_mut()[j].particle.trackId=1000+j;
    };
    prepare(emission,sources,secondary,source,rng,uniform);prepare(referenceEmission,referenceSources,referenceSecondary,referenceSource,referenceRng,referenceUniform);
    check(secondary.size()==std::size_t(groups&2?12:0) && source.size()==std::size_t(groups&4?20:0),"explicit acceptance fixture produces domain SEE and five tagged sources");
    if(failure==1){secondary.push_back(particle);referenceSecondary.push_back(particle);}
    if(failure==2){source.Push("electron",particle);referenceSource.Push("electron",particle);}
    const auto beforeRng=rng;Manager::ExplicitResult result,reference;
    const Manager::ExplicitStage stage{mesh,actual,rng,uniform,2.,{1,2,3,4},active};
    const Manager::ExplicitStage ref{mesh,expected,referenceRng,referenceUniform,2.,{1,2,3,4},active};
    std::string error,referenceError;
    try{result=manager.acceptExplicit(stage,collision,emission,secondary,sources,source);}catch(const std::logic_error& e){error=e.what();}
    if(failure!=3)try{reference=referenceExplicitAcceptance(manager,config,ref,referenceCollision,referenceEmission,referenceSecondary,referenceSources,referenceSource,weights,names);}catch(const std::logic_error& e){referenceError=e.what();}
    if(failure==3){check(!error.empty() && rng==beforeRng,"wrong-mode explicit stage rejects before mutation");}
    else check(error==referenceError && (failure==0)==error.empty(),"explicit acceptance preserves failure point and diagnostic");
    check(rng==referenceRng,"explicit acceptance changed RNG consumption");sameAcceptedSpecies(actual,expected);
    const auto sameParticles=[](const auto& a,const auto& b){std::vector<ParticleType> x{{"p",1.,1.}},y=x;x[0].particle=a;y[0].particle=b;sameAcceptedSpecies(x,y);};
    const auto sameBuffer=[&](const auto& a,const auto& b){check(a.size()==b.size(),"retained named spawn count");std::vector<Particle> x,y;
        for(std::size_t j=0;j<a.size();++j){check(a.items()[j].species_name==b.items()[j].species_name,"retained spawn species order");x.push_back(a.items()[j].particle);y.push_back(b.items()[j].particle);}sameParticles(x,y);};
    sameBuffer(collision,referenceCollision);sameBuffer(source,referenceSource);sameParticles(secondary,referenceSecondary);
    const auto& a=emission.stats();const auto& b=referenceEmission.stats();
    check(a.ionImpactMacroCathode==b.ionImpactMacroCathode && a.ionImpactPhysCathode==b.ionImpactPhysCathode && a.ionImpactEnergyJCathode==b.ionImpactEnergyJCathode
        && a.emittedMacro==b.emittedMacro && a.emittedPhys==b.emittedPhys && a.emittedEnergyJ==b.emittedEnergyJ && a.emittedEnergyJ2==b.emittedEnergyJ2
        && a.injectedCharge==b.injectedCharge,"retained SEE charge/energy/impact accounting changed");
    if(failure){sameAcceptedSpecies(actual,initial);return;}
    sameSpawnStatistics(result.collision,reference.collision);sameSpawnStatistics(result.secondary,reference.secondary);sameSpawnStatistics(result.source,reference.source);
    check(result.keepProbability==reference.keepProbability && result.injectionEnergyJ==reference.injectionEnergyJ && result.sourceChargeC==reference.sourceChargeC
        && result.touchedSpecies==reference.touchedSpecies,"explicit shared probability, accounting and touched species");
    check(actual[2].number==77 && result.touchedSpecies[2]==0,"untouched metadata stays untouched");
    if(flags==7 && active && groups==7)check(result.collision.drop>0 && result.secondary.drop>0 && result.source.drop>0
        && result.collision.kept>0 && result.secondary.kept>0 && result.source.kept>0,"joint fixture exercises all three thinning groups");
    // Retained reports and tags stay valid after commit, before external culling.
    const auto again=sources.accountRetained(source,actual,weights,names,2.,{1,2,3,4});
    check(again.domainChargeC==result.sourceChargeC,"source provenance survives acceptance");
    emission.accountRetained(secondary,actual,weights);check(emission.stats().injectedCharge==b.injectedCharge,"SEE provenance survives acceptance");
}
void sameSpawnBalance(const PicStepParticleBalance& a,const PicStepParticleBalance& b)
{
    check(a.current().macroBeginC==b.current().macroBeginC && a.current().sourceC==b.current().sourceC && a.current().secondaryC==b.current().secondaryC && a.injectionEnergyJ()==b.injectionEnergyJ(),"spawn charge/energy accounting");
    for(auto kind:{PicStepParticleBalance::Spawn::Collision,PicStepParticleBalance::Spawn::Secondary,PicStepParticleBalance::Spawn::Source}) {
        const auto& x=a.spawns(kind);const auto& y=b.spawns(kind);
        check(x.energyJ==y.energyJ && x.species.size()==y.species.size(),"spawn energy and species scope");
        for(std::size_t si=0;si<x.species.size();++si)check(x.species[si].macro==y.species[si].macro && x.species[si].physical==y.species[si].physical,"spawn physical/macro accounting");
    }
}
// Original driver operations, kept as an independent sequence for new stage tests.
void referenceImplicitAcceptance(const PicSpawnManagement& management,const ParticleManagementConfig& config,
    const PicSpawnManagement::ImplicitStage& stage,ParticleSpawnBuffer& pending,PicSpawnManagement::ImplicitKind kind,
    const std::vector<double>& weights,const std::unordered_map<std::string,std::size_t>& names,PicSpawnManagement::Stats* statistics)
{
    using Kind=PicSpawnManagement::ImplicitKind;using Spawn=PicStepParticleBalance::Spawn;
    if(kind==Kind::Secondary && pending.empty())return;
    const bool thinning=stage.managementActive && config.thinningEnabled && (kind==Kind::Secondary?config.thinningSeeSpawn:config.thinningMccSpawn);
    const double probability=management.keepProbability(stage.species,static_cast<long long>(pending.size()),thinning);
    const auto local=management.thin(pending,thinning,probability,stage.rng,stage.uniform);
    if(statistics){if(kind==Kind::Collision)statistics->accumulate(local);else *statistics=local;}
    switch(kind) {
    case Kind::Collision:
        stage.balance.addSpawns(Spawn::Collision,pending,stage.species,weights,names,stage.velocityScale);
        management.commit(pending,stage.mesh,stage.species,nullptr,true);break;
    case Kind::Secondary:
        stage.balance.addSpawns(Spawn::Secondary,pending,stage.species,weights,names,stage.velocityScale);
        management.commit(pending,stage.mesh,stage.species);
        stage.balance.includeSecondarySource();break;
    case Kind::Source:
        management.commit(pending,stage.mesh,stage.species);
        stage.balance.addSpawns(Spawn::Source,pending,stage.species,weights,names,stage.velocityScale);break;
    }
}
void implicitSpawnAcceptance()
{
    using Manager=PicSpawnManagement;using Kind=Manager::ImplicitKind;
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=2;dictionary.maxCellNum=16;dictionary.maxInitialIterations=0;
    QuadtreeMesh mesh(dictionary);const std::unordered_map<std::string,std::size_t> names{{"electron",0},{"ion",1}};
    const std::vector<double> weights{7.};std::vector<ParticleType> seed;
    seed.emplace_back("electron",electron().m,electron().q);seed.emplace_back("ion",6.6335209e-26,-electron().q);
    Particle particle;particle.position={.25,.25};particle.velocity={2.,3.};particle.vz=4.;particle.w=7.;
    particle.inCellL=mesh.findCellLabel(particle.position);particle.inBackgroundCellL=mesh.findBackgroundCellLabel(particle.position);
    for(auto& sp:seed){sp.particle.push_back(particle);sp.number=77;}
    for(auto kind:{Kind::Collision,Kind::Secondary,Kind::Source})for(bool active:{false,true})for(bool empty:{false,true})for(int flags=0;flags<8;++flags) {
        ParticleManagementConfig config;config.enabled=true;config.softCapTotal=8;config.minKeepProb=.2;
        config.thinningEnabled=(flags&1)!=0;config.thinningMccSpawn=(flags&2)!=0;config.thinningSeeSpawn=(flags&4)!=0;
        Manager management(Manager::Mode::Implicit,config,names,weights);
        auto actual=seed,reference=seed;ParticleSpawnBuffer pending,referencePending;
        if(!empty)for(int j=0;j<25;++j) {
            auto p=particle;p.position=j%4==0?Vector2D(2.,.5):Vector2D(.25,.25);p.inCellL=-1;p.inBackgroundCellL=-1;
            p.w=j%2?0.:3.;p.trackId=j+1;p.rngState=100+j;p.levelSetValues={.2,.3};p.weight[2]=.5;
            const std::string name=j%3==0?"electron":j%3==1?"ion":"unknown";
            pending.Push(name,p);referencePending.Push(name,p);
        }
        PicStepParticleBalance balance(actual,weights),referenceBalance(reference,weights);
        std::mt19937_64 rng(971),referenceRng(971);std::uniform_real_distribution<double> uniform,referenceUniform;
        Manager::Stats statistics{3,2,1,5.,4.,1.,2},referenceStatistics=statistics;
        Manager::ImplicitStage stage{mesh,actual,rng,uniform,balance,3.,active},ref{mesh,reference,referenceRng,referenceUniform,referenceBalance,3.,active};
        struct Capture {std::ostringstream text;std::streambuf* old=std::cerr.rdbuf(text.rdbuf());~Capture(){std::cerr.rdbuf(old);}} capture;
        management.acceptImplicit(stage,pending,kind,&statistics);const auto actualWarnings=capture.text.str();capture.text.str("");
        referenceImplicitAcceptance(management,config,ref,referencePending,kind,weights,names,&referenceStatistics);
        check(capture.text.str()==actualWarnings,"spawn commit warning sequence");
        sameAcceptedSpecies(actual,reference);sameSpawnStatistics(statistics,referenceStatistics);sameSpawnBalance(balance,referenceBalance);
        check(rng==referenceRng && pending.size()==referencePending.size(),"spawn thinning RNG and retained buffer size");
        for(std::size_t j=0;j<pending.size();++j) {
            check(pending.items()[j].species_name==referencePending.items()[j].species_name,"retained buffer order");
            check(pending.items()[j].particle.w==referencePending.items()[j].particle.w && pending.items()[j].particle.trackId==referencePending.items()[j].particle.trackId && pending.items()[j].particle.inCellL==referencePending.items()[j].particle.inCellL,"retained buffer is not repaired in place or cleared");
        }
        if(empty)check(actual[0].number==77,"empty retained buffers never normalize live metadata");
        auto silent=seed;PicStepParticleBalance silentBalance(silent,weights);std::mt19937_64 silentRng(971);std::uniform_real_distribution<double> silentUniform;
        ParticleSpawnBuffer silentPending;
        // Disabled statistics must also work for an empty buffer in all three phases.
        management.acceptImplicit({mesh,silent,silentRng,silentUniform,silentBalance,3.,active},silentPending,kind);
    }
    // A deliberately undersized balance exposes accounting vs commit order:
    // source has committed when accounting throws, collision/SEE have not.
    for(auto kind:{Kind::Collision,Kind::Secondary,Kind::Source}) {
        ParticleManagementConfig config;Manager management(Manager::Mode::Implicit,config,names,weights);
        auto actual=seed;ParticleSpawnBuffer pending;pending.Push("electron",particle);
        PicStepParticleBalance balance(0.,0);std::mt19937_64 rng(2);std::uniform_real_distribution<double> uniform;Manager::Stats statistics;
        bool thrown=false;try{management.acceptImplicit({mesh,actual,rng,uniform,balance,3.,false},pending,kind,&statistics);}catch(const std::out_of_range&){thrown=true;}
        check(thrown && actual[0].particle.size()==(kind==Kind::Source?2:1) && statistics.in==1 && statistics.kept==1 && pending.size()==1,"failure retains original accounting/commit partial state");
    }
    ParticleManagementConfig config;Manager wrong(Manager::Mode::Explicit,config,names,weights);
    auto actual=seed;PicStepParticleBalance balance(actual,weights);ParticleSpawnBuffer pending;pending.Push("electron",particle);
    std::mt19937_64 rng(2),before(rng);std::uniform_real_distribution<double> uniform;bool rejected=false;
    try{wrong.acceptImplicit({mesh,actual,rng,uniform,balance,3.,true},pending,Kind::Source);}catch(const std::logic_error&){rejected=true;}
    check(rejected && rng==before && pending.size()==1,"wrong management mode rejected before RNG/buffer mutation");sameAcceptedSpecies(actual,seed);
}

void completeImplicitHalfStep(CollisionsRngMode mode,bool active,bool recordStatistics)
{
    using Manager=PicSpawnManagement;using Stage=PicCollisions::Stage;
    auto ion=electron();ion.name="ion";ion.m=6.6335209e-26;ion.q=-ion.q;ion.weight=11.;
    const auto collisionInputs=inputs(mode,true);std::string error;
    auto owner=PicCollisions::create(collisionInputs,{electron(),ion},3.,error),referenceOwner=PicCollisions::create(collisionInputs,{electron(),ion},3.,error);
    check(owner && referenceOwner,"complete collision half-step fixture");
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=2;dictionary.maxCellNum=16;dictionary.maxInitialIterations=0;QuadtreeMesh mesh(dictionary);
    const std::unordered_map<std::string,std::size_t> names{{"electron",0},{"ion",1}};const std::vector<double> weights{7.,11.};
    Stage stage(owner.get(),names,weights),referenceStage(referenceOwner.get(),names,weights);
    std::vector<ParticleType> actual;actual.emplace_back("electron",electron().m,electron().q);actual.emplace_back("ion",ion.m,ion.q);
    for(int si=0;si<2;++si)for(int j=0;j<(si?512:256);++j) {
        Particle p;p.position=si?Vector2D(.75,.75):Vector2D(.25,.25);p.velocity={si?1000.:2000000.,0.};p.vz=si?0.:100000.;p.w=weights[si];
        p.inCellL=mesh.findCellLabel(p.position);p.inBackgroundCellL=mesh.findBackgroundCellLabel(p.position);p.trackId=1000*si+j;actual[si].particle.push_back(p);
    }
    for(auto& sp:actual)sp.number=static_cast<int>(sp.particle.size());auto reference=actual;
    ParticleManagementConfig config;config.thinningEnabled=true;config.thinningMccSpawn=true;config.softCapTotal=500;config.minKeepProb=.2;
    Manager management(Manager::Mode::Implicit,config,names,weights);
    PicStepParticleBalance balance(actual,weights),referenceBalance(reference,weights);
    std::mt19937_64 rng(86420),referenceRng(86420);std::uniform_real_distribution<double> uniform,referenceUniform;
    Manager::Stats statistics,referenceStatistics;
    const Manager::ImplicitStage context{mesh,actual,rng,uniform,balance,3.,active},ref{mesh,reference,referenceRng,referenceUniform,referenceBalance,3.,active};
    for(int half=0;half<2;++half) {
        check(stage.advanceImplicitHalf(management,context,1e-10,recordStatistics?&statistics:nullptr),"complete collision half-step applied");
        check(referenceStage.applyImplicitHalf(reference,1e-10,referenceRng),"reference kernel half-step applied");
        ParticleSpawnBuffer retained;for(const auto& item:referenceStage.pending().items())retained.Push(item.species_name,item.particle);
        referenceImplicitAcceptance(management,config,ref,retained,Manager::ImplicitKind::Collision,weights,names,recordStatistics?&referenceStatistics:nullptr);
        sameAcceptedSpecies(actual,reference);sameSpawnStatistics(statistics,referenceStatistics);sameSpawnBalance(balance,referenceBalance);
        check(rng==referenceRng && owner->rngState().applyCount==referenceOwner->rngState().applyCount,"complete half-step preserves RNG and collision counter");
        check(stage.batchStats().n_collisions>0 && stage.batchStats().n_collisions==referenceStage.batchStats().n_collisions && stage.batchStats().per_channel==referenceStage.batchStats().per_channel,"actual multi-projectile collision statistics");
        check(stage.pending().size()==retained.size() && !retained.empty(),"half-step retains spawned buffer after commit");
        check(stage.accumulatedStats().n_collisions==referenceStage.accumulatedStats().n_collisions,"complete half-step must not publish batch totals early");
        stage.accumulateBatch();referenceStage.accumulateBatch();
        const auto savedRng=rng;const auto savedCounter=owner->rngState();const auto saved=statistics;const auto count=stage.pending().size();
        for(double dt:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
            check(!stage.advanceImplicitHalf(management,context,dt,&statistics),"invalid half-step is inactive");
        Stage disabled(nullptr,names,weights);check(!disabled.advanceImplicitHalf(management,context,1e-10,&statistics),"disabled half-step is inactive");
        check(rng==savedRng && owner->rngState().applyCount==savedCounter.applyCount && stage.pending().size()==count,"inactive half-step preserves previous work");sameSpawnStatistics(statistics,saved);
    }
    if(active && recordStatistics)check(statistics.drop>0 && statistics.kept>0,"complete half-step exercises probabilistic thinning");
}

void collisionReports()
{
    auto data=inputs(CollisionsRngMode::kGlobal,true);
    auto reactions=std::make_shared<ReactionDatabase>(*data.reactions);
    auto set=*reactions->Get("ion-Ar");set.channels.clear();
    for(const auto& name:{"ELASTIC_1","elastic_2","ElAsTiC_3","bounce","swap","disabled"}) {
        ReactionChannel channel;channel.name=name;
        channel.cross_section_key=std::string(name)=="swap"?"Ar+-Ar:CHARGE_EXCHANGE":"Ar+-Ar:ELASTIC";
        channel.enabled=std::string(name)!="disabled";set.channels.push_back(channel);
    }
    std::string error;check(reactions->Register(set,true,&error),"report reactions");data.reactions=reactions;
    auto e=electron();e.m=2.*PhysicalConstants::eV_J;e.numParticles=0;
    auto ion=e;ion.name="ion";ion.q=-e.q;
    auto neutral=e;neutral.name="neutral";neutral.q=0.;neutral.weight=3.;
    auto service=PicCollisions::create(data,{e,ion,neutral},2.,error);
    check(service!=nullptr,error.c_str());const auto& op=service->operation();
    for(const auto& name:{"ELASTIC_1","elastic_2","ElAsTiC_3","bounce"})
        check(op.ChannelKind("ion",name)=="elastic","model kind independent of literal report spelling");
    check(op.ChannelKind("ion","swap")=="charge_exchange" && op.ChannelKind("ion","disabled").empty()
        && op.ChannelKind("ion","missing").empty() && op.ChannelKind("missing","elastic").empty(),
        "only enabled initialized channels are classified");
    ReactionChannel heuristic;heuristic.products={"electron","Ar+"};
    check(PostCollisionModel::InferKind(set,heuristic)=="ionization","legacy product heuristic");
    heuristic.products.clear();heuristic.name="CX";
    check(PostCollisionModel::InferKind(set,heuristic)=="charge_exchange","legacy CX alias");

    DictionaryMesh dict;dict.Nx=dict.Ny=1;dict.maxCellNum=16;
    QuadtreeMesh mesh(dict);PicSpeciesState particles(3,17);
    particles.initialize({e,ion,neutral},NormalizationScales{},mesh,true);
    particles.rebuildIndex();
    Particle p;p.velocity={1.,0.};p.vz=100.;p.w=2.;particles.species()[0].particle.push_back(p);
    p.velocity.x=2.;p.w=-1.;particles.species()[0].particle.push_back(p);
    p.velocity.x=std::numeric_limits<double>::quiet_NaN();particles.species()[0].particle.push_back(p);
    auto energy=summarizePicCollisionElectronEnergy(particles.species(),particles.weights(),2.);
    check(energy.macroCount==2 && energy.physicalCount==9. && std::abs(energy.meanEV-120./9.)<1e-12
        && energy.tailFraction==7./9.,"planar weighted electron proxy/invalid-velocity filtering");
    auto invalid=particles.species();invalid[0].particle.resize(1);invalid[0].particle[0].w=-1.;
    energy=summarizePicCollisionElectronEnergy(invalid,{std::numeric_limits<double>::infinity()},2.);
    check(energy.macroCount==1 && energy.physicalCount==0. && energy.meanEV==0.,"macro count precedes invalid weight filter");
    p.velocity={1.,2.};p.vz=2.;p.w=2.;particles.species()[2].particle.push_back(p);
    p.velocity={0.,0.};p.vz=0.;p.w=0.;particles.species()[2].particle.push_back(p);
    const auto neutralEnergy=summarizePicNeutralEnergy(particles.species()[2],3.,2.);
    check(neutralEnergy.macroCount==2 && neutralEnergy.physicalCount==5.
        && std::abs(neutralEnergy.meanEV-72./5.)<1e-12,"three-component neutral energy/weight fallback");

    MCCCollisionOperator::StepStats stats;
    auto& ions=stats.per_projectile["ion"];ions.per_channel={{"ELASTIC_1",2},{"elastic_2",3},{"ElAsTiC_3",5},{"bounce",7},{"swap",11}};
    ions.n_collisions=28;ions.energy_loss_J_total=.125;ions.max_nu_dt_majorant=.25;
    stats.per_projectile["electron"]=ions;stats.spawn_per_species["neutral"]=11;
    const auto path=std::filesystem::current_path()/("collision_reports_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(path);
    auto text=[&](const char* name){std::ifstream in(path/name);return std::string(std::istreambuf_iterator<char>(in),{});};
    auto rows=[&](const char* name){std::istringstream in(text(name));std::vector<std::string> r;std::string line;
        while(std::getline(in,line))if(!line.empty() && line[0]!='#')r.push_back(line);return r;};
    const auto rng=particles.rng();const auto callCount=service->rngState().applyCount;
    {
        ExplicitPicCollisionReports writer(path,2);
        writer.writeStep(1,.5,particles,2.,&op,stats);
        writer.writeStep(2,1.,particles,2.,&op,stats);
        writer.writeStep(4,2.,particles,2.,nullptr,stats);
        writer.append(2,1.,particles,2.,stats,stats);
        writer.append(4,2.,particles,2.,stats,stats);
    }
    const auto ionsOut=rows("ion_report.txt"),neutralsOut=rows("neutral_report.txt");
    check(ionsOut==std::vector<std::string>{"2,1,ion,17,11,28,0.125,0.25"},"ion cadence, summed kinds, positive-charge filter and MCC-off gating");
    check(neutralsOut.size()==2 && neutralsOut[0].rfind("2,1,neutral,2,5,",0)==0
        && neutralsOut[0].substr(neutralsOut[0].find_last_of(',')+1)=="11"
        && neutralsOut[1].substr(neutralsOut[1].find_last_of(',')+1)=="0","neutral report remains active without MCC");
    const auto report=text("collisions_report.txt");
    check(report.find("=== MCC report step=2 t=1 ===")!=std::string::npos
        && report.find("=== MCC report step=4 t=2 ===")!=std::string::npos
        && report.find("electron_n_macro=2")!=std::string::npos,"append reports and electron reduction");
    for(int every:{0,-1}) {
        ExplicitPicCollisionReports writer(path,every);writer.writeStep(2,1.,particles,2.,&op,stats);
    }
    check(rows("ion_report.txt").empty() && rows("neutral_report.txt").empty()
        && text("collisions_report.txt")==report,"disabled cadence and legacy startup truncation");
    const auto blocked=path/"blocked";std::filesystem::create_directories(blocked/"ion_report.txt");
    std::filesystem::create_directories(blocked/"neutral_report.txt");
    std::filesystem::create_directories(blocked/"collisions_report.txt");
    ExplicitPicCollisionReports failed(blocked,1);failed.writeStep(1,.5,particles,2.,&op,stats);
    failed.append(1,.5,particles,2.,stats,stats);
    check(particles.rng()==rng && service->rngState().applyCount==callCount && ions.n_collisions==28
        && particles.species()[0].particle[0].vz==100.,"diagnostics mutate no RNG, physics or statistics");
}
}
int main()
{
    static_assert(!std::is_move_constructible_v<PicCollisions>);
    static_assert(!std::is_copy_constructible_v<PicCollisions>);
    static_assert(std::is_same_v<decltype(std::declval<PicCollisions&>().operation()), const MCCCollisionOperator&>);
    try
    {
        rejectedInitialization();
        collisionReports();
        implicitSpawnAcceptance();
        for(int flags=0;flags<8;++flags)for(bool active:{false,true})for(int groups=0;groups<8;++groups)
            explicitSpawnAcceptance(flags,active,groups);
        for(int failure:{1,2,3})explicitSpawnAcceptance(7,true,7,failure);
        for(auto mode:{CollisionsRngMode::kGlobal,CollisionsRngMode::kPerParticle,CollisionsRngMode::kCounterBased})
            for(bool active:{false,true})for(bool statistics:{false,true})completeImplicitHalfStep(mode,active,statistics);
        for (auto mode : {CollisionsRngMode::kGlobal, CollisionsRngMode::kPerParticle, CollisionsRngMode::kCounterBased})
            for (double scale : {1.0, 3.0, 8.0}) compareAdvance(mode, scale);
        for(auto mode:{CollisionsRngMode::kGlobal,CollisionsRngMode::kPerParticle,CollisionsRngMode::kCounterBased})
            for(bool implicit:{false,true})for(bool unknown:{false,true})collisionStage(mode,implicit,unknown);
        std::cout << "PIC collision ownership, sequential RNG and unit conversion passed\n";
        return 0;
    }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
