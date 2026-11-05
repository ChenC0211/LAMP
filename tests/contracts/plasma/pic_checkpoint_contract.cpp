#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Modules/Plasma/IO/PicCheckpointWriter.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/IO/ImplicitPicFieldCheckpoint.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/MeshCacheIO.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/IO/FileManagement.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include "Modules/Plasma/PicStepSession.h"
#include "Modules/Plasma/PicCaseInput.h"
#include "Modules/Plasma/Physics/CollisionsConfig.h"
#include "Modules/Plasma/Physics/PicCollisions.h"
#include "Modules/Plasma/Physics/BuiltInCrossSections.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Core/Units/NormalizationContext.h"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <array>
#include <cstddef>
#include <new>

namespace {
void check(bool value,const char* why) { if(!value) throw std::runtime_error(why); }
std::string bytes(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary);return {std::istreambuf_iterator<char>(in),{}};
}
// Selection fixtures need only nonempty regular files, not serialized payloads.
void touch(const std::filesystem::path& path) { std::ofstream file(path); file << 'x'; }
template<class Action> void rejects(Action action, const char* why) {
    bool failed=false; try { action(); } catch(const std::exception&) { failed=true; }
    check(failed,why);
}
std::array<DiscCoeff*,8> coefficients(Vert& vertex) {
    return {&vertex.lapCoeff,&vertex.selfCoeff,&vertex.nCoeff,&vertex.tCoeff,
            &vertex.xCoeff,&vertex.yCoeff,&vertex.nnCoeff,&vertex.ttCoeff};
}
bool coefficientMember(std::size_t offset) {
    return offset<sizeof(int) ||
        (offset>=offsetof(DiscCoeff,discOpr) && offset<offsetof(DiscCoeff,discOpr)+6*sizeof(double)) ||
        (offset>=offsetof(DiscCoeff,rhs) && offset<offsetof(DiscCoeff,rhs)+sizeof(double));
}
void seedMeshCoefficients(BackgroundMesh& mesh) {
    for(int vi=0;vi<mesh.vertLabelCapacity();++vi) {
        auto values=coefficients(mesh.verts(vi));
        for(std::size_t ci=0;ci<values.size();++ci) {
            auto& value=*values[ci];value.discNum=(vi+static_cast<int>(ci))%7;
            for(int j=0;j<6;++j)value.discOpr[j]=vi+ci*.125+j*.25;
            value.rhs=vi+ci*.5+1.25; // Includes a meaningful nonzero RHS when count is zero.
            auto* representation=reinterpret_cast<unsigned char*>(&value);
            for(std::size_t offset=0;offset<sizeof(value);++offset)
                if(!coefficientMember(offset))representation[offset]=0xA5;
        }
    }
}
void verifyMeshCoefficients(BackgroundMesh& mesh,bool zeroPadding=false) {
    for(int vi=0;vi<mesh.vertLabelCapacity();++vi) {
        auto values=coefficients(mesh.verts(vi));
        for(std::size_t ci=0;ci<values.size();++ci) {
            const auto& value=*values[ci];check(value.discNum==(vi+static_cast<int>(ci))%7,"legacy coefficient count");
            for(int j=0;j<6;++j)check(value.discOpr[j]==vi+ci*.125+j*.25,"legacy operator coefficient");
            check(value.rhs==vi+ci*.5+1.25,"legacy RHS including zero-count stencil");
            if(zeroPadding)for(std::size_t offset=0;offset<sizeof(value);++offset)
                if(!coefficientMember(offset))check(reinterpret_cast<const unsigned char*>(&value)[offset]==0,
                    "mesh writer retained native coefficient padding");
        }
    }
}
void meshFixtures(const std::filesystem::path& output, const std::filesystem::path& input={},bool zeroPadding=false) {
    namespace fs=std::filesystem;
    const auto destination=fs::absolute(output);
    const auto source=input.empty()?fs::path{}:fs::absolute(input);
    DictionaryMesh dictionary;dictionary.Nx=dictionary.Ny=2;dictionary.maxCellNum=16;
    dictionary.maxInitialIterations=0;dictionary.levelSetCount=2;
    QuadtreeMesh quadtree(dictionary);
    if(source.empty())seedMeshCoefficients(quadtree);
    else { quadtree.readCache((source/"qmc").generic_string()+"/");verifyMeshCoefficients(quadtree,zeroPadding); }
    fs::create_directories(destination/"qmc");quadtree.writeCache((destination/"qmc").generic_string()+"/");
    BackgroundMesh background(dictionary);
    if(source.empty())seedMeshCoefficients(background);
    else {
        struct Restore { fs::path saved=fs::current_path();~Restore(){fs::current_path(saved);} } restore;
        fs::current_path(source/"bmc");background.readCache();verifyMeshCoefficients(background,zeroPadding);
    }
    fs::create_directories(destination/"bmc/cache");background.writeCache((destination/"bmc/cache").generic_string()+"/");
}
void checkCoefficientSerialization(const std::filesystem::path& root) {
    for(unsigned char fill:{0x5A,0xCC}) {
        alignas(DiscCoeff) std::array<unsigned char,sizeof(DiscCoeff)> storage;storage.fill(fill);
        auto* value=::new(storage.data()) DiscCoeff;
        check(value->discNum==0 && value->rhs==0.,"default stencil retained uninitialized RHS");
        for(double coefficient:value->discOpr)check(coefficient==0.,"default stencil coefficient");
        value->~DiscCoeff();
    }
    DiscCoeff value;value.discNum=0;
    for(int j=0;j<6;++j)value.discOpr[j]=j+.125;
    const std::uint64_t negativeZero=0x8000000000000000ULL,nanPayload=0x7FF8000000001234ULL;
    static_assert(sizeof(double)==sizeof(negativeZero),"native cache double width");
    std::memcpy(&value.discOpr[0],&negativeZero,sizeof(double));
    std::memcpy(&value.rhs,&nanPayload,sizeof(double));
    auto* representation=reinterpret_cast<unsigned char*>(&value);
    for(std::size_t offset=0;offset<sizeof(value);++offset)
        if(!coefficientMember(offset))representation[offset]=0xA5;
    const std::string original(reinterpret_cast<const char*>(&value),sizeof(value));
    std::ostringstream output(std::ios::binary);MeshCacheIO::writeDiscCoeff(output,value);
    const auto encoded=output.str();check(encoded.size()==sizeof(value),"native record width changed");
    for(std::size_t offset=0;offset<encoded.size();++offset)
        check(encoded[offset]==(coefficientMember(offset)?original[offset]:0),"member bits or zero padding changed");
    check(original==std::string(reinterpret_cast<const char*>(&value),sizeof(value)),"writer mutated live coefficient");
    DiscCoeff restored;std::istringstream input(encoded,std::ios::binary);MeshCacheIO::readRaw(input,restored);
    check(std::memcmp(&restored.discOpr[0],&negativeZero,sizeof(double))==0 &&
          std::memcmp(&restored.rhs,&nanPayload,sizeof(double))==0,"legacy raw reader lost signed zero or NaN payload");
    const auto first=root/"coefficient_roundtrip",second=root/"coefficient_roundtrip_twice";
    meshFixtures(first,LAMP_MESH_COEFFICIENT_FIXTURE_DIR);
    meshFixtures(second,first,true);
    for(const auto file:{"qmc/mesh_cache.bin","bmc/cache/mesh_cache.bin"})
        check(bytes(first/file)==bytes(second/file),"mesh coefficient roundtrip must be byte stable");
}
DictionaryMesh labelDictionary(int capacity=256) {
    DictionaryMesh d;d.Nx=d.Ny=2;d.maxCellNum=capacity;d.maxInitialIterations=0;
    d.levelSetCount=2;d.minLevel_neg=d.minLevel_pos={0,0};d.maxLevel_local={2,2};return d;
}
void adaptLabels(QuadtreeMesh& mesh,bool coarsenQuadrant) {
    CellField<double> gradient(mesh,0.),target(mesh,2.);
    for(int i=0;i<mesh.cellIndexNum();++i) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
        if(coarsenQuadrant && p.x<.5 && p.y<.5)target.value[i]=1.;
    }
    std::vector<VertField<double>> levelSets;
    std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;
    std::vector<VEdgeField<double>*> vertical;
    std::vector<CellField<double>*> cells{&target};
    mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&target);
}
struct LabelTail { std::size_t begin,end; };
std::array<LabelTail,4> labelTails(const QuadtreeMesh& mesh,const std::string& data,bool zero) {
    check(data.substr(0,4)=="QMC2","label fixture format");
    const auto integer=[&](std::size_t at) {
        check(at+sizeof(int)<=data.size(),"label fixture truncated");
        int value;std::memcpy(&value,data.data()+at,sizeof(value));return value;
    };
    check(integer(4)==2,"QMC2 version changed");
    std::size_t header=28;
    for(int kind=0;kind<4;++kind)header+=sizeof(int)*integer(12+kind*sizeof(int));
    std::size_t offset=header+9*sizeof(int);
    std::array<LabelTail,4> tails;
    const auto inspect=[&](int kind,int capacity,const auto& entity) {
        check(integer(header+(kind+1)*sizeof(int))==capacity,"label capacity changed");
        int free=0;
        for(int label=0;label<capacity;++label)
            if(entity(label).label()==-1)
                check(integer(offset+sizeof(int)*free++)==label,"free-label prefix changed");
        tails[kind]={offset+sizeof(int)*free,offset+sizeof(int)*capacity};
        if(zero)for(auto at=tails[kind].begin;at<tails[kind].end;++at)
            check(data[at]==0,"unused free-label suffix is not zero");
        offset=tails[kind].end;
    };
    inspect(0,mesh.vertLabelCapacity(),[&](int i)->const auto&{return mesh.verts(i);});
    inspect(1,mesh.hEdgeLabelCapacity(),[&](int i)->const auto&{return mesh.hEdges(i);});
    inspect(2,mesh.vEdgeLabelCapacity(),[&](int i)->const auto&{return mesh.vEdges(i);});
    inspect(3,mesh.cellLabelCapacity(),[&](int i)->const auto&{return mesh.cells(i);});
    return tails;
}
std::string writeLabels(const QuadtreeMesh& mesh,const std::filesystem::path& dir,bool zero) {
    std::filesystem::create_directories(dir);mesh.writeCache(dir.generic_string()+"/");
    auto data=bytes(dir/"mesh_cache.bin");const auto tails=labelTails(mesh,data,zero);
    std::ofstream layout(dir/"label_tails.txt");
    for(auto tail:tails)layout<<tail.begin<<' '<<tail.end<<'\n';
    return data;
}
void requireLabelHoles(const QuadtreeMesh& mesh) {
    int parents=0,holes=0,last=0;
    for(int i=0;i<mesh.cellLabelCapacity();++i)if(mesh.cells(i).label()!=-1) {
        last=i;parents+=mesh.cells(i).son[0]!=-1;
    }
    for(int i=0;i<last;++i)holes+=mesh.cells(i).label()==-1;
    check(parents>0 && holes>0 && mesh.cellIndexNum()==52,"fixture needs live parents and recycled holes");
}
void labelFixtures(const std::filesystem::path& output,const std::filesystem::path& input={}) {
    QuadtreeMesh mesh(labelDictionary());
    if(input.empty()) { adaptLabels(mesh,false);adaptLabels(mesh,true); }
    else mesh.readCache((input/"holes").generic_string()+"/");
    requireLabelHoles(mesh);writeLabels(mesh,output/"holes",false);
    adaptLabels(mesh,false);check(mesh.cellIndexNum()==64,"reuse holes in old/new continuation");
    writeLabels(mesh,output/"refined",false);
    adaptLabels(mesh,true);requireLabelHoles(mesh);writeLabels(mesh,output/"coarsened",false);
}
void checkLabelSerialization(const std::filesystem::path& root) {
    QuadtreeMesh full(labelDictionary(9)); // No free vertex slots: an empty valid prefix.
    writeLabels(full,root/"full_vertex_storage",true);
    QuadtreeMesh mesh(labelDictionary());
    writeLabels(mesh,root/"fresh_labels",true);
    adaptLabels(mesh,false);check(mesh.cellIndexNum()==64,"real two-level refinement");
    writeLabels(mesh,root/"refined_labels",true);
    adaptLabels(mesh,true);requireLabelHoles(mesh);
    const auto original=writeLabels(mesh,root/"hole_labels",true);
    check(original==writeLabels(mesh,root/"hole_labels_twice",true),"repeated cache write changed bytes");
    const auto tails=labelTails(mesh,original,true);
    std::string resumed;
    for(unsigned char pattern:{0xA5,0x5A}) {
        auto poisoned=original;
        for(auto tail:tails)std::fill(poisoned.begin()+tail.begin,poisoned.begin()+tail.end,static_cast<char>(pattern));
        const auto source=root/("legacy_labels_"+std::to_string(pattern));
        std::filesystem::create_directories(source);
        {std::ofstream file(source/"mesh_cache.bin",std::ios::binary);file.write(poisoned.data(),poisoned.size());}
        QuadtreeMesh loaded(labelDictionary());loaded.readCache(source.generic_string()+"/");
        requireLabelHoles(loaded);
        check(writeLabels(loaded,source/"rewritten",true)==original,"legacy tail influenced rewritten state");
        check(bytes(source/"mesh_cache.bin")==poisoned,"reader/writer modified legacy input file");
        auto cloned=loaded.clone();
        writeLabels(*cloned,source/"cloned",true);
        loaded.swapStorage(*cloned);writeLabels(loaded,source/"swapped",true);
        adaptLabels(loaded,false);check(loaded.cellIndexNum()==64,"restored free labels cannot refine");
        adaptLabels(loaded,true);requireLabelHoles(loaded);
        const auto continued=writeLabels(loaded,source/"continued",true);
        if(resumed.empty())resumed=continued;
        else check(resumed==continued,"unused legacy tail affected AMR continuation");
    }
    check(original==writeLabels(mesh,root/"source_unchanged",true),"serialization/clone changed source state");
}
void checkOwnedFields(const DictionaryMesh& input, const std::filesystem::path& fixture,
                      const std::string& encoded)
{
    static_assert(!std::is_copy_constructible_v<PicFieldState>);
    static_assert(!std::is_move_constructible_v<PicFieldState>);
    static_assert(!std::is_copy_constructible_v<ImplicitPicFieldHistory>);
    static_assert(!std::is_move_constructible_v<ImplicitPicFieldHistory>);
    static_assert(std::is_same_v<decltype(std::declval<PicFieldState&>().potential().value[0]), const double&>);
    static_assert(std::is_same_v<decltype(std::declval<ImplicitPicFieldHistory&>().nextElectricField().value[0]), const Vector2D&>);
    QuadtreeMesh mesh(input);
    PicFieldState fields(mesh);
    ImplicitPicFieldHistory history(fields);
    check(history.readCheckpoint(fixture), "owned fields read legacy seven-field fixture");
    const auto roundtrip=fixture.parent_path()/"owned_fields";
    std::filesystem::create_directories(roundtrip);
    history.writeCheckpoint(roundtrip);
    check(bytes(roundtrip/PicCheckpointStore::fieldName)==encoded, "owner preserves IPF1 order and bytes");
    for(int i=0;i<mesh.cellIndexNum();++i) {
        check(fields.potential().value[i]==i+.125 && fields.poissonRhs().value[i]==-2*i-.25
            && fields.chargeDensity().value[i]==3*i+.5 && fields.chargeAmount().value[i]==i*.01,
            "owner maps charge amount, density and RHS separately");
    }
    // Constructing history after the initial solve must not silently seed it.
    ImplicitPicFieldHistory fresh(fields);
    check(fresh.nextElectricField().value[0].x==0. && fresh.halfElectricField().value[0].y==0.,
          "fresh history starts at zero independently of current electric field");
    fresh.seedFromCurrent();
    check(fresh.nextElectricField().value[0].x==.2 && fresh.halfElectricField().value[0].y==-.3,
          "restart/AMR rebuild seeds both layers from the current field");
    const auto* currentStorage=fields.electricField().value.get();
    const auto* nextStorage=history.nextElectricField().value.get();
    history.formHalfStep();
    for(int i=0;i<mesh.cellIndexNum();++i)
        check(history.halfElectricField().value[i].x==((i+.2)+(2*i+.4))*.5
            && history.halfElectricField().value[i].y==((-i-.3)+(-3*i-.5))*.5,
            "half-step uses both distinct time layers");
    const CellField<Vector2D> halfBefore(history.halfElectricField());
    history.acceptNext();
    check(fields.electricField().value.get()==currentStorage && history.nextElectricField().value.get()==nextStorage,
          "acceptance preserves current and next storage identity");
    for(int i=0;i<mesh.cellIndexNum();++i)
        check(fields.electricField().value[i].x==2*i+.4 && fields.electricField().value[i].y==-3*i-.5
            && history.halfElectricField().value[i].x==halfBefore.value[i].x
            && history.halfElectricField().value[i].y==halfBefore.value[i].y,
            "acceptance keeps the half-step used by the particle push");
    check(!history.readCheckpoint(roundtrip/"missing"), "missing history reports rebuild requirement");
    history.writeCheckpoint(roundtrip);
    const auto accepted=bytes(roundtrip/PicCheckpointStore::fieldName);
    {std::ofstream bad(roundtrip/PicCheckpointStore::fieldName,std::ios::binary);bad<<"bad header";}
    check(!history.readCheckpoint(roundtrip), "incompatible header is rejected through owner");
    history.writeCheckpoint(roundtrip);
    check(bytes(roundtrip/PicCheckpointStore::fieldName)==accepted, "incompatible header preserves live fields");

    const int oldCells=mesh.cellIndexNum();
    DictionaryMesh::RefineBox box;box.min={0.,0.};box.max={.49,.49};box.level=1;
    mesh.applyRefineBoxes({box});
    check(mesh.cellIndexNum()>oldCells, "owned field contract changes real mesh topology");
    rejects([&]{history.formHalfStep();}, "stale histories must not index the larger mesh");
    rejects([&]{fields.resetAfterMeshChange(&fields.potential());}, "stale remap cell count rejected");
    check(fields.potential().size()==oldCells && fields.potential().value[0]==.125,
          "bad remap rejected before changing live fields");
    CellField<double> remapped(mesh,7.25);
    const int label=remapped.label(0);remapped.setLabel(0,-1);
    rejects([&]{fields.resetAfterMeshChange(&remapped);}, "wrong remap cell order rejected");
    remapped.setLabel(0,label);
    fields.resetAfterMeshChange(&remapped);
    rejects([&]{history.acceptNext();}, "current reset alone does not rebuild implicit histories");
    history.resetAfterMeshChange();
    for(int i=0;i<mesh.cellIndexNum();++i)
        check(fields.potential().value[i]==7.25 && fields.potential().label(i)==mesh.cellIndexToLabel[i]
            && fields.chargeDensity().value[i]==0. && fields.chargeAmount().value[i]==0.
            && fields.poissonRhs().value[i]==0. && fields.electricField().value[i].x==0.
            && history.nextElectricField().value[i].x==0. && history.halfElectricField().value[i].y==0.,
            "fast remap retains potential while rebuilding other fields on new labels");
    remapped.value[0]=99.;
    check(fields.potential().value[0]==7.25, "remapped field storage is independently owned");
    fields.resetAfterMeshChange();
    check(fields.potential().value[0]==0., "ordinary AMR reset clears potential");
}
}
namespace {
void checkRunRequirements(const std::filesystem::path& root) {
    namespace fs=std::filesystem;
    // Cases deliberately distinguish presence from enablement, named explicit
    // geometry from implicit geometry, and all three refinement rule kinds.
    struct GeometryCase { bool have,embedded,named,present,physics,enabled; int rule; bool explicitNeed,implicitNeed; };
    const GeometryCase cases[]={{false,true,true,true,false,true,1,false,false},
        {true,false,true,false,false,false,0,false,false},{true,true,false,false,false,false,0,false,true},
        {true,true,true,false,false,false,0,true,true},{true,false,true,false,false,true,1,false,false},
        {true,false,true,false,true,true,1,true,true},{true,false,true,true,false,true,2,true,true},
        {true,false,true,true,false,true,3,true,true},{true,false,true,true,false,false,1,false,false}};
    int serial=0;
    for(const auto& c:cases)for(auto mode:{PicIntegrator::Explicit,PicIntegrator::Implicit})
    for(int collision=0;collision<5;++collision)for(bool saveRng:{false,true}) {
        PicCaseInput input;input.haveGeometry=c.have;input.geometry.levelSetName=c.named?std::vector<std::string>{"wall"}:std::vector<std::string>{};
        if(c.embedded){DictionaryGeometry::BoundaryEntity boundary;boundary.id=1;boundary.kind=GeometricBoundary::fixed;input.geometry.boundaries={boundary};}
        auto& refinement=input.mesh.geometryRefinement;refinement.present=c.present;refinement.enabled=c.enabled;
        input.mesh.physicsRefinement.present=c.physics;
        if(c.rule==1)refinement.boundaryRules.resize(1);
        if(c.rule==2)refinement.phaseRules.resize(1);
        if(c.rule==3)refinement.regionRules.resize(1);
        input.pic.species.resize(2);input.pic.species[0].name="electron";input.pic.species[1].name="ion";
        input.pic.checkpoint.path=(root/std::to_string(serial++)).string();input.pic.checkpoint.saveRng=saveRng;
        input.movingBoundary.enabled=collision!=0;input.movingBoundary.restart=collision!=1;
        input.movingBoundary.restartPath="custom_wall.bin";
        CollisionsConfig config;config.enabled=collision!=1;
        config.rng_mode=collision==2?CollisionsRngMode::kGlobal:collision==3?CollisionsRngMode::kPerParticle:CollisionsRngMode::kCounterBased;
        const auto store=PicCheckpointStore::forRun(mode,input,collision==0?nullptr:&config);
        check(!fs::exists(input.pic.checkpoint.path),"prepared checkpoint policy must not access disk");
        const bool geometry=mode==PicIntegrator::Explicit?c.explicitNeed:c.implicitNeed;
        const std::vector<std::pair<std::string,bool>> required{{"mesh_cache.bin",true},{"particle_electron_cache.bin",true},
            {"particle_ion_cache.bin",true},{"iGeometry_cache.bin",geometry},{PicCheckpointStore::fieldName,mode==PicIntegrator::Implicit},
            {"custom_wall.bin",geometry && collision>=2},{PicCheckpointStore::collisionName,saveRng && collision==4}};
        // Mutations after construction cannot change component names or paths.
        input.pic.species.clear();input.haveGeometry=false;input.pic.checkpoint.path="unused";input.movingBoundary.restartPath="changed";config.enabled=false;
        fs::path dir;check(store.createDirectory(2,dir),"run policy directory");
        check(PicCheckpointStore::writeMetadata(dir,2,.5,.25),"run policy metadata");
        for(const auto& item:required)touch(dir/item.first);
        check(store.latest()==dir,"complete prepared run checkpoint");
        for(const auto& item:required){fs::remove(dir/item.first);check(store.usable(dir)!=item.second,"prepared component requirement changed");touch(dir/item.first);}
    }
}
void checkExplicitPreflight(const std::filesystem::path& root) {
    const PicCheckpointStore store(root,{});std::filesystem::path dir;
    const auto resumeStep = [&](const CollisionsConfig* config, bool restoreRng) {
        const auto prepared=PicCheckpointReader::prepareExplicitRestart(store,4,config,restoreRng);
        return prepared?std::optional<int>(prepared->outputResumeStep()):std::nullopt;
    };
    check(!PicCheckpointReader::prepareExplicitRestart(store,-1,nullptr,true),"missing preflight must fail");
    check(store.createDirectory(4,dir),"preflight directory");touch(dir/"mesh_cache.bin");
    {std::ofstream file(dir/PicCheckpointStore::metaName);file<<"version future\nstep 7\ndt_code 0.5\n";}
    auto prepared=PicCheckpointReader::prepareExplicitRestart(store,4,nullptr,true);
    check(prepared && prepared->outputResumeStep()==7,"preflight keeps permissive output metadata reader");
    double dt=.25;auto reader=prepared->read(dt);
    check(reader && reader->nextStep()==0 && dt==.25,"preflight must not replace full reader semantics");
    CollisionsConfig config;config.enabled=true;config.rng_mode=CollisionsRngMode::kCounterBased;config.rng_seed=42;
    rejects([&]{resumeStep(&config,true);},"counter state validated before output resume");
    check(resumeStep(&config,false)==7,"disabled RNG restore bypasses counter validation");
    {std::ofstream file(dir/PicCheckpointStore::collisionName);file<<"LAMP_MCC_RNG 1\n2 42 9\n";}
    check(resumeStep(&config,true)==7,"matching counter preflight");
    // Previously saved checkpoints retain their wire header; new writes use LAMP.
    {std::ofstream file(dir/PicCheckpointStore::collisionName);file<<"QUANTUM_MCC_RNG 1\n2 42 9\n";}
    check(resumeStep(&config,true)==7,"legacy counter header remains readable");
    BackgroundGas gas("Ar",1000.,300.);gas.SetParticleMassKg(6.6335209e-26);
    auto xs=std::make_shared<CrossSectionDB>();std::string error;
    check(BuiltInCrossSections::RegisterArV1(xs.get(),"",false,&error),"counter fixture cross sections");
    auto reactions=std::make_shared<ReactionDatabase>();ReactionSet set;
    set.name="electron-Ar";set.projectile="electron";set.target="Ar";
    ReactionChannel channel;channel.name="elastic";channel.cross_section_key="e-Ar:ELASTIC";channel.scatter_model="isotropic";
    set.channels.push_back(channel);check(reactions->Register(set,false,&error),"counter fixture reaction set");
    auto collisionConfig=std::make_shared<CollisionsConfig>(config);
    collisionConfig->parallel=false;collisionConfig->mappings.push_back({"electron","electron-Ar"});
    PlasmaCollisionInputs collisionInputs{std::make_shared<GasState>(gas.Evaluate()),xs,reactions,collisionConfig};
    SpeciesConfig electron;electron.name="electron";electron.m=9.1093837015e-31;electron.q=-1.602176634e-19;electron.weight=1.;
    auto collisions=PicCollisions::create(collisionInputs,{electron},1.,error,&*reader,true);
    check(collisions && collisions->rngState().seed==42 && collisions->rngState().applyCount==9,
        "legacy header restores actual collision RNG state");
    std::mt19937_64 rng(42);PicCheckpointWriter(store,rng,true,0,"ContractPIC").writeCollisions(dir,*collisions);
    check(bytes(dir/PicCheckpointStore::collisionName).rfind("LAMP_MCC_RNG 1",0)==0,
        "new collision checkpoint writer emits LAMP header");
    check(resumeStep(&config,true)==7,"new writer preserves readable counter state");
    config.rng_seed=43;rejects([&]{resumeStep(&config,true);},"counter seed mismatch preflight");
    check(PicCheckpointStore::writeMetadata(dir,-2,0.,.125),"preflight negative metadata");
    check(resumeStep(nullptr,false)==0,"preflight clamps negative output step");
    // A new latest checkpoint between output preparation and startup must not
    // redirect the physical state to another source or another time step.
    check(PicCheckpointStore::writeMetadata(dir,4,.5,.125),"selected metadata");
    prepared=PicCheckpointReader::prepareExplicitRestart(store,-1,nullptr,false);
    check(prepared && prepared->directory()==dir && prepared->outputResumeStep()==4,"latest directory not retained");
    std::filesystem::path newer;check(store.createDirectory(8,newer),"new latest directory");
    touch(newer/"mesh_cache.bin");check(PicCheckpointStore::writeMetadata(newer,8,4.,.5),"new latest metadata");
    check(store.latest()==newer,"fixture did not change latest checkpoint");
    dt=.25;reader=prepared->read(dt);
    check(reader && reader->directory()==dir && reader->nextStep()==4 && dt==.125,
        "output preparation and physical restore selected different checkpoints");
    store.beginWrite(dir);
    check(!prepared->read(dt) && dt==.125,"newly incomplete selected checkpoint fell back to latest");
    check(!resumeStep(nullptr,false),"incomplete explicit preflight");
    std::filesystem::remove(dir/PicCheckpointStore::incompleteName);
    std::filesystem::remove(dir/PicCheckpointStore::metaName);
    check(!prepared->read(dt) && dt==.125,"missing selected metadata fell back to latest or changed dt");
}
void checkLiveCheckpoint(const std::filesystem::path& root) {
    namespace fs=std::filesystem;
    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=32;input.maxInitialIterations=0;input.levelSetCount=1;
    DictionaryGeometry dict;dict.levelSetName={"solid"};dict.innerPoints={{.5,.5}};dict.moving={true};
    dict.linesOnInterface={{{0,false},{1,false},{2,false},{3,false}}};
    dict.phaseName={"solid","gas"};dict.levelSetOperation={{LevelSetOperation::non},{LevelSetOperation::non}};
    dict.levelSetIndexA={{{0,false}},{{0,true}}};dict.levelSetIndexB=dict.levelSetIndexA;
    dict.boundaryLines={{{0,false},{1,false},{2,false},{3,false}},{{0,true},{1,true},{2,true},{3,true}}};
    const std::array<Vector2D,4> corners{{{.3,.3},{.7,.3},{.7,.7},{.3,.7}}};
    for(int i=0;i<4;++i){Line line;line.points[0]=corners[i];line.points[1]=corners[(i+1)%4];line.levelSetIndex=0;dict.lines.push_back(line);}
    auto units=std::make_shared<const UnitsContext>(UnitsContext::MakeSI(UnitsScales{}));
    const auto norm=NormalizationContext::MakeFromUnits(units,NormalizationConfig{});
    for(bool implicit:{false,true}) {
        const auto mode=implicit?PicDomain::Mode::Implicit:PicDomain::Mode::Explicit;
        PicDomain domain(input,input,mode);bool embedded=true,refine=false;
        check(domain.initializeGeometry(&dict,embedded,refine,{}),"live checkpoint geometry");
        PicMovingGeometry moving(mode);MovingEBConfig motion;motion.enabled=true;
        const auto base=root/(implicit?"implicit":"explicit");fs::create_directories(base);
        motion.metricsPath=(base/"metrics.csv").string();motion.erosionMetricsPath=(base/"regions.csv").string();motion.impactStatsPath=(base/"impacts.csv").string();
        check(moving.initialize(domain.mesh(),domain.geometry(),&dict,motion,units.get(),&norm,std::nullopt) && moving.controller(),"live checkpoint wall");
        PicFieldState fields(domain.mesh());ImplicitPicFieldHistory history(fields);
        std::vector<ParticleType> species;species.emplace_back("electron",1.,-1.);species.emplace_back("ion",2.,1.);
        for(auto& sp:species){sp.particle.resize(1);sp.particle[0].position={.1,.2};sp.particle[0].w=3.;sp.number=1;}
        std::mt19937_64 rng(123);const auto original=rng;
        PicCheckpointStore::Requirements requirements{{"electron","ion"},true,implicit,true,{}};
        const PicCheckpointStore actual(base/"actual",requirements),reference(base/"reference",requirements);
        const PicCheckpointWriter writer(actual,rng,true,0,"ContractPIC");
        const auto direct=[&](const fs::path& dir,const std::string& prefix){domain.mesh().writeCache(prefix);domain.writeGeometryCache(prefix);
            moving.writeCache(dir);for(const auto& sp:species)sp.writeCache(prefix);if(implicit)history.writeCheckpoint(dir);};
        for(int step:{1,2}) {
            if(step==2){rng.discard(3);species[0].particle[0].velocity.x=9.;}
            writer.write(step,step*.25,.25,{domain,moving,species,nullptr,implicit?&history:nullptr});
            PicCheckpointWriter(reference,rng,true,0,"ContractPIC").write(step,step*.25,.25,direct);
            std::size_t count=0;
            for(const auto& file:fs::directory_iterator(reference.directory(step))){++count;check(bytes(file.path())==bytes(actual.directory(step)/file.path().filename()),"live checkpoint differs from original component sequence");}
            check(count==static_cast<std::size_t>(std::distance(fs::directory_iterator(actual.directory(step)),fs::directory_iterator{})),"live checkpoint component set changed");
            check(actual.usable(actual.directory(step)),"live checkpoint completion");
        }
        auto expected=original;expected.discard(3);check(rng==expected,"checkpoint serialization consumed RNG");
        std::vector<std::string> order{"mesh_cache.bin","iGeometry_cache.bin","moving_eb_cache.bin","particle_electron_cache.bin","particle_ion_cache.bin"};
        if(implicit)order.push_back(PicCheckpointStore::fieldName);
        for(std::size_t blocked=0;blocked<order.size();++blocked) {
            const int step=10+static_cast<int>(blocked);fs::create_directories(actual.directory(step)/order[blocked]);
            rejects([&]{writer.write(step,step*.25,.25,{domain,moving,species,nullptr,implicit?&history:nullptr});},"live component failure must propagate");
            check(fs::exists(actual.directory(step)/PicCheckpointStore::incompleteName) && actual.select(step).empty(),"failed live write published checkpoint");
            for(std::size_t i=0;i<order.size();++i)if(i!=blocked)
                check(fs::is_regular_file(actual.directory(step)/order[i])==(i<blocked),"live component failure order changed");
            check(actual.latest()==actual.directory(2),"failed write displaced previous checkpoint");
        }
    }
}
}
int main(int argc,char** argv) try
{
    if(argc==3 && std::string(argv[1])=="--write-mesh-fixtures") { meshFixtures(argv[2]);return 0; }
    if(argc==4 && std::string(argv[1])=="--roundtrip-mesh-fixtures") { meshFixtures(argv[3],argv[2]);return 0; }
    if(argc==3 && std::string(argv[1])=="--write-label-fixtures") { labelFixtures(argv[2]);return 0; }
    if(argc==4 && std::string(argv[1])=="--roundtrip-label-fixtures") { labelFixtures(argv[3],argv[2]);return 0; }
    const auto root=std::filesystem::current_path()/
        ("checkpoint_contract_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    if(argc==2 && std::string(argv[1])=="--rng-format-only") {
        checkExplicitPreflight(root/"explicit_preflight");
        std::cout<<"PIC collision RNG legacy-read, LAMP-write and preflight checks passed\n";
        return 0;
    }
    checkCoefficientSerialization(root);
    checkLabelSerialization(root);
    checkRunRequirements(root/"run_requirements");
    checkExplicitPreflight(root/"explicit_preflight");
    checkLiveCheckpoint(root/"live_components");
    PicCheckpointStore::Requirements requirements{{"electron","ion"},true,true,false,{}};
    const PicCheckpointStore store(root,requirements);
    requirements.species.push_back("caller_mutation");
    for(int step:{1,2,9}) {
        std::filesystem::path dir;check(store.createDirectory(step,dir),"create checkpoint directory");
        check(PicCheckpointStore::writeMetadata(dir,step,1.25,.125),"write metadata");
        for(const auto* name:{"mesh_cache.bin","iGeometry_cache.bin","particle_electron_cache.bin",PicCheckpointStore::fieldName})touch(dir/name);
        if(step<9)touch(dir/"particle_ion_cache.bin");
    }
    check(store.latest()==store.directory(2),"latest skips incomplete species set");
    check(store.select(9)==store.directory(9),"explicit step preserves direct selection");
    check(bytes(store.directory(2)/PicCheckpointStore::metaName)=="version 1\nnextStep 2\nt_code 1.25\ndt_code 0.125\n" ||
          bytes(store.directory(2)/PicCheckpointStore::metaName)=="version 1\r\nnextStep 2\r\nt_code 1.25\r\ndt_code 0.125\r\n","legacy metadata layout");
    PicCheckpointStore::Metadata meta;
    check(PicCheckpointStore::readMetadata(store.directory(2),meta) && meta.version==1 && meta.nextStep==2 && meta.t_code==1.25 && meta.dt_code==.125,"metadata values");
    { std::ofstream out(store.directory(1)/PicCheckpointStore::metaName);out<<"version future\nstep 7\nunknown ignored tail\ndt_code 0.5\n"; }
    meta={};check(PicCheckpointStore::readOutputMetadata(store.directory(1),meta) && meta.nextStep==7 && meta.dt_code==.5,"explicit output preflight retains version skip and old step alias");
    meta={};check(PicCheckpointStore::readMetadata(store.directory(1),meta) && meta.nextStep==0,"legacy full reader keeps permissive stream behavior");
    std::mt19937_64 original(12345);original.discard(17);auto expected=original;
    check(PicCheckpointStore::writeRng(store.directory(2),original),"write RNG");
    std::mt19937_64 restored(99);check(PicCheckpointStore::readRng(store.directory(2),restored),"read RNG");
    for(int i=0;i<100;++i)check(restored()==expected(),"RNG continuation");
    check(!PicCheckpointStore::readRng(store.directory(9),restored),"missing RNG reports failure");

    auto explicitRequirements=requirements;explicitRequirements.species={"electron","ion"};
    explicitRequirements.implicitFields=false;explicitRequirements.movingBoundary=true;explicitRequirements.movingBoundaryCache="moving_eb_cache.bin";
    PicCheckpointStore moving(root,explicitRequirements);
    check(moving.latest().empty(),"explicit moving boundary cache required");
    touch(store.directory(2)/"moving_eb_cache.bin");check(moving.latest()==store.directory(2),"explicit moving cache selected");
    explicitRequirements.movingBoundaryCache.clear();
    const PicCheckpointStore relative(std::filesystem::relative(root),explicitRequirements);
    check(relative.latest()==relative.directory(2),"default moving cache resolves once under relative base");
    explicitRequirements.movingBoundaryCache=store.directory(2)/"moving_eb_cache.bin";
    check(PicCheckpointStore(root,explicitRequirements).usable(store.directory(2)),"absolute moving cache path");
    {std::ofstream empty(store.directory(2)/"particle_ion_cache.bin");}
    check(!store.usable(store.directory(2)),"empty component rejected");
    touch(store.directory(2)/"particle_ion_cache.bin");

    PicCheckpointStore::Requirements collisionRequirements{{"electron","ion"},true,true,false,{}};
    collisionRequirements.collisionState=true;
    const PicCheckpointStore collisions(root,collisionRequirements);
    check(collisions.latest().empty(),"counterBased latest selection rejects missing collision state");
    touch(store.directory(2)/PicCheckpointStore::collisionName);
    check(collisions.latest()==store.directory(2),"counterBased checkpoint requires its operator state component");
    {std::ofstream empty(store.directory(2)/PicCheckpointStore::collisionName);}
    check(!collisions.usable(store.directory(2)),"empty collision state rejected");
    store.beginWrite(store.directory(2));
    rejects([&]{collisions.completeWrite(store.directory(2),true);},"missing collision state cannot publish completed checkpoint");
    check(std::filesystem::exists(store.directory(2)/PicCheckpointStore::incompleteName),"collision component failure keeps incomplete marker");
    touch(store.directory(2)/PicCheckpointStore::collisionName);
    collisions.completeWrite(store.directory(2),true);

    const auto checkedPath=root/"checked_output.bin";
    auto checked=FileIO::openOutputFile(checkedPath,std::ios::binary|std::ios::out);
    checked << "payload"; FileIO::finishOutputFile(checked,checkedPath);
    check(bytes(checkedPath)=="payload","checked close preserves bytes");
    checked=FileIO::openOutputFile(checkedPath,std::ios::binary|std::ios::out);
    checked << "buffered"; checked.setstate(std::ios::badbit);
    rejects([&]{FileIO::finishOutputFile(checked,checkedPath);},"post-open write failure must throw");
    checked=FileIO::openOutputFile(checkedPath,std::ios::binary|std::ios::out);
    checked.setstate(std::ios::badbit);
    check(!FileIO::closeOutputFile(checked),"boolean close preserves write failure");

    DictionaryMesh input;input.Nx=input.Ny=2;input.maxCellNum=16;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    CellField<double> phi(mesh,0.),rhs(mesh,0.),rho(mesh,0.),q(mesh,0.);
    CellField<Vector2D> e(mesh,{}),en(mesh,{}),eh(mesh,{});
    for(int i=0;i<mesh.cellIndexNum();++i) {
        phi.value[i]=i+.125;rhs.value[i]=-2*i-.25;rho.value[i]=3*i+.5;q.value[i]=i*.01;
        e.value[i]={i+.2,-i-.3};en.value[i]={2*i+.4,-3*i-.5};eh.value[i]={-i-.6,4*i+.7};
    }
    const auto dir=store.directory(2);
    ImplicitPicFieldCheckpoint::write(mesh,dir,phi,rhs,rho,q,e,en,eh);
    const auto encoded=bytes(dir/PicCheckpointStore::fieldName);
    check(encoded.substr(0,4)=="IPF1" && encoded.size()==12+10*sizeof(double)*mesh.cellIndexNum(),"legacy IPF1 header and seven-field layout");
    checkOwnedFields(input,dir,encoded);
    CellField<double> wrongSize;
    rejects([&]{ImplicitPicFieldCheckpoint::write(mesh,dir,wrongSize,rhs,rho,q,e,en,eh);},"invalid field size must throw");
    check(bytes(dir/PicCheckpointStore::fieldName)==encoded,"field validation precedes truncation");
    const auto blockedComponents=root/"blocked_components";
    std::filesystem::create_directories(blockedComponents);
    for(const auto* name:{"mesh_cache.bin","iGeometry_cache.bin","particle_electron_cache.bin","moving_eb_cache.bin",PicCheckpointStore::fieldName})
        std::filesystem::create_directory(blockedComponents/name);
    const auto blockedPrefix=PicCheckpointStore::pathPrefix(blockedComponents);
    rejects([&]{mesh.writeCache(blockedPrefix);},"mesh open failure must throw");
    ParticleType blockedSpecies("electron",1.,-1.);
    rejects([&]{blockedSpecies.writeCache(blockedPrefix);},"particle open failure must throw");
    DictionaryGeometry geometryConfig;
    ImmersedGeometry<QuadtreeMesh> blockedGeometry(geometryConfig,mesh);
    rejects([&]{blockedGeometry.wirteCache(blockedPrefix);},"geometry open failure must throw");
    MovingEBConfig movingConfig; movingConfig.enabled=true;
    MovingEBController blockedMoving(geometryConfig,movingConfig);
    rejects([&]{blockedMoving.writeCache(blockedComponents);},"moving boundary open failure must throw");
    rejects([&]{ImplicitPicFieldCheckpoint::write(mesh,blockedComponents,phi,rhs,rho,q,e,en,eh);},"implicit field open failure must throw");
    CellField<double> p2(mesh,99.),r2(mesh,99.),n2(mesh,99.),q2(mesh,99.);
    CellField<Vector2D> e2(mesh,{}),en2(mesh,{}),eh2(mesh,{});
    check(ImplicitPicFieldCheckpoint::read(mesh,dir,p2,r2,n2,q2,e2,en2,eh2),"read implicit field history");
    for(int i=0;i<mesh.cellIndexNum();++i) {
        check(p2.value[i]==phi.value[i] && r2.value[i]==rhs.value[i] && n2.value[i]==rho.value[i] && q2.value[i]==q.value[i],"scalar state restored");
        check(e2.value[i].x==e.value[i].x && e2.value[i].y==e.value[i].y && en2.value[i].x==en.value[i].x && en2.value[i].y==en.value[i].y && eh2.value[i].x==eh.value[i].x && eh2.value[i].y==eh.value[i].y,"three vector time layers restored");
    }
    {std::ofstream out(dir/PicCheckpointStore::fieldName,std::ios::binary);out.write(encoded.data(),encoded.size()-1);}
    check(!ImplicitPicFieldCheckpoint::read(mesh,dir,p2,r2,n2,q2,e2,en2,eh2),"truncated field payload rejected");
    std::filesystem::create_directories(root/"unrelated");
    store.keepLatest(0);check(std::filesystem::exists(store.directory(1)),"keep all policy");
    store.keepLatest(2);check(std::filesystem::exists(store.directory(1)),"incomplete checkpoint must not consume retention slot");
    store.keepLatest(1);check(!std::filesystem::exists(store.directory(1)) && std::filesystem::exists(store.directory(2)) && std::filesystem::exists(store.directory(9)) && std::filesystem::exists(root/"unrelated"),"retention counts only usable checkpoints");

    const PicCheckpointStore written(root/"written", {});
    const PicCheckpointWriter writer(written, original, true, 2, "ContractPIC");
    for (int step : {1, 2, 3})
    {
        writer.write(step, step * .125, .125, [&](const std::filesystem::path& dir, const std::string& prefix) {
            check(std::filesystem::exists(dir/PicCheckpointStore::incompleteName)
                && !written.usable(dir) && written.select(step).empty(),"pending write cannot be selected");
            PicCheckpointStore::Metadata metadata;
            check(PicCheckpointStore::readMetadata(dir, metadata) && metadata.nextStep == step
                && metadata.t_code == step * .125 && metadata.dt_code == .125,
                "metadata must precede physical components");
            std::mt19937_64 loaded;
            check(PicCheckpointStore::readRng(dir, loaded) && loaded == original,
                "current RNG must precede physical components without advancing it");
            check(prefix == PicCheckpointStore::pathPrefix(dir), "component cache prefix changed");
            touch(dir/"mesh_cache.bin");
            check(!written.usable(dir),"all components present still requires completion");
        });
        check(written.usable(written.directory(step))
            && !std::filesystem::exists(written.directory(step)/PicCheckpointStore::incompleteName),"completed write removes marker");
    }
    check(!std::filesystem::exists(written.directory(1)) && std::filesystem::exists(written.directory(2))
        && std::filesystem::exists(written.directory(3)), "retention after component writes");
    bool componentFailed = false;
    try { writer.write(4, .5, .125, [](const auto&, const auto&) { throw std::runtime_error("component"); }); }
    catch (const std::runtime_error&) { componentFailed = true; }
    check(componentFailed && std::filesystem::exists(written.directory(2)),
        "component exception must propagate before retention");
    check(std::filesystem::exists(written.directory(4)/PicCheckpointStore::metaName)
        && !std::filesystem::exists(written.directory(4)/"mesh_cache.bin"),
        "partial files remain available for diagnosis");
    check(written.latest()==written.directory(3) && written.select(4).empty()
        && std::filesystem::exists(written.directory(4)/PicCheckpointStore::incompleteName),"failed checkpoint excluded from auto and explicit selection");
    rejects([&]{writer.write(5,.625,.125,[](const auto&,const auto&){});},"missing component prevents completion");
    rejects([&]{writer.write(6,.75,.125,[&](const auto& path,const auto&){
        touch(path/"mesh_cache.bin");
        std::filesystem::remove(path/PicCheckpointStore::incompleteName);
        std::filesystem::create_directory(path/PicCheckpointStore::incompleteName);
        touch(path/PicCheckpointStore::incompleteName/"block_removal");
    });},"marker removal failure must propagate");
    check(written.latest()==written.directory(3) && std::filesystem::exists(written.directory(2)),"completion failures preserve older usable checkpoints");

    const PicCheckpointStore warned(root/"warned", {});
    std::filesystem::path warningDir;
    check(warned.createDirectory(1, warningDir), "warning directory");
    std::filesystem::create_directory(warningDir/PicCheckpointStore::metaName);
    std::filesystem::create_directory(warningDir/PicCheckpointStore::rngName);
    int componentCalls = 0;
    const auto writeBlocked=[&]{PicCheckpointWriter(warned, original, true, 0, "ContractPIC").write(1, .125, .125,
        [&](const auto&, const auto&) { ++componentCalls; });};
    rejects(writeBlocked,"metadata failure must throw");
    check(componentCalls==0,"components ran after metadata failure");
    std::filesystem::remove(warningDir/PicCheckpointStore::metaName);
    rejects(writeBlocked,"RNG failure must throw");
    check(componentCalls==0 && warned.select(1).empty(),"RNG failure accepted checkpoint or ran components");
    touch(root/"blocked");
    const PicCheckpointStore blocked(root/"blocked", {});
    rejects([&]{PicCheckpointWriter(blocked, original, true, 0, "ContractPIC").write(1, .125, .125,
        [&](const auto&, const auto&) { ++componentCalls; });},"directory creation failure must throw");
    check(componentCalls == 0, "components ran after directory creation failed");

    rejects([&]{drivePicStepSession(0,2,.125,[](const PicStepTime&){return 0;},[&](PicStepSession& session){
        rejects([&]{session.advance(true);},"checkpoint exception must escape advance");
        check(session.nextStep()==0,"failed checkpoint advanced accepted step");
        rejects([&]{session.advance();},"failed session allowed another physical step");
        return 0;
    },[&](const PicStepTime&,bool){writeBlocked();return 0;});},"suppressed checkpoint failure was reported as success");

    const PicCheckpointStore noRng(root/"no_rng", {});
    std::filesystem::create_directories(noRng.directory(1));
    touch(noRng.directory(1)/PicCheckpointStore::collisionName);
    PicCheckpointWriter(noRng, original, false, 0, "ContractPIC").write(1, .125, .125,
        [&](const auto& dir, const auto&) {
            check(!std::filesystem::exists(dir/PicCheckpointStore::rngName), "disabled RNG save wrote a file");
            check(!std::filesystem::exists(dir/PicCheckpointStore::collisionName), "reused checkpoint kept stale collision RNG state");
            touch(dir/"mesh_cache.bin");
        });
    const PicCheckpointStore selected(root/"selected", {{"electron", "ion"}, false, false, false, {}});
    std::vector<ParticleType> saved;
    saved.emplace_back("electron", 2., -3.);
    saved.emplace_back("ion", 4., 3.);
    Particle inside, outside;
    inside.position = {.75, .75}; inside.velocity = {.25, -.5}; inside.vz = 1.5;
    inside.w = 7.; inside.rngState = 8642; inside.inCellL = -1; inside.inBackgroundCellL = 99;
    outside = inside; outside.position = {2., 2.};
    saved[0].particle = {inside, outside}; saved[0].number = 2;
    saved[1].particle = {inside}; saved[1].number = 1;
    const auto writeComponents = [&](const auto&, const std::string& prefix) {
        mesh.writeCache(prefix);
        for (const auto& sp : saved) sp.writeCache(prefix);
    };
    PicCheckpointWriter(selected, original, true, 0, "ContractPIC").write(7, .875, .125, writeComponents);
    std::filesystem::path incomplete;
    check(selected.createDirectory(99, incomplete), "incomplete checkpoint directory");
    check(PicCheckpointStore::writeMetadata(incomplete, 99, 0., .25), "incomplete metadata");
    double restoredDt = .5;
    auto reader = PicCheckpointReader::select(selected, -1, restoredDt, "ContractPIC");
    check(reader && reader->directory() == selected.directory(7) && reader->nextStep() == 7
        && reader->metadata().t_code == .875 && restoredDt == .125, "selected metadata and restored dt");
    PicCheckpointWriter(selected, original, true, 0, "ContractPIC").write(8, 1., .125, writeComponents);
    check(selected.latest() == selected.directory(8) && reader->nextStep() == 7,
        "selected checkpoint must not drift to a newly created directory");

    QuadtreeMesh restoredMesh(input);
    DictionaryMesh::RefineBox restoreBox;
    restoreBox.min = {0., .5}; restoreBox.max = {1., 1.}; restoreBox.level = 1;
    restoredMesh.applyRefineBoxes({restoreBox});
    check(restoredMesh.cellIndexNum() != mesh.cellIndexNum(), "distinct pre-restore mesh");
    check(reader->readMesh(restoredMesh, PicCheckpointReader::MeshCache::Required)
        && restoredMesh.cellIndexNum() == mesh.cellIndexNum(), "actual mesh restored");
    std::vector<ParticleType> loaded;
    loaded.emplace_back("electron", 9., 9.); loaded.emplace_back("ion", 9., 9.);
    std::mt19937_64 loadedRng(987);
    check(reader->readParticles(restoredMesh, loaded, loadedRng, true) && loadedRng == original,
        "live RNG restored without advancing it");
    check(loaded[0].number == 2 && loaded[1].number == 1 && loaded[0].mass() == 2.
        && loaded[0].charge() == -3., "particle payload and species restored");
    const auto& rebound = loaded[0].particle[0];
    check(rebound.inCellL == restoredMesh.findCellLabel(inside.position)
        && rebound.inBackgroundCellL == restoredMesh.findBackgroundCellLabel(inside.position)
        && rebound.velocity.x == .25 && rebound.velocity.y == -.5 && rebound.vz == 1.5
        && rebound.w == 7. && rebound.rngState == 8642, "particle rebound without changing its state");
    check(loaded[0].particle[1].inCellL == -1 && loaded[0].particle[1].inBackgroundCellL == -1,
        "outside particle retained for the original physical culling stage");
    loadedRng.seed(543); const auto untouchedRng = loadedRng;
    check(reader->readParticles(restoredMesh, loaded, loadedRng, false) && loadedRng == untouchedRng,
        "disabled RNG restore preserved caller state");

    auto partial = PicCheckpointReader::select(selected, 99, restoredDt, "ContractPIC");
    check(partial && partial->nextStep() == 99, "explicit selection retains incomplete-directory behavior");
    const int cellsBefore = restoredMesh.cellIndexNum();
    check(partial->readMesh(restoredMesh, PicCheckpointReader::MeshCache::Optional)
        && restoredMesh.cellIndexNum() == cellsBefore, "explicit optional mesh behavior");
    check(!partial->readMesh(restoredMesh, PicCheckpointReader::MeshCache::Required),
        "implicit missing mesh rejection");
    saved[0].writeCache(PicCheckpointStore::pathPrefix(incomplete));
    check(PicCheckpointStore::writeRng(incomplete, original), "partial checkpoint RNG");
    loadedRng.seed(543);
    check(!partial->readParticles(restoredMesh, loaded, loadedRng, true) && loadedRng == original
        && loaded[0].particle[0].inCellL == -1, "RNG and earlier species restore precede missing species; no relabel yet");
    saved[1].writeCache(PicCheckpointStore::pathPrefix(incomplete));

    std::filesystem::path noRngDir;
    check(selected.createDirectory(6, noRngDir), "missing RNG directory");
    check(PicCheckpointStore::writeMetadata(noRngDir, -2, 0., 0.), "legacy nonpositive metadata");
    writeComponents(noRngDir, PicCheckpointStore::pathPrefix(noRngDir));
    restoredDt = .75;
    auto missingRng = PicCheckpointReader::select(selected, 6, restoredDt, "ContractPIC");
    check(missingRng && missingRng->nextStep() == 0 && restoredDt == .75,
        "legacy negative step clamp and nonpositive dt preserved");
    loadedRng.seed(543);
    check(missingRng->readParticles(restoredMesh, loaded, loadedRng, true) && loadedRng == untouchedRng,
        "missing RNG warning still allows particle restore");
    restoredDt = .125 + 5e-17; const double withinTolerance = restoredDt;
    check(PicCheckpointReader::select(selected, 7, restoredDt, "ContractPIC").has_value()
        && restoredDt == withinTolerance, "original dt override threshold");
    check(!PicCheckpointReader::select(selected, 1234, restoredDt, "ContractPIC")
        && restoredDt == withinTolerance, "missing selection leaves dt unchanged");
    reader->report();
    std::cout<<"PIC checkpoint selection, restoration, metadata, RNG, IPF1 history and retention checks passed\n";
    return 0;
}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
