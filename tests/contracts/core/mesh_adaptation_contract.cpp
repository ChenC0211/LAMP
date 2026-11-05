// Core mesh contract: region refinement, QMC2 roundtrip, coarsening and repeated
// refinement on the native mesh, with field transfer, 2:1 balance and face closure
// checked against analytic references. No dictionary or App is involved.
#include "Core/Field/CellField.h"
#include "Core/Field/NativeTargetLevelField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/TecplotMeshOutput.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "Core/Mesh/MeshCacheIO.h" // Task40 full-capacity member oracle.
#include <sstream>
#include <cfenv>
#include <cerrno>
#include <thread>
#include <exception>
#include <type_traits>
#include <atomic>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace
{
namespace fs = std::filesystem;

void require(bool value, const std::string& message)
{
    if (!value) throw std::runtime_error(message);
}

struct Setup
{
    const char* name;
    Vector2D min, max;
    int nx, ny, level;
    Vector2D regionMin, regionMax;
    int refinedCells; // exact leaf count after the first refinement; 0 = not asserted
};

DictionaryMesh makeInput(const Setup& s)
{
    DictionaryMesh input;
    input.minCoord = s.min; input.maxCoord = s.max;
    input.Nx = s.nx; input.Ny = s.ny;
    input.maxCellNum = 4096; input.maxInitialIterations = 0;
    input.maxLevel_local.assign(1, s.level);
    auto& refinement = input.geometryRefinement;
    refinement.present = refinement.enabled = true;
    DictionaryMesh::GeometryRefinementConfig::RegionRule rule;
    rule.name = "region"; rule.shape = DictionaryMesh::GeometryRefinementConfig::RegionShape::kBox;
    rule.min = s.regionMin; rule.max = s.regionMax; rule.maxLevel = s.level;
    refinement.regionRules.push_back(rule);
    return input;
}

double reference(const QuadtreeMesh& mesh, Vector2D point)
{
    const double x = (point.x-mesh.minCoord.x)/(mesh.maxCoord.x-mesh.minCoord.x);
    const double y = (point.y-mesh.minCoord.y)/(mesh.maxCoord.y-mesh.minCoord.y);
    return 1. + 2.*x - 3.*y + 4.*x*y;
}

struct Probe
{
    std::vector<VertField<double>> nodal;
    CellField<double> constant;
    explicit Probe(QuadtreeMesh& mesh) : nodal{VertField<double>(mesh,0.)}, constant(mesh,2.)
    {
        for (int i=0;i<mesh.vertIndexNum();++i)
            nodal[0].value[i] = reference(mesh,mesh.verts(mesh.vertIndexToLabel[i]).coord());
    }
};

int adapt(QuadtreeMesh& mesh, Probe& probe, QuadtreeMeshRefiner& refiner, bool coarsen)
{
    int sweeps = 0;
    for (;sweeps <= mesh.maxLevel()+2;++sweeps) {
        CellField<double> target(mesh,0.), gradient(mesh,0.);
        if (!coarsen) require(refiner.buildTargetLevels(mesh,target), "configured region produced no targets");
        if (!mesh.mayAdaptToTargets(target)) return sweeps;
        std::vector<VertField<double>> levelSets;
        std::vector<VertField<double>*> vertices{&probe.nodal[0]};
        std::vector<HEdgeField<double>*> horizontal; std::vector<VEdgeField<double>*> vertical;
        std::vector<CellField<double>*> cells{&target,&probe.constant};
        mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&target);
    }
    throw std::runtime_error("refinement targets did not settle within the refinement-level bound");
}

void measure(const QuadtreeMesh& mesh, const Probe& probe)
{
    require(probe.nodal[0].size() == mesh.vertIndexNum() && probe.constant.size() == mesh.cellIndexNum(),
        "stale transferred probe fields");
    double area = 0., nodalError = 0., cellError = 0., closureError = 0.;
    for (int i=0;i<mesh.vertIndexNum();++i) {
        const auto& vertex = mesh.verts(mesh.vertIndexToLabel[i]);
        require(vertex.index() == i && probe.nodal[0].label(i) == vertex.label(), "invalid vertex/field identity");
        const double error = std::abs(probe.nodal[0].value[i]-reference(mesh,vertex.coord()));
        require(std::isfinite(error), "nonfinite transferred nodal probe");
        nodalError = std::max(nodalError,error);
    }
    FaceSegmentGraph graph; graph.build(mesh);
    for (int i=0;i<mesh.cellIndexNum();++i) {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        require(cell.index() == i && cell.son[0] == -1 && probe.constant.label(i) == cell.label(),
            "invalid leaf/field identity");
        area += cell.area();
        const double error = std::abs(probe.constant.value[i]-2.);
        require(std::isfinite(error), "nonfinite transferred cell probe");
        cellError = std::max(cellError,error);
        Vector2D closure{0.,0.}; const auto range = graph.segmentRangeForCellIndex(i);
        for (auto j=range.first;j<range.second;++j) {
            const auto& segment = graph.segments()[j];
            require(segment.K == cell.label() && segment.Af > 0. && segment.df > 0., "invalid face segment geometry");
            closure = closure + segment.n*segment.Af;
            if (segment.L >= 0)
                require(std::abs(cell.level()-mesh.cells(segment.L).level()) <= 1, "mesh violates 2:1 face-neighbor balance");
        }
        closureError = std::max(closureError,std::hypot(closure.x,closure.y)/(cell.size().x+cell.size().y));
    }
    const double expectedArea = (mesh.maxCoord.x-mesh.minCoord.x)*(mesh.maxCoord.y-mesh.minCoord.y);
    require(std::abs(area/expectedArea-1.) <= 1e-12, "mesh area changed");
    require(nodalError <= 1e-11, "nodal field transfer reference failed");
    require(cellError <= 1e-11, "cell field transfer reference failed");
    require(closureError <= 1e-12, "face closure reference failed");
}

std::vector<char> bytes(const fs::path& path)
{
    std::ifstream file(path,std::ios::binary); require(bool(file), "missing mesh cache");
    return {std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
}

void identityMutations(QuadtreeMesh& mesh)
{
    const NativeMeshIdentity identity(mesh);
    const auto equivalent=mesh.clone();identity.validate(*equivalent);
    const auto rejectedMutation=[&](auto&& alter,auto&& restore,const char* name){
        alter();bool rejected=false;
        try{identity.validate(mesh);}catch(const std::invalid_argument&){rejected=true;}
        restore();require(rejected,std::string("native identity missed ")+name);identity.validate(mesh);
    };
    auto& leaf=mesh.cells(mesh.cellIndexToLabel[mesh.cellIndexNum()-1]);
    const auto center=leaf.centerCoord();
    rejectedMutation([&]{leaf.setCenterCoord({std::nextafter(center.x,INFINITY),center.y});},
        [&]{leaf.setCenterCoord(center);},"one-ULP leaf geometry change");
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[mesh.vertIndexNum()-1]);const auto coordinate=vertex.coord();
    rejectedMutation([&]{vertex.setCoord({coordinate.x,std::nextafter(coordinate.y,INFINITY)});},
        [&]{vertex.setCoord(coordinate);},"late native vertex change");
    auto& edge=mesh.hEdges(mesh.hEdgeIndexToLabel[mesh.hEdgeIndexNum()-1]);const auto endpoint=edge.edgeToVert[1];
    rejectedMutation([&]{edge.edgeToVert[1]=edge.edgeToVert[0];},[&]{edge.edgeToVert[1]=endpoint;},"late native edge connectivity change");
    int parent=-1;
    for(int i=0;i<mesh.cellIndexNum()&&parent<0;++i)parent=mesh.cells(mesh.cellIndexToLabel[i]).parent;
    require(parent>=0,"native identity test requires the existing real AMR ancestor");
    auto& ancestor=mesh.cells(parent);const auto son=ancestor.son[0];
    require(son>=0,"real AMR ancestor has no child");
    rejectedMutation([&]{ancestor.son[0]=-1;},[&]{ancestor.son[0]=son;},"ancestor-only child connectivity change");
    const int cache=leaf.supCellIndex;leaf.supCellIndex=cache==0?1:0;identity.validate(mesh);leaf.supCellIndex=cache;
    identity.validate(mesh);
}

void verify(const Setup& s, const fs::path& scratch)
{
    const std::string context = std::string(s.name) + ": ";
    try {
        const auto input = makeInput(s);
        const int background = input.Nx*input.Ny;
        QuadtreeMesh mesh(input); Probe probe(mesh); QuadtreeMeshRefiner refiner(input);
        require(refiner.attachGeometry(nullptr,nullptr), "region-only refinement did not initialize");
        require(mesh.cellIndexNum() == background, "initial mesh is not the background lattice");
        measure(mesh,probe);

        adapt(mesh,probe,refiner,false);
        require(mesh.cellIndexNum() > background, "refinement did not change the mesh");
        require(s.refinedCells == 0 || mesh.cellIndexNum() == s.refinedCells,
            "unexpected refined leaf count " + std::to_string(mesh.cellIndexNum()));
        measure(mesh,probe);

        identityMutations(mesh);

        const auto refined = scratch/"refined", roundtrip = scratch/"roundtrip";
        fs::create_directories(refined); fs::create_directories(roundtrip);
        mesh.writeCache(refined.generic_string()+"/");
        QuadtreeMesh restored(input); restored.readCache(refined.generic_string()+"/");
        restored.writeCache(roundtrip.generic_string()+"/");
        require(bytes(refined/"mesh_cache.bin") == bytes(roundtrip/"mesh_cache.bin"), "QMC2 roundtrip changed bytes");

        // QMC2 stores the mesh only; probe fields are copied with the same native indexing.
        Probe continued = probe; measure(restored,continued);
        adapt(restored,continued,refiner,true);
        require(restored.cellIndexNum() == background, "coarsening did not return to the background leaves");
        measure(restored,continued);
        adapt(restored,continued,refiner,false);
        measure(restored,continued);

        const auto original = TecplotMeshOutput::nodalZone(mesh,{"bilinear_probe"},probe.nodal);
        const auto repeated = TecplotMeshOutput::nodalZone(restored,{"bilinear_probe"},continued.nodal);
        require(original.connectivity == repeated.connectivity && original.fields[0].values == repeated.fields[0].values
            && original.points.size() == repeated.points.size(), "refine/coarsen/refine changed nodal topology or values");
        for (std::size_t i=0;i<original.points.size();++i)
            require(original.points[i].x == repeated.points[i].x && original.points[i].y == repeated.points[i].y,
                "repeated refinement changed coordinates");
        std::cout << "mesh_adaptation " << s.name << ": " << background << " -> " << mesh.cellIndexNum()
                  << " -> " << background << " -> " << restored.cellIndexNum() << " leaves\n";
    } catch (const std::exception& error) {
        throw std::runtime_error(context + error.what());
    }
}
// A passive observer of the REAL original CellField update. It invokes
// that production implementation unchanged; only intermediate facts are kept.
struct LegacyTargetTrace : CellField<double>
{
    Vector2D peakLower,peakUpper;bool observe=false;
    int dilutedChildren=0;double minimumChildTarget=std::numeric_limits<double>::infinity();
    double maximumChildTarget=-std::numeric_limits<double>::infinity();
    LegacyTargetTrace(const QuadtreeMesh& mesh,double initial,Vector2D lower,Vector2D upper,bool record)
        :CellField<double>(mesh,initial),peakLower(lower),peakUpper(upper),observe(record){}
    void updateWithMesh(const QuadtreeMesh& mesh)
    {
        CellField<double>::updateWithMesh(mesh);
        if(!observe)return;
        for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);const auto x=c.centerCoord();
            if(c.level()==3&&x.x>peakLower.x&&x.x<peakUpper.x&&x.y>peakLower.y&&x.y<peakUpper.y){
                require(std::isfinite(value[i]),"legacy target witness became nonfinite");
                minimumChildTarget=std::min(minimumChildTarget,value[i]);maximumChildTarget=std::max(maximumChildTarget,value[i]);
                dilutedChildren+=value[i]<2.5;
            }
        }
    }
};
void nativeTargetTransfer()
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=4;
    input.maxCellNum=4096;input.maxInitialIterations=0;input.maxLevel_local={3};
    QuadtreeMesh source(input);
    DictionaryMesh::RefineBox box;box.min={0,0};box.max={1,1};box.level=2;source.applyRefineBoxes({box});
    for(int i=0;i<source.cellIndexNum();++i)require(source.cells(source.cellIndexToLabel[i]).level()==2,"target fixture is not real uniform L2");
    const int peakLabel=source.findCellLabel({.40625,.40625},0.);
    const auto& peak=source.cells(peakLabel);const auto center=peak.centerCoord(),size=peak.size();
    const Vector2D lower{center.x-size.x/2,center.y-size.y/2},upper{center.x+size.x/2,center.y+size.y/2};
    require(peak.level()==2&&lower.x>0&&lower.y>0&&upper.x<1&&upper.y<1,"target fixture peak is not an interior real leaf");
    const NativeMeshIdentity sourceIdentity(source);
    // The oracle refers ONLY to the ORIGINAL physical peak rectangle. No new
    // target field, transferred value or provenance helper computes expected.
    const auto frozenDesired=[&](const Cell& c){const auto x=c.centerCoord(),h=c.size();
        const bool overlap=std::max(x.x-h.x/2,lower.x)<std::min(x.x+h.x/2,upper.x)
            &&std::max(x.y-h.y/2,lower.y)<std::min(x.y+h.y/2,upper.y);
        return overlap?3.:1.;
    };
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
    auto legacy=source.clone();
    LegacyTargetTrace legacyTargets(*legacy,1.,lower,upper,true),legacyGradient(*legacy,0.,lower,upper,false);
    for(int i=0;i<legacy->cellIndexNum();++i)legacyTargets.value[i]=frozenDesired(legacy->cells(legacy->cellIndexToLabel[i]));
    require(legacy->mayAdaptToTargets(legacyTargets),"legacy sharp native target did not request adaptation");
    std::vector<LegacyTargetTrace*> legacyFields{&legacyTargets};
    legacy->adapt(levelSets,legacyGradient,vertices,horizontal,vertical,legacyFields,&legacyTargets);
    const int legacyPeakLevel=legacy->cells(legacy->findCellLabel(center,0.)).level();
    require(legacyTargets.dilutedChildren>=4&&legacyTargets.minimumChildTarget<2.5&&legacyTargets.maximumChildTarget<2.5,
        "old real CellField witness did not dilute all actual target3 children toward2");
    require(legacyPeakLevel<3,"old real CellField witness unexpectedly retained the requested peak L3");
    std::cout<<"NATIVE_TARGET_LEGACY diluted_children="<<legacyTargets.dilutedChildren
        <<" child_target_min="<<legacyTargets.minimumChildTarget<<" child_target_max="<<legacyTargets.maximumChildTarget
        <<" final_peak_level="<<legacyPeakLevel<<" leaves="<<legacy->cellIndexNum()<<'\n';
    sourceIdentity.validate(source);

    auto fixed=source.clone();
    LAMP::Mesh::NativeTargetLevelField targets(*fixed,1.),gradient(*fixed,0.);
    for(int i=0;i<fixed->cellIndexNum();++i)targets.value[i]=frozenDesired(fixed->cells(fixed->cellIndexToLabel[i]));
    std::vector<LAMP::Mesh::NativeTargetLevelField*> fields{&targets};
    const int bound=fixed->maxLevel()+2;int sweeps=0;
    while(fixed->mayAdaptToTargets(targets)){
        require(sweeps<bound,"discrete native target did not settle within original bound");
        fixed->adapt(levelSets,gradient,vertices,horizontal,vertical,fields,&targets);++sweeps;
        for(int i=0;i<fixed->cellIndexNum();++i)
            require(targets.label(i)==fixed->cellIndexToLabel[i]&&targets.value[i]==frozenDesired(fixed->cells(fixed->cellIndexToLabel[i])),
                "target-only transfer differs from original-coordinate maximum");
    }
    int peakLeaves=0;std::array<int,4> levels{};
    for(int i=0;i<fixed->cellIndexNum();++i){const auto& c=fixed->cells(fixed->cellIndexToLabel[i]);
        require(c.level()>=1&&c.level()<=3,"fixed native target escaped real1..3 levels");++levels[c.level()];
        if(frozenDesired(c)==3){require(c.level()==3,"fixed target lost actual original-rectangle L3");++peakLeaves;}
    }
    require(sweeps>0&&peakLeaves==4&&levels[1]>0&&levels[2]>0&&levels[3]==4,"fixed native target lacks actual1/2/3 graded topology");
    // Existing analytic geometry/2:1/closure gates apply to the real result.
    Probe geometryProbe(*fixed);measure(*fixed,geometryProbe);
    const NativeMeshIdentity fixedIdentity(*fixed);
    LAMP::Mesh::NativeTargetLevelField repeated(*fixed,0.);
    for(int i=0;i<fixed->cellIndexNum();++i)repeated.value[i]=frozenDesired(fixed->cells(fixed->cellIndexToLabel[i]));
    require(!fixed->mayAdaptToTargets(repeated),"original-coordinate frozen target is not stable after rebind");
    fixedIdentity.validate(*fixed);sourceIdentity.validate(source);
    std::cout<<"NATIVE_TARGET_ONLY_PASS sweeps="<<sweeps<<" levels1/2/3="<<levels[1]<<'/'<<levels[2]<<'/'<<levels[3]
        <<" actual_peak_leaves="<<peakLeaves<<" source_unchanged=1\n";
}
// These observers inspect the real pre-transfer records, then invoke the
// unchanged public edge consumer. They do not manufacture samples or values.
struct ThinHProbe : HEdgeField<double>
{
    int twoPointNormal=0,onePointBoundary=0;
    ThinHProbe(QuadtreeMesh& mesh,double value):HEdgeField<double>(mesh,value){}
    void updateWithMesh(const QuadtreeMesh& mesh)
    {
        for(int i=0;i<mesh.hEdgeIndexNum();++i){const auto& e=mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
            if(e.preIndex()!=-1)continue;const auto& c=e.relatedPreCoord;
            if(c.size()==2&&c[0].x==c[1].x&&c[0].y!=c[1].y&&e.midpointCoord().y==(c[0].y+c[1].y)/2)++twoPointNormal;
            if(c.size()==1&&e.midpointCoord().y==c[0].y)++onePointBoundary;
        }
        HEdgeField<double>::updateWithMesh(mesh);
    }
};
struct ThinVProbe : VEdgeField<double>
{
    int twoPointNormal=0,onePointBoundary=0;
    ThinVProbe(QuadtreeMesh& mesh,double value):VEdgeField<double>(mesh,value){}
    void updateWithMesh(const QuadtreeMesh& mesh)
    {
        for(int i=0;i<mesh.vEdgeIndexNum();++i){const auto& e=mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
            if(e.preIndex()!=-1)continue;const auto& c=e.relatedPreCoord;
            if(c.size()==2&&c[0].y==c[1].y&&c[0].x!=c[1].x&&e.midpointCoord().x==(c[0].x+c[1].x)/2)++twoPointNormal;
            if(c.size()==1&&e.midpointCoord().x==c[0].x)++onePointBoundary;
        }
        VEdgeField<double>::updateWithMesh(mesh);
    }
};
void thinMeshTransfer(int nx,int ny)
{
    DictionaryMesh input;input.Nx=nx;input.Ny=ny;input.minCoord={11./64,5./16};
    input.maxCoord={input.minCoord.x+nx/64.,input.minCoord.y+ny/64.};
    input.maxCellNum=512;input.maxInitialIterations=0;input.maxLevel_local={2};
    QuadtreeMesh mesh(input);const int background=nx*ny;
    const Vector2D seed{input.minCoord.x+(nx==1?.5:1.5)/64,input.minCoord.y+(ny==1?.5:1.5)/64};
    const auto& root=mesh.cells(mesh.findCellLabel(seed,0));const auto center=root.centerCoord(),size=root.size();
    const Vector2D rootLower{center.x-size.x/2,center.y-size.y/2},rootUpper{center.x+size.x/2,center.y+size.y/2};
    const Vector2D peakLower=center,peakUpper=rootUpper;
    ThinHProbe hc(mesh,2.),hn(mesh,0.);ThinVProbe vc(mesh,-3.),vn(mesh,0.);
    const auto hValue=[](Vector2D p){return 1.25+.5*p.y;};const auto vValue=[](Vector2D p){return -.75+.25*p.x;};
    for(int i=0;i<mesh.hEdgeIndexNum();++i)hn.value[i]=hValue(mesh.hEdges(mesh.hEdgeIndexToLabel[i]).midpointCoord());
    for(int i=0;i<mesh.vEdgeIndexNum();++i)vn.value[i]=vValue(mesh.vEdges(mesh.vEdgeIndexToLabel[i]).midpointCoord());
    double maximumError=0;int hanging=0;
    const auto check=[&]{
        require(hc.size()==mesh.hEdgeIndexNum()&&hn.size()==mesh.hEdgeIndexNum()&&vc.size()==mesh.vEdgeIndexNum()&&vn.size()==mesh.vEdgeIndexNum(),"thin edge field dimensions differ");
        for(int i=0;i<mesh.hEdgeIndexNum();++i){const auto& e=mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
            require(hc.label(i)==e.label()&&hn.label(i)==e.label(),"thin H field labels differ");
            const auto error=std::max(std::abs(hc.value[i]-2.),std::abs(hn.value[i]-hValue(e.midpointCoord())));
            require(std::isfinite(error)&&error<=1e-11,"thin H constant/identifiable normal-linear transfer failed");maximumError=std::max(maximumError,error);}
        for(int i=0;i<mesh.vEdgeIndexNum();++i){const auto& e=mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
            require(vc.label(i)==e.label()&&vn.label(i)==e.label(),"thin V field labels differ");
            const auto error=std::max(std::abs(vc.value[i]+3.),std::abs(vn.value[i]-vValue(e.midpointCoord())));
            require(std::isfinite(error)&&error<=1e-11,"thin V constant/identifiable normal-linear transfer failed");maximumError=std::max(maximumError,error);}
        Probe geometry(mesh);measure(mesh,geometry);
        FaceSegmentGraph graph;graph.build(mesh);for(const auto& face:graph.segments())if(face.L>=0){
            const int jump=std::abs(mesh.cells(face.K).level()-mesh.cells(face.L).level());require(jump<=1,"thin mesh lost native 2:1");hanging+=jump==1;}
    };
    const auto settle=[&](bool coarsen){
        int sweeps=0;while(true){LAMP::Mesh::NativeTargetLevelField target(mesh,0.),gradient(mesh,0.);
            for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);const auto p=c.centerCoord(),h=c.size();
                const auto intersects=[&](Vector2D lo,Vector2D hi){return std::max(p.x-h.x/2,lo.x)<std::min(p.x+h.x/2,hi.x)&&std::max(p.y-h.y/2,lo.y)<std::min(p.y+h.y/2,hi.y);};
                target.value[i]=coarsen?0:intersects(peakLower,peakUpper)?2:intersects(rootLower,rootUpper)?1:0;}
            if(!mesh.mayAdaptToTargets(target))break;require(sweeps<mesh.maxLevel()+2,"thin native target failed original finite settle bound");
            std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
            std::vector<ThinHProbe*> horizontal{&hc,&hn};std::vector<ThinVProbe*> vertical{&vc,&vn};
            std::vector<LAMP::Mesh::NativeTargetLevelField*> cells{&target};
            mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&target);++sweeps;check();
        }
        require(sweeps>0,"thin mesh did not perform actual requested topology change");return sweeps;
    };
    check();const int refine=settle(false);const int refined=mesh.cellIndexNum();require(refined>background,"thin refinement kept background topology");
    const int coarsen=settle(true);require(mesh.cellIndexNum()==background,"thin coarsening did not recover background");
    const int repeat=settle(false);require(mesh.cellIndexNum()==refined,"thin repeated refinement changed real leaf count");
    require(hanging>0,"thin fixture did not exercise actual native 2:1 faces");
    if(nx==1)require(hn.twoPointNormal>0&&hn.onePointBoundary>0,"one-column fixture missed true H two/one-point source stencils");
    if(ny==1)require(vn.twoPointNormal>0&&vn.onePointBoundary>0,"one-row fixture missed true V two/one-point source stencils");
    std::cout<<"THIN_MESH_PASS nx="<<nx<<" ny="<<ny<<" leaves="<<background<<"->"<<refined<<"->"<<background<<"->"<<mesh.cellIndexNum()
        <<" sweeps="<<refine<<'/'<<coarsen<<'/'<<repeat<<" H2/H1="<<hn.twoPointNormal<<'/'<<hn.onePointBoundary
        <<" V2/V1="<<vn.twoPointNormal<<'/'<<vn.onePointBoundary<<" hanging_visits="<<hanging<<" max_field_error="<<maximumError<<'\n';
}
void thinMeshTransferCases(){thinMeshTransfer(1,1);thinMeshTransfer(1,3);thinMeshTransfer(3,1);}
} // namespace

namespace
{
void refinementPrerequisiteCases()
{
    for (const int requestedLevel : {4,6})
    {
        DictionaryMesh input;input.minCoord={0,0};input.maxCoord={4,4};input.Nx=input.Ny=4;
        input.maxCellNum=32768;input.maxInitialIterations=0;input.maxLevel_local={requestedLevel};
        QuadtreeMesh original(input);const NativeMeshIdentity originalIdentity(original);
        const auto inOriginalBox=[](const Cell& c){const auto x=c.centerCoord(),h=c.size();
            return std::max(x.x-h.x/2,1.)<std::min(x.x+h.x/2,2.)
                && std::max(x.y-h.y/2,1.)<std::min(x.y+h.y/2,2.);};
        const auto fill=[&](const QuadtreeMesh& mesh,auto& target){
            for(int i=0;i<mesh.cellIndexNum();++i)
                target.value[i]=inOriginalBox(mesh.cells(mesh.cellIndexToLabel[i]))?requestedLevel:0.;};
        std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        // Real old target-only path: the four L1 children all touch a coarse
        // base neighbor. False preflight here does not prove original-box coverage.
        auto legacy=original.clone();int oldSweeps=0;
        for(;;){
            LAMP::Mesh::NativeTargetLevelField targets(*legacy,0.),gradient(*legacy,0.);fill(*legacy,targets);
            if(!legacy->mayAdaptToTargets(targets))break;
            require(oldSweeps<requestedLevel+2,"old prerequisite witness exceeded original level bound");
            std::vector<LAMP::Mesh::NativeTargetLevelField*> fields{&targets};
            legacy->adapt(levelSets,gradient,vertices,horizontal,vertical,fields,&targets);++oldSweeps;
        }
        int oldPeak=0;for(int i=0;i<legacy->cellIndexNum();++i){const auto& c=legacy->cells(legacy->cellIndexToLabel[i]);
            if(inOriginalBox(c))oldPeak=std::max(oldPeak,c.level());}
        require(oldSweeps>0&&oldPeak==1&&oldPeak<requestedLevel,
            "old real narrow-ROI path did not witness false-settled incomplete target");
        auto fixed=original.clone();int sweeps=0;std::size_t totalRaised=0;
        LAMP::Mesh::NativeTargetLevelField constant(*fixed,2.);
        for(;;){
            LAMP::Mesh::NativeTargetLevelField targets(*fixed,0.),gradient(*fixed,0.);fill(*fixed,targets);
            std::vector<double> before(fixed->cellIndexNum());
            for(int i=0;i<fixed->cellIndexNum();++i)before[i]=targets.value[i];
            const NativeMeshIdentity beforePlan(*fixed);
            totalRaised+=fixed->completeRefinementPrerequisites(targets);
            beforePlan.validate(*fixed);
            for(int i=0;i<fixed->cellIndexNum();++i)
                require(targets.value[i]>=before[i],"prerequisite completion lowered an original target");
            std::vector<double> once(fixed->cellIndexNum());for(int i=0;i<fixed->cellIndexNum();++i)once[i]=targets.value[i];
            require(fixed->completeRefinementPrerequisites(targets)==0,"prerequisite closure is not idempotent");
            for(int i=0;i<fixed->cellIndexNum();++i)require(targets.value[i]==once[i],"repeated closure altered target bits");
            beforePlan.validate(*fixed);
            if(!fixed->mayAdaptToTargets(targets)){
                for(int i=0;i<fixed->cellIndexNum();++i){const auto& c=fixed->cells(fixed->cellIndexToLabel[i]);
                    if(inOriginalBox(c))require(c.level()==requestedLevel,"false preflight still hides unmet original physical ROI");}
                break;
            }
            require(sweeps<requestedLevel+2,"completed target did not settle within the original level bound");
            std::vector<LAMP::Mesh::NativeTargetLevelField*> fields{&targets,&constant};
            fixed->adapt(levelSets,gradient,vertices,horizontal,vertical,fields,&targets);++sweeps;
            for(int i=0;i<constant.size();++i)require(constant.value[i]==2.,"constant discrete field changed during actual adapt");
        }
        std::vector<int> levels(requestedLevel+1);int covered=0;
        for(int i=0;i<fixed->cellIndexNum();++i){const auto& c=fixed->cells(fixed->cellIndexToLabel[i]);
            require(c.level()>=0&&c.level()<=requestedLevel,"prerequisite target escaped allocated levels");++levels[c.level()];
            if(inOriginalBox(c)){require(c.level()==requestedLevel,"original-box level not actually attained");++covered;}}
        require(totalRaised>0&&covered==(1<<(2*requestedLevel))&&levels[0]>0,
            "narrow ROI lost requested coverage or expanded to uniform refinement");
        Probe actualGeometry(*fixed);measure(*fixed,actualGeometry); // original 2:1/area/face-closure gates
        const NativeMeshIdentity finalIdentity(*fixed);
        LAMP::Mesh::NativeTargetLevelField noop(*fixed,0.);
        for(int i=0;i<fixed->cellIndexNum();++i)noop.value[i]=fixed->cells(fixed->cellIndexToLabel[i]).level();
        require(fixed->completeRefinementPrerequisites(noop)==0&&!fixed->mayAdaptToTargets(noop),"original current-level no-op changed");
        finalIdentity.validate(*fixed);originalIdentity.validate(original);
        const auto rejects=[&](auto& field){bool rejected=false;
            try{fixed->completeRefinementPrerequisites(field);}catch(const std::invalid_argument&){rejected=true;}
            require(rejected,"invalid/stale prerequisite target was accepted");finalIdentity.validate(*fixed);};
        Field<double> shortField(fixed->cellIndexNum()-1,0.);rejects(shortField);
        CellField<double> bad(*fixed,0.);bad.value[0]=requestedLevel;
        bad.value[bad.size()-1]=std::numeric_limits<double>::quiet_NaN();rejects(bad);
        require(bad.value[0]==requestedLevel&&std::isnan(bad.value[bad.size()-1]),"nonfinite rejection partly wrote targets");
        bad.value[bad.size()-1]=std::numeric_limits<double>::infinity();rejects(bad);
        bad.value[bad.size()-1]=static_cast<double>(std::numeric_limits<int>::max())+1.;rejects(bad);
        bad.value[bad.size()-1]=0.;const int oldLabel=bad.label(0);bad.setLabel(0,-1);rejects(bad);bad.setLabel(0,oldLabel);
        finalIdentity.validate(*fixed);
        std::cout<<"REFINEMENT_PREREQUISITES_PASS target="<<requestedLevel<<" old_peak="<<oldPeak
            <<" actual_ROI_leaves="<<covered<<" source_leaves="<<original.cellIndexNum()<<" target_leaves="<<fixed->cellIndexNum()
            <<" sweeps="<<sweeps<<" raised="<<totalRaised<<" coarse_base_leaves="<<levels[0]<<" guards=5\n";
    }
}
}

namespace
{
// Explicit storage opt-in only: the legacy protected setter still clamps to one.
struct ZeroLsSetterWitness : BackgroundMesh
{
    explicit ZeroLsSetterWitness(const DictionaryMesh& input):BackgroundMesh(input){}
    void legacySet(int count){setLevelSetCount(count);ensureLevelSetStorage();}
};
struct ZeroLsSlots {int parents=0,free=0;};
ZeroLsSlots checkLsWidth(const BackgroundMesh& mesh,int count)
{
    require(mesh.levelSetCount()==count,"level-set count differs from requested storage width");
    const auto width=static_cast<std::size_t>(count);ZeroLsSlots slots;
    for(int label=0;label<mesh.vertLabelCapacity();++label){const auto& v=mesh.verts(label);
        require(v.onCutEdge.size()==width&&v.onCutCell.size()==width&&v.nearestNodeCoord.size()==width
            &&v.nearestNodeLabel.size()==width,"vertex LS width differs in allocated storage");slots.free+=v.label()==-1;}
    for(int label=0;label<mesh.hEdgeLabelCapacity();++label){const auto& e=mesh.hEdges(label);
        require(e.isCut.size()==width,"horizontal edge LS width differs in allocated storage");slots.free+=e.label()==-1;}
    for(int label=0;label<mesh.vEdgeLabelCapacity();++label){const auto& e=mesh.vEdges(label);
        require(e.isCut.size()==width,"vertical edge LS width differs in allocated storage");slots.free+=e.label()==-1;}
    for(int label=0;label<mesh.cellLabelCapacity();++label){const auto& c=mesh.cells(label);
        require(c.isCut.size()==width,"cell LS width differs in allocated storage");slots.free+=c.label()==-1;
        if(c.label()!=-1&&c.son[0]>=0){++slots.parents;
            require(std::find(mesh.cellIndexToLabel,mesh.cellIndexToLabel+mesh.cellIndexNum(),label)
                ==mesh.cellIndexToLabel+mesh.cellIndexNum(),"occupied parent appeared as native leaf");}}
    return slots;
}
std::vector<char> zeroLsCache(const QuadtreeMesh& mesh,const fs::path& folder)
{
    fs::create_directories(folder);mesh.writeCache(folder.generic_string()+"/");
    return bytes(folder/"mesh_cache.bin");
}
int zeroLsEncodedCount(const std::vector<char>& data)
{
    static_assert(sizeof(int)==4&&sizeof(std::uint32_t)==4,"QMC2 native four-byte count contract");
    require(data.size()>=12&&std::string(data.data(),4)=="QMC2","missing true QMC2 header");
    std::uint32_t version=0;int count=-1;
    std::memcpy(&version,data.data()+4,sizeof(version));std::memcpy(&count,data.data()+8,sizeof(count));
    require(version==2,"zero-LS test did not serialize QMC2 version2");return count;
}
void zeroLsWrite(const fs::path& folder,const std::vector<char>& data)
{
    fs::create_directories(folder);std::ofstream out(folder/"mesh_cache.bin",std::ios::binary);
    require(bool(out),"cannot open zero-LS bad-cache fixture");
    out.write(data.data(),static_cast<std::streamsize>(data.size()));out.close();
    require(bool(out),"cannot finish zero-LS bad-cache fixture");
}
// The requested ROI is frozen in background coordinates, independently of any
// transferred target. Empty physical LS input is used for every real adapt().
int zeroLsSettle(QuadtreeMesh& mesh,bool coarsen)
{
    LAMP::Mesh::NativeTargetLevelField constant(mesh,2.);
    VertField<double> nodal(mesh,0.);
    for(int i=0;i<mesh.vertIndexNum();++i)nodal.value[i]=reference(mesh,mesh.verts(mesh.vertIndexToLabel[i]).coord());
    int sweeps=0;
    for(;;){
        LAMP::Mesh::NativeTargetLevelField targets(mesh,0.),gradient(mesh,0.);
        for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);const auto x=c.centerCoord(),h=c.size();
            const bool overlap=std::max(x.x-h.x/2,1.)<std::min(x.x+h.x/2,2.)
                &&std::max(x.y-h.y/2,1.)<std::min(x.y+h.y/2,2.);
            targets.value[i]=!coarsen&&overlap?3.:0.;}
        mesh.completeRefinementPrerequisites(targets);
        if(!mesh.mayAdaptToTargets(targets))break;
        require(sweeps<mesh.maxLevel()+2,"zero-LS native adaptation exceeded original finite level bound");
        std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices{&nodal};
        std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
        std::vector<LAMP::Mesh::NativeTargetLevelField*> fields{&targets,&constant};
        mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,fields,&targets);++sweeps;
        checkLsWidth(mesh,0);
        require(constant.size()==mesh.cellIndexNum()&&nodal.size()==mesh.vertIndexNum(),"zero-LS transferred field dimensions differ");
        for(int i=0;i<constant.size();++i)
            require(constant.label(i)==mesh.cellIndexToLabel[i]&&constant.value[i]==2.,"zero-LS actual constant native transfer changed");
        for(int i=0;i<nodal.size();++i){const int label=mesh.vertIndexToLabel[i];const double error=std::abs(nodal.value[i]-reference(mesh,mesh.verts(label).coord()));
            require(nodal.label(i)==label&&std::isfinite(error)&&error<=1e-11,"zero-LS actual bilinear vertex transfer changed");}
        Probe geometry(mesh);measure(mesh,geometry);
    }
    require(sweeps>0,"zero-LS path did not perform an actual topology change");
    std::array<int,4> levels{};int roi=0;
    for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);const auto x=c.centerCoord(),h=c.size();
        require(c.level()>=0&&c.level()<=3,"zero-LS mesh escaped allocated levels");++levels[c.level()];
        if(!coarsen&&std::max(x.x-h.x/2,1.)<std::min(x.x+h.x/2,2.)&&std::max(x.y-h.y/2,1.)<std::min(x.y+h.y/2,2.)){
            require(c.level()==3,"zero-LS preflight falsely settled below the original ROI level");++roi;}}
    if(coarsen)require(mesh.cellIndexNum()==mesh.Nx*mesh.Ny&&levels[0]==mesh.Nx*mesh.Ny,"zero-LS coarsening did not recover background");
    else require(roi==64&&levels[0]>0&&levels[1]>0&&levels[2]>0&&levels[3]==64,"zero-LS fixture lost actual mixed L0/1/2/3 topology");
    return sweeps;
}
void zeroLsCases(const fs::path& scratch)
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={4,4};input.Nx=input.Ny=4;
    input.maxCellNum=1024;input.maxInitialIterations=0;
    // Equal auxiliary vector lengths make the LS1-minus-LS0 file length an
    // independent exact count of the seven serialized cut-vector categories.
    input.minLevel_neg={0};input.minLevel_pos={0};input.maxLevel_local={3};
    require(input.levelSetCount==1,"dictionary default LS count changed");
    QuadtreeMesh oldDefault(input);checkLsWidth(oldDefault,1);
    auto negativeInput=input;negativeInput.levelSetCount=-2;
    QuadtreeMesh negative(negativeInput);checkLsWidth(negative,1);
    auto zeroInput=input;zeroInput.levelSetCount=0;
    ZeroLsSetterWitness setter(zeroInput);checkLsWidth(setter,0);
    setter.legacySet(0);checkLsWidth(setter,1);setter.legacySet(-2);checkLsWidth(setter,1);
    QuadtreeMesh mesh(zeroInput);checkLsWidth(mesh,0);
    const NativeMeshIdentity background(mesh);background.validate(oldDefault);background.validate(negative);
    const auto zeroBackground=zeroLsCache(mesh,scratch/"zero_background");
    const auto positiveBackground=zeroLsCache(oldDefault,scratch/"positive_background");
    const std::size_t expectedExtra=static_cast<std::size_t>(mesh.vertLabelCapacity())
        *(2*sizeof(unsigned char)+sizeof(Vector2D)+sizeof(int))
        +static_cast<std::size_t>(mesh.hEdgeLabelCapacity()+mesh.vEdgeLabelCapacity()+mesh.cellLabelCapacity())*sizeof(unsigned char);
    require(zeroLsEncodedCount(zeroBackground)==0&&zeroLsEncodedCount(positiveBackground)==1
        &&positiveBackground.size()==zeroBackground.size()+expectedExtra,"zero-LS QMC2 did not remove the exact cut-array payload");
    auto twoInput=input;twoInput.levelSetCount=2;twoInput.minLevel_neg={0,0};twoInput.minLevel_pos={0,0};twoInput.maxLevel_local={3,3};
    QuadtreeMesh two(twoInput);checkLsWidth(two,2);
    const auto twoBytes=zeroLsCache(two,scratch/"positive_two");
    require(zeroLsEncodedCount(twoBytes)==2,"positive two-LS QMC2 encoding changed");
    QuadtreeMesh positiveRestored(zeroInput);positiveRestored.readCache((scratch/"positive_two").generic_string()+"/");
    checkLsWidth(positiveRestored,2);
    require(zeroLsCache(positiveRestored,scratch/"positive_two_roundtrip")==twoBytes,"positive two-LS roundtrip changed native records");

    const int refine=zeroLsSettle(mesh,false),refinedLeaves=mesh.cellIndexNum();
    const auto mixedSlots=checkLsWidth(mesh,0);
    require(mixedSlots.parents>0&&mixedSlots.free>0,"zero-LS mixed fixture missed occupied parents or allocated free slots");
    const NativeMeshIdentity mixed(mesh);auto candidate=mesh.clone();mixed.validate(*candidate);
    checkLsWidth(*candidate,0); // Includes nonindexed parents and the entire free capacity.
    const auto mixedBytes=zeroLsCache(mesh,scratch/"mixed");
    require(zeroLsEncodedCount(mixedBytes)==0,"mixed zero-LS cache silently became one-LS");
    QuadtreeMesh restored(twoInput);restored.readCache((scratch/"mixed").generic_string()+"/");
    checkLsWidth(restored,0);mixed.validate(restored);
    require(zeroLsCache(restored,scratch/"mixed_roundtrip")==mixedBytes,"zero-LS QMC2 full-capacity roundtrip changed bytes");
    const auto reject=[&](std::vector<char> bad,const char* name){
        const auto folder=scratch/name;zeroLsWrite(folder,bad);bool rejected=false;
        try{restored.readCache(folder.generic_string()+"/");}catch(const std::runtime_error&){rejected=true;}
        require(rejected,std::string("zero-LS invalid QMC2 accepted: ")+name);
        require(bytes(folder/"mesh_cache.bin")==bad,"rejected QMC2 input file changed");
        mixed.validate(restored);checkLsWidth(restored,0);
        require(zeroLsCache(restored,scratch/(std::string(name)+"_state"))==mixedBytes,"rejected QMC2 partially committed native storage/metadata");
    };
    auto bad=mixedBytes;const int negativeCount=-1;std::memcpy(bad.data()+8,&negativeCount,sizeof(negativeCount));reject(bad,"negative_count");
    bad=mixedBytes;const int wrongWidth=1;std::memcpy(bad.data()+8,&wrongWidth,sizeof(wrongWidth));reject(bad,"false_positive_width");
    bad=mixedBytes;bad.pop_back();reject(bad,"truncated_tail");

    const int coarsen=zeroLsSettle(*candidate,true);zeroLsSettle(restored,true);
    const auto coarseSlots=checkLsWidth(*candidate,0);
    require(coarseSlots.free>mixedSlots.free,"actual coarsening did not create retired/free storage");
    auto holes=candidate->clone();checkLsWidth(*holes,0);
    const NativeMeshIdentity coarse(*candidate);coarse.validate(*holes);coarse.validate(restored);
    // A cloned free pool is exercised by real subsequent allocation, rather
    // than comparing unspecified retired-record contents in a cache file.
    const int repeat=zeroLsSettle(*candidate,false);zeroLsSettle(*holes,false);zeroLsSettle(restored,false);
    mixed.validate(*candidate);mixed.validate(*holes);mixed.validate(restored);mixed.validate(mesh);
    checkLsWidth(*candidate,0);checkLsWidth(*holes,0);checkLsWidth(restored,0);
    Probe original(mesh),repeated(*holes);measure(mesh,original);measure(*holes,repeated);
    const auto first=TecplotMeshOutput::nodalZone(mesh,{"bilinear_probe"},original.nodal);
    const auto last=TecplotMeshOutput::nodalZone(*holes,{"bilinear_probe"},repeated.nodal);
    require(first.connectivity==last.connectivity&&first.fields[0].values==last.fields[0].values,
        "zero-LS refine/coarsen/clone/refine changed physical nodal topology or values");
    std::cout<<"ZERO_LS_PASS background="<<input.Nx*input.Ny<<" mixed="<<refinedLeaves
        <<" sweeps="<<refine<<'/'<<coarsen<<'/'<<repeat<<" occupied_parents="<<mixedSlots.parents
        <<" free_slots="<<mixedSlots.free<<"->"<<coarseSlots.free<<" omitted_payload_bytes="<<expectedExtra
        <<" rejected_cache_guards=3 old_default_negative_setter=1 positive_two_roundtrip=1\n";
}
}


namespace NativeTreeWorkspaceContract
{
using Native=LAMP::Mesh::NativeTargetLevelField;
// The reference retains the REAL discrete target transfer, including all
// current ancestry and previous-family/source checks. Only its C++ exact type
// is different, so optional built-in maintenance cannot enter this path.
struct LiteralTarget : Native
{
    explicit LiteralTarget(const QuadtreeMesh& mesh,double value=0.):Native(mesh,value){}
    void updateWithMesh(const QuadtreeMesh& mesh){Native::updateWithMesh(mesh);}
};
std::size_t pairs=0,slotRecords=0,fieldValues=0,fastCalls=0,unsupportedCalls=0,
    fallbackCalls=0,negativePairs=0,flagsCompared=0,changedCalls=0,noWorkCalls=0;
struct Floating
{
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    Floating(){require(std::fegetenv(&native)==0,"tree oracle cannot capture native FP state");}
    void restore()const{require(std::fesetenv(&native)==0,"tree oracle cannot restore native FP state");
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
void sameFloating(const Floating& a,const Floating& b)
{
    require(std::memcmp(&a.native,&b.native,sizeof a.native)==0&&a.error==b.error,
        "tree workspace changed original full native FP state/errno");
#if defined(_M_X64) || defined(__SSE2__)
    require(a.csr==b.csr,"tree workspace changed original raw MXCSR");
#endif
    ++flagsCompared;
}
bool sameDouble(double a,double b){return std::memcmp(&a,&b,sizeof a)==0;}
template<class A,class B>void sameField(const A& a,const B& b)
{
    require(a.size()==b.size(),"tree original/native target field size differs");
    for(int i=0;i<a.size();++i){require(a.label(i)==b.label(i)&&sameDouble(a.value[i],b.value[i]),
        "tree original/native target labels or value bits differ");++fieldValues;}
}
std::string record(const Vert& v,int ls){std::ostringstream out(std::ios::out|std::ios::binary);MeshCacheIO::writeVert(out,v,ls);return out.str();}
std::string record(const Edge& e,int ls){std::ostringstream out(std::ios::out|std::ios::binary);MeshCacheIO::writeEdge(out,e,ls);return out.str();}
std::string record(const Cell& c,int ls){std::ostringstream out(std::ios::out|std::ios::binary);MeshCacheIO::writeCell(out,c,ls);return out.str();}
struct CachePrefix
{
    std::array<int,4> count{},capacity{},ptr{};std::array<std::vector<int>,4> labels;
};
CachePrefix prefix(const std::vector<char>& data)
{
    std::size_t offset=8;
    const auto get=[&](){int result=0;require(offset+sizeof result<=data.size(),"tree QMC2 prefix truncated");
        std::memcpy(&result,data.data()+offset,sizeof result);offset+=sizeof result;return result;};
    require(data.size()>=8&&std::string(data.data(),4)=="QMC2","tree oracle lacks QMC2 native header");
    require(get()>=0,"tree negative actual LS count");CachePrefix r;
    for(auto& n:r.count)n=get();
    for(const auto n:r.count){require(n>0,"tree empty actual indexed pool");offset+=static_cast<std::size_t>(n)*sizeof(int);}
    require(get()>=0,"tree invalid actual max level");
    for(auto& n:r.capacity)n=get();for(auto& n:r.ptr)n=get();
    for(int k=0;k<4;++k){require(r.capacity[k]>0&&r.ptr[k]>=0&&r.ptr[k]<=r.capacity[k],"tree private pool limits differ");
        r.labels[k].reserve(r.capacity[k]);for(int i=0;i<r.capacity[k];++i)r.labels[k].push_back(get());}
    return r;
}
template<class Entity>void checkPrefix(const CachePrefix& p,int k,Entity entity)
{
    int free=0,last=-1;
    for(int i=0;i<p.capacity[k];++i){const auto& t=entity(i);
        if(t.label()==-1)require(p.labels[k][free++]==i,"tree complete ascending free prefix differs from actual ALL-slot occupancy");
        else{require(t.label()==i,"tree occupied nonindexed slot lost label identity");last=i;}}
    require(last<p.ptr[k],"tree actual live slot escaped original private ptr");
    // This is the writer's canonical serialized suffix, NOT a claim to have
    // inspected private unused in-memory label-space tail bytes.
    for(int i=free;i<p.capacity[k];++i)require(p.labels[k][i]==0,"tree QMC2 canonical unused suffix changed");
}
void sameMesh(const QuadtreeMesh& a,const QuadtreeMesh& b,const fs::path& scratch)
{
    require(a.levelSetCount()==b.levelSetCount()&&a.maxLevel()==b.maxLevel(),"tree LS/max-level metadata differs");
    const std::array<int,4> ac{{a.vertLabelCapacity(),a.hEdgeLabelCapacity(),a.vEdgeLabelCapacity(),a.cellLabelCapacity()}},
        bc{{b.vertLabelCapacity(),b.hEdgeLabelCapacity(),b.vEdgeLabelCapacity(),b.cellLabelCapacity()}};
    require(ac==bc,"tree allocated four-pool capacities differ");
    const std::array<int,4> an{{a.vertIndexNum(),a.hEdgeIndexNum(),a.vEdgeIndexNum(),a.cellIndexNum()}},
        bn{{b.vertIndexNum(),b.hEdgeIndexNum(),b.vEdgeIndexNum(),b.cellIndexNum()}};
    require(an==bn,"tree current dense four-pool counts differ");
    const std::array<const int*,4> ai{{a.vertIndexToLabel,a.hEdgeIndexToLabel,a.vEdgeIndexToLabel,a.cellIndexToLabel}},
        bi{{b.vertIndexToLabel,b.hEdgeIndexToLabel,b.vEdgeIndexToLabel,b.cellIndexToLabel}};
    for(int k=0;k<4;++k)for(int i=0;i<an[k];++i)require(ai[k][i]==bi[k][i],"tree complete dense ascending index arrays differ");
    for(int i=0;i<ac[0];++i){require(record(a.verts(i),a.levelSetCount())==record(b.verts(i),b.levelSetCount()),"tree vertex record differs, including live parents/free capacity");++slotRecords;}
    for(int i=0;i<ac[1];++i){require(record(a.hEdges(i),a.levelSetCount())==record(b.hEdges(i),b.levelSetCount()),"tree H record/ordered transfer history differs over ALL capacity");++slotRecords;}
    for(int i=0;i<ac[2];++i){require(record(a.vEdges(i),a.levelSetCount())==record(b.vEdges(i),b.levelSetCount()),"tree V record/ordered transfer history differs over ALL capacity");++slotRecords;}
    for(int i=0;i<ac[3];++i){require(record(a.cells(i),a.levelSetCount())==record(b.cells(i),b.levelSetCount()),"tree cell record/ordered transfer history differs, including unindexed parents/free capacity");++slotRecords;}
    const auto ab=zeroLsCache(a,scratch/"actual"),bb=zeroLsCache(b,scratch/"literal");
    require(ab==bb,"tree full QMC2 bytes differ (ptr/free prefix/per-level/private metadata included)");
    const auto p=prefix(ab);require(p.count==an&&p.capacity==ac,"tree serialized pool metadata disagrees with live objects");
    checkPrefix(p,0,[&](int i)->const Vert&{return a.verts(i);});checkPrefix(p,1,[&](int i)->const Edge&{return a.hEdges(i);});
    checkPrefix(p,2,[&](int i)->const Edge&{return a.vEdges(i);});checkPrefix(p,3,[&](int i)->const Cell&{return a.cells(i);});
    require(a.negFaceVerts==b.negFaceVerts&&a.posFaceVerts==b.posFaceVerts,"tree phase boundary vertex lists differ");
    NativeMeshIdentity identity(a);identity.validate(b);checkLsWidth(a,a.levelSetCount());checkLsWidth(b,b.levelSetCount());
    Probe pa(const_cast<QuadtreeMesh&>(a)),pb(const_cast<QuadtreeMesh&>(b));measure(a,pa);measure(b,pb);
    ++pairs;
}
void report(const QuadtreeMesh& mesh,int expected,std::size_t entryConstructedTransferSlots=0)
{
    const auto r=mesh.nativeAdaptationMaintenanceReport();
    require(r.workspaceLimitBytes==64u*1024u*1024u,"tree workspace public budget differs");
    if(expected==0){require(!r.requested&&!r.eligible&&!r.actuallyUsed&&!r.fallback&&r.journaledSlots==0&&r.changedMembershipSlots==0
        &&r.clearedTransferSlots==0&&r.exportedIndexEntries==0&&r.exportedFreeEntries==0&&r.workspaceBytes==0,
        "tree opaque/generic path fabricated workspace work");++unsupportedCalls;return;}
    if(expected==2){require(r.requested&&!r.eligible&&!r.actuallyUsed&&r.fallback&&r.workspaceBytes==0,
        "tree inconsistent public binding did not use unchanged literal fallback");++fallbackCalls;return;}
    require(r.requested&&r.eligible&&r.actuallyUsed&&!r.fallback&&r.workspaceBytes>0&&r.workspaceBytes<=r.workspaceLimitBytes,
        "tree real small native call failed to consume qualified call-local workspace");
    require(entryConstructedTransferSlots>0&&r.clearedTransferSlots>=entryConstructedTransferSlots,
        "tree first entry-constructed-slot transfer cleanup was skipped");
    require(r.exportedIndexEntries>0&&r.exportedFreeEntries>0,"tree dense indices/complete free prefix were not actually exported");++fastCalls;
}
struct Desired
{
    int mask=0;Vector2D min,max;
    double operator()(const Cell& c)const
    {
        const auto x=c.centerCoord(),h=c.size();const double dx=(max.x-min.x)/4,dy=(max.y-min.y)/4;
        const auto overlap=[&](int bx,int by){const Vector2D lo{min.x+bx*dx,min.y+by*dy},hi{lo.x+dx,lo.y+dy};
            return std::max(x.x-h.x/2,lo.x)<std::min(x.x+h.x/2,hi.x)&&std::max(x.y-h.y/2,lo.y)<std::min(x.y+h.y/2,hi.y);};
        return ((mask&1)&&overlap(0,0))||((mask&2)&&overlap(3,3))?3.:0.;
    }
};
template<class T>void fill(const QuadtreeMesh& mesh,T& t,const Desired& desired,bool noop)
{for(int i=0;i<mesh.cellIndexNum();++i)t.value[i]=noop?mesh.cells(mesh.cellIndexToLabel[i]).level():desired(mesh.cells(mesh.cellIndexToLabel[i]));}
template<class T>void call(QuadtreeMesh& mesh,T& target,T& constant,T& gradient,VertField<double>* nodal=nullptr)
{
    std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;if(nodal)vertices.push_back(nodal);
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
    std::vector<T*> fields{&target,&constant};mesh.adapt(levelSets,gradient,vertices,horizontal,vertical,fields,&target);
}
bool pair(QuadtreeMesh& actual,QuadtreeMesh& literal,const Desired& desired,const fs::path& scratch,
    bool noop=false,bool corruptIndex=false,bool opaqueVertex=false,bool force=false)
{
    Native a(actual,0.),ac(actual,2.),ag(actual,0.);LiteralTarget b(literal,0.),bc(literal,2.),bg(literal,0.);
    fill(actual,a,desired,noop);fill(literal,b,desired,noop);
    const auto raisedA=actual.completeRefinementPrerequisites(a),raisedB=literal.completeRefinementPrerequisites(b);
    require(raisedA==raisedB,"tree actual-neighbor F10 closure changed between real target modes");sameField(a,b);
    const bool neededA=actual.mayAdaptToTargets(a),neededB=literal.mayAdaptToTargets(b);require(neededA==neededB,"tree original preflight differs");
    if(!neededA&&!force)return false;
    std::unique_ptr<VertField<double>> av,bv;
    if(opaqueVertex){av=std::make_unique<VertField<double>>(actual,0.);bv=std::make_unique<VertField<double>>(literal,0.);
        for(int i=0;i<actual.vertIndexNum();++i)av->value[i]=reference(actual,actual.verts(actual.vertIndexToLabel[i]).coord());
        for(int i=0;i<literal.vertIndexNum();++i)bv->value[i]=reference(literal,literal.verts(literal.vertIndexToLabel[i]).coord());}
    const NativeMeshIdentity before(actual);const int oldLeaves=actual.cellIndexNum();
    if(corruptIndex){actual.cells(actual.cellIndexToLabel[0]).setIndex(-7);literal.cells(literal.cellIndexToLabel[0]).setIndex(-7);}
    // Actual object stock is sampled BEFORE adapt; lazy default tail slots are not objects to clear.
    const auto entryPools=actual.entityPoolStorageReport();
    const auto entryConstructedTransferSlots=entryPools[1].constructedSlots+entryPools[2].constructedSlots+entryPools[3].constructedSlots;
    const Floating entry;call(actual,a,ac,ag,av.get());const Floating afterActual;entry.restore();call(literal,b,bc,bg,bv.get());const Floating afterLiteral;
    sameFloating(afterActual,afterLiteral);sameField(a,b);sameField(ac,bc);sameField(ag,bg);
    for(int i=0;i<ac.size();++i)require(ac.value[i]==2.,"tree NativeTarget original constant source transfer changed");
    if(opaqueVertex){sameField(*av,*bv);for(int i=0;i<av->size();++i){const auto error=std::abs(av->value[i]-reference(actual,actual.verts(av->label(i)).coord()));
        require(std::isfinite(error)&&error<=1e-11,"tree opaque-path original bilinear vertex transfer failed");}}
    report(actual,(opaqueVertex||actual.levelSetCount()!=0)?0:corruptIndex?2:1,entryConstructedTransferSlots);report(literal,0);sameMesh(actual,literal,scratch);
    if(noop){before.validate(actual);before.validate(literal);++noWorkCalls;}else if(oldLeaves!=actual.cellIndexNum())++changedCalls;
    return true;
}
void dirtyPublicTransfer(QuadtreeMesh& a,QuadtreeMesh& b)
{
    const auto dirty=[](auto& t){t.setPreIndex(791);t.relatedPreIndex={711,712};t.relatedPreCoord={{.125,.25},{.5,.75}};};
    for(int i=0;i<a.cellLabelCapacity();++i)if(a.cells(i).label()==-1||a.cells(i).son[0]>=0){dirty(a.cells(i));dirty(b.cells(i));}
    for(int i=0;i<a.hEdgeLabelCapacity();++i)if(a.hEdges(i).label()==-1||a.hEdges(i).son[0]>=0){dirty(a.hEdges(i));dirty(b.hEdges(i));}
    for(int i=0;i<a.vEdgeLabelCapacity();++i)if(a.vEdges(i).label()==-1||a.vEdges(i).son[0]>=0){dirty(a.vEdges(i));dirty(b.vEdges(i));}
}
void cleanTransfer(const QuadtreeMesh& mesh)
{
    const auto check=[&](const auto& t){require(t.relatedPreIndex.empty()&&t.relatedPreCoord.empty(),"tree no-work final-level cleanup retained old public metadata");
        const bool leaf=t.label()!=-1&&t.son[0]==-1;require(t.preIndex()==(leaf?t.index():-1),"tree no-work final preIndex differs from original indexed/nonindexed rule");};
    for(int i=0;i<mesh.cellLabelCapacity();++i)check(mesh.cells(i));
    for(int i=0;i<mesh.hEdgeLabelCapacity();++i)check(mesh.hEdges(i));
    for(int i=0;i<mesh.vEdgeLabelCapacity();++i)check(mesh.vEdges(i));
}
template<class F>void invalidPair(F action,const QuadtreeMesh& a,const QuadtreeMesh& b,const fs::path& scratch)
{
    bool rejectedA=false,rejectedB=false;try{action(a,false);}catch(const std::invalid_argument&){rejectedA=true;}
    try{action(b,true);}catch(const std::invalid_argument&){rejectedB=true;}
    require(rejectedA&&rejectedB,"tree literal wrapper bypassed an original target source/F10 guard");sameMesh(a,b,scratch);++negativePairs;
}
void guards(const QuadtreeMesh& a,const QuadtreeMesh& b,const fs::path& scratch)
{
    invalidPair([](const auto& mesh,bool generic){if(generic){LiteralTarget t(mesh);t.setLabel(0,-1);mesh.completeRefinementPrerequisites(t);}
        else{Native t(mesh);t.setLabel(0,-1);mesh.completeRefinementPrerequisites(t);}},a,b,scratch);
    for(double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),static_cast<double>(std::numeric_limits<int>::max())+1.})
        invalidPair([&](const auto& mesh,bool generic){if(generic){LiteralTarget t(mesh);t.value[t.size()-1]=bad;mesh.mayAdaptToTargets(t);}
            else{Native t(mesh);t.value[t.size()-1]=bad;mesh.mayAdaptToTargets(t);}},a,b,scratch);
    invalidPair([](const auto& mesh,bool generic){if(generic){LiteralTarget t(mesh);t.setLabel(0,-1);t.updateWithMesh(mesh);}
        else{Native t(mesh);t.setLabel(0,-1);t.updateWithMesh(mesh);}},a,b,scratch);
    invalidPair([](const auto& mesh,bool generic){if(generic){LiteralTarget t(mesh);t.value[t.size()-1]=.5;t.updateWithMesh(mesh);}
        else{Native t(mesh);t.value[t.size()-1]=.5;t.updateWithMesh(mesh);}},a,b,scratch);
    invalidPair([](const auto& mesh,bool generic){
        const auto loseActualParentSlot=[&](auto& target){int child=-1;
            for(int i=0;i<mesh.cellIndexNum()&&child<0;++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);if(c.parent>=0)child=c.label();}
            require(child>=0,"tree source guard fixture has no actual non-root child");auto& parent=mesh.cells(mesh.cells(child).parent);int slot=-1;
            for(int k=0;k<4;++k)if(parent.son[k]==child)slot=k;
            require(slot>=0,"tree source guard fixture has no actual parent slot");parent.son[slot]=-1;
            try{target.updateWithMesh(mesh);}catch(...){parent.son[slot]=child;throw;}parent.son[slot]=child;
        };
        if(generic){LiteralTarget t(mesh);loseActualParentSlot(t);}else{Native t(mesh);loseActualParentSlot(t);}
    },a,b,scratch);
}
void coverage(const QuadtreeMesh& mesh,const Desired& desired)
{
    int fine=0,coarse=0;
    for(int i=0;i<mesh.cellIndexNum();++i){const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);
        if(desired(c)==3.){require(c.level()==3,"tree F10 false-settled below actual frozen original ROI");++fine;}coarse+=c.level()==0;}
    const int boxes=!!(desired.mask&1)+!!(desired.mask&2);require(fine==64*boxes,"tree actual original ROI leaf coverage differs");
    if(boxes)require(coarse>0,"tree localized targets were replaced by global uniform refinement");
    else require(mesh.cellIndexNum()==16,"tree all-zero coarsening did not recover original background");
}
void scene(const char* name,Vector2D lower,Vector2D upper,const fs::path& scratch)
{
    DictionaryMesh input;input.minCoord=lower;input.maxCoord=upper;input.Nx=input.Ny=4;input.levelSetCount=0;
    input.maxCellNum=2048;input.maxInitialIterations=0;input.maxLevel_local={3};
    QuadtreeMesh source(input);const NativeMeshIdentity sourceIdentity(source);const auto original=zeroLsCache(source,scratch/"accepted_source");
    auto actual=source.clone(),literal=source.clone();
    const auto settle=[&](int mask){const Desired desired{mask,lower,upper};int sweeps=0;
        while(pair(*actual,*literal,desired,scratch/"pair")){require(++sweeps<=source.maxLevel()+2,"tree exact target failed original bounded F10 settle schedule");}
        coverage(*actual,desired);coverage(*literal,desired);return sweeps;};
    require(settle(1)>0,"tree first local refinement did not change real topology");
    int parents=0;for(int i=0;i<actual->cellLabelCapacity();++i)parents+=actual->cells(i).label()!=-1&&actual->cells(i).son[0]>=0;
    require(parents>0&&actual->cellIndexNum()>16,"tree fixture lacks actual unindexed live ancestors");
    // Move the desired region across the domain: true removals plus births,
    // not a positive net-leaf count used as an oracle for changed work.
    require(settle(2)>0,"tree simultaneous distant coarsen/refine was not exercised");
    dirtyPublicTransfer(*actual,*literal);require(pair(*actual,*literal,{2,lower,upper},scratch/"mutated_noop",true,false,false,true),"tree explicit no-work adapt omitted");
    cleanTransfer(*actual);cleanTransfer(*literal);
    // Public mutable dense index is corruptible. The entry must rediscover it;
    // the unchanged generic count() repairs it before NativeTarget consumes it.
    require(pair(*actual,*literal,{2,lower,upper},scratch/"inconsistent_entry",true,true,false,true),"tree inconsistent entry was not exercised");
    auto cloned=actual->clone();report(*cloned,0);NativeMeshIdentity(*actual).validate(*cloned);
    const auto accepted=zeroLsCache(*actual,scratch/"roundtrip_source");QuadtreeMesh restored(input);restored.readCache((scratch/"roundtrip_source").generic_string()+"/");
    require(zeroLsCache(restored,scratch/"roundtrip_restored")==accepted,"tree actual post-workspace QMC2 roundtrip changed complete records");
    guards(*actual,*literal,scratch/"guards");require(settle(0)>0,"tree coarsening did not retire/recycle actual slots");
    require(settle(3)>0,"tree actual recycled-label repeated distant refinement was not exercised");
    // Real external vertex consumers disqualify the specialization and retain
    // the original generic transfer; no trust flag makes them eligible.
    require(pair(*actual,*literal,{0,lower,upper},scratch/"opaque",false,false,true,true),"tree real opaque fallback not entered");
    sourceIdentity.validate(source);require(zeroLsCache(source,scratch/"accepted_source_after")==original,"tree candidate maintenance mutated original accepted storage");
    std::cout<<"NATIVE_TREE_WORKSPACE_SCENE,"<<name<<",leaves,"<<actual->cellIndexNum()<<",occupied_parents,"<<parents<<",all_capacity_and_transfer_QMC2,exact\n";
}
void positiveLsFallback(const fs::path& scratch)
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={4,4};input.Nx=input.Ny=4;input.levelSetCount=1;
    input.maxCellNum=128;input.maxInitialIterations=0;input.maxLevel_local={3};QuadtreeMesh source(input);auto a=source.clone(),b=source.clone();
    require(pair(*a,*b,{0,input.minCoord,input.maxCoord},scratch,true,false,false,true),"tree positive-LS original fallback not executed");
    require(a->levelSetCount()==1&&b->levelSetCount()==1,"tree unsupported-LS call silently changed storage policy");
}
void allocationVeto(const fs::path& scratch)
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={4,4};input.Nx=input.Ny=4;input.levelSetCount=0;
    input.maxCellNum=32;input.maxInitialIterations=0;input.maxLevel_local={3};QuadtreeMesh source(input);auto a=source.clone(),b=source.clone();
    require(pair(*a,*b,{3,input.minCoord,input.maxCoord},scratch,false,false,false,true),"tree original capacity-veto adapt did not execute");
    require(a->cellIndexNum()<128,"tree bounded-allocation fixture unexpectedly completed two fine ROIs");
    const auto p=prefix(zeroLsCache(*a,scratch/"veto_limits"));
    require(p.capacity[0]-p.ptr[0]<5||p.capacity[1]-p.ptr[1]<6||p.capacity[2]-p.ptr[2]<6||p.capacity[3]-p.ptr[3]<4,
        "tree capacity-veto fixture did not reach an original conservative allocation bound");
}
void run(const fs::path& scratch)
{
    scene("XY",{-2,-1},{2,3},scratch/"xy");
    scene("axis_coordinates",{0,-2},{4,2},scratch/"axis");
    scene("annulus_coordinates",{.125,-.25},{4.125,3.75},scratch/"annulus");
    positiveLsFallback(scratch/"positive_ls");allocationVeto(scratch/"veto");
    require(fastCalls>0&&unsupportedCalls>0&&fallbackCalls==3&&negativePairs==21&&flagsCompared>0&&noWorkCalls>=6&&changedCalls>0,
        "tree workspace finite actual scenario/guard coverage incomplete");
    std::cout<<"NATIVE_TREE_WORKSPACE_PASS,mesh_pairs,"<<pairs<<",all_capacity_record_bits,"<<slotRecords<<",field_value_bits,"<<fieldValues
        <<",native_workspace_calls,"<<fastCalls<<",literal_or_opaque_calls,"<<unsupportedCalls<<",actual_binding_fallbacks,"<<fallbackCalls
        <<",paired_original_rejections,"<<negativePairs<<",full_FP_errno_pairs,"<<flagsCompared
        <<",private_unused_tail_dynamic_unavailable,optional_resource_failure_not_injected,no_speed_claim\n";
}
}


namespace EntityPoolStorageContract
{
using NativeTreeWorkspaceContract::Floating;
using NativeTreeWorkspaceContract::LiteralTarget;
using NativeTreeWorkspaceContract::Desired;
using NativeTreeWorkspaceContract::record;
using NativeTreeWorkspaceContract::sameFloating;
using NativeTreeWorkspaceContract::sameMesh;
using Stock=std::array<NativeEntityPoolStorageReport,4>;
std::size_t scenes=0,denseReferences=0,parentRecords=0,writerPairs=0,firstClears=0,
    highLabels=0,concurrentPairs=0,cacheRejections=0,positiveCases=0,environmentPairs=0,initialSources=0,initialWriterPairs=0;
Stock stock(const QuadtreeMesh& mesh)
{
    const Floating before;const auto result=mesh.entityPoolStorageReport();const Floating after;sameFloating(before,after);
    const std::array<std::size_t,4> capacity{{static_cast<std::size_t>(mesh.vertLabelCapacity()),
        static_cast<std::size_t>(mesh.hEdgeLabelCapacity()),static_cast<std::size_t>(mesh.vEdgeLabelCapacity()),
        static_cast<std::size_t>(mesh.cellLabelCapacity())}};
    const std::array<std::size_t,4> sizes{{sizeof(Vert),sizeof(Edge),sizeof(Edge),sizeof(Cell)}};
    for(std::size_t k=0;k<result.size();++k){const auto& r=result[k];
        require(r.logicalCapacity==capacity[k]&&r.constructedSlots<=r.logicalCapacity,
            "pool report does not describe actual logical/constructed storage");
        require(r.entityPayloadBytes==r.constructedSlots*sizes[k],"pool actual entity payload byte count differs");
        if(r.denseMode)require(r.constructedSlots==r.logicalCapacity&&r.publishedChunks==0,
            "pool dense storage fabricated a sparse chunk report");
        else require(r.publishedChunks<=r.constructedSlots,"pool published chunks exceed constructed objects");}
    return result;
}
void sameStock(const Stock& a,const Stock& b)
{
    for(std::size_t k=0;k<a.size();++k)require(a[k].logicalCapacity==b[k].logicalCapacity
        &&a[k].constructedSlots==b[k].constructedSlots&&a[k].publishedChunks==b[k].publishedChunks
        &&a[k].entityPayloadBytes==b[k].entityPayloadBytes&&a[k].denseMode==b[k].denseMode,
        "pool passive operation changed actual entity stock");
}
void dense(const QuadtreeMesh& mesh)
{for(const auto& r:stock(mesh))require(r.denseMode,"pool independent constructor/cache reference is not dense");++denseReferences;}
void sparse(const QuadtreeMesh& mesh)
{for(const auto& r:stock(mesh))require(!r.denseMode&&r.constructedSlots>0&&r.constructedSlots<r.logicalCapacity
    &&r.publishedChunks>0,"pool LS0 clone did not retain untouched default chunks");}
// This is the ORIGINAL clone free-record rule on the independent dense-adopt
// reference below, never a new clone() used as its own cleanup reference.
void resetFreeReference(QuadtreeMesh& mesh)
{
    dense(mesh);const auto width=static_cast<std::size_t>(mesh.levelSetCount());
    const auto reset=[&](auto& value){using Entity=std::decay_t<decltype(value)>;
        if(value.label()==-1){value=Entity{};value.setLabel(-1);value.ensureLevelSetSize(width);}};
    for(int i=0;i<mesh.vertLabelCapacity();++i)reset(mesh.verts(i));
    for(int i=0;i<mesh.hEdgeLabelCapacity();++i)reset(mesh.hEdges(i));
    for(int i=0;i<mesh.vEdgeLabelCapacity();++i)reset(mesh.vEdges(i));
    for(int i=0;i<mesh.cellLabelCapacity();++i)reset(mesh.cells(i));
}
std::vector<char> passiveCache(const QuadtreeMesh& mesh,const fs::path& folder)
{
    const auto before=stock(mesh);const auto result=zeroLsCache(mesh,folder);sameStock(before,stock(mesh));return result;
}
int settle(QuadtreeMesh& actual,QuadtreeMesh& literal,const Desired& desired,const fs::path& scratch)
{
    int sweeps=0;
    while(NativeTreeWorkspaceContract::pair(actual,literal,desired,scratch))
        require(++sweeps<=actual.maxLevel()+2,"pool actual NativeTarget/generic settle exceeded original bound");
    NativeTreeWorkspaceContract::coverage(actual,desired);NativeTreeWorkspaceContract::coverage(literal,desired);return sweeps;
}
void dirtyFree(QuadtreeMesh& mesh)
{
    for(int i=0;i<mesh.vertLabelCapacity();++i)if(mesh.verts(i).label()==-1){auto& t=mesh.verts(i);
        t.setCoord({.125,-.25});t.setPhase(17);t.supCell=73;}
    const auto dirty=[](auto& t){if(t.label()==-1){t.setPreIndex(19);t.relatedPreIndex={2,3};t.relatedPreCoord={{.125,.25},{.5,.75}};}};
    for(int i=0;i<mesh.hEdgeLabelCapacity();++i)dirty(mesh.hEdges(i));
    for(int i=0;i<mesh.vEdgeLabelCapacity();++i)dirty(mesh.vEdges(i));
    for(int i=0;i<mesh.cellLabelCapacity();++i)dirty(mesh.cells(i));
}
void cloneEnvironment(const QuadtreeMesh& source,const std::vector<char>& expected,const fs::path& scratch)
{
    const Floating original;
    const auto one=[&](int mode,const char* name){require(std::fesetround(mode)==0,"pool cannot set directed clone fixture");
        errno=EDOM;const Floating before;auto candidate=source.clone();const Floating after;sameFloating(before,after);
        original.restore();require(passiveCache(*candidate,scratch/name)==expected,"pool clone controls changed original record bits");++environmentPairs;};
    try{one(FE_UPWARD,"upward");
#if defined(_M_X64) || defined(__SSE2__)
        original.restore();_mm_setcsr(_mm_getcsr()|0x8040u);errno=ERANGE;
        const Floating before;auto candidate=source.clone();const Floating after;sameFloating(before,after);
        original.restore();require(passiveCache(*candidate,scratch/"ftz_daz")==expected,"pool FTZ/DAZ clone changed original record bits");++environmentPairs;
#endif
        original.restore();
    }catch(...){original.restore();throw;}
}
void firstMaintenance(QuadtreeMesh& actual,QuadtreeMesh& literal,const fs::path& scratch)
{
    using Native=NativeTreeWorkspaceContract::Native;
    Native a(actual,0.),ac(actual,2.),ag(actual,0.);LiteralTarget b(literal,0.),bc(literal,2.),bg(literal,0.);
    const Desired desired{2,actual.minCoord,actual.maxCoord};
    NativeTreeWorkspaceContract::fill(actual,a,desired,true);NativeTreeWorkspaceContract::fill(literal,b,desired,true);
    require(actual.completeRefinementPrerequisites(a)==literal.completeRefinementPrerequisites(b),
        "pool first-clean original F10 closure differs");NativeTreeWorkspaceContract::sameField(a,b);
    require(!actual.mayAdaptToTargets(a)&&!literal.mayAdaptToTargets(b),"pool first-clean fixture is not a true no-op");
    const NativeMeshIdentity identity(actual);
    // The obligation is frozen before adapt, not inferred from its new report.
    const auto entryStock=stock(actual);const auto first=entryStock[1].constructedSlots+entryStock[2].constructedSlots+entryStock[3].constructedSlots;
    const Floating entry;NativeTreeWorkspaceContract::call(actual,a,ac,ag);const Floating afterActual;
    // Do not visit full capacity before this observation: inspection/cleanup
    // must not publish absent chunks merely to inspect default records.
    sparse(actual);const auto afterStock=stock(actual);sameStock(entryStock,afterStock);
    entry.restore();NativeTreeWorkspaceContract::call(literal,b,bc,bg);const Floating afterLiteral;
    sameFloating(afterActual,afterLiteral);NativeTreeWorkspaceContract::sameField(a,b);
    NativeTreeWorkspaceContract::sameField(ac,bc);NativeTreeWorkspaceContract::sameField(ag,bg);
    for(int i=0;i<ac.size();++i)require(ac.value[i]==2.,"pool first-clean changed original constant transfer");
    NativeTreeWorkspaceContract::report(actual,1,first);NativeTreeWorkspaceContract::report(literal,0);
    require(actual.nativeAdaptationMaintenanceReport().clearedTransferSlots>=first,
        "pool first clear failed real entry constructed-slot obligation");
    sameMesh(actual,literal,scratch);identity.validate(actual);identity.validate(literal);++firstClears;
}
void mutableHighLabels(QuadtreeMesh& actual,QuadtreeMesh& literal,const fs::path& scratch)
{
    for(const int label:{1023,1024,4096}){
        require(label<actual.cellLabelCapacity()&&label<actual.vertLabelCapacity(),"pool high/short-chunk fixture capacity missing");
        const int first=(label/256)*256,last=std::min(first+256,literal.cellLabelCapacity());
        for(int i=first;i<last;++i)require(literal.verts(i).label()==-1&&literal.hEdges(i).label()==-1
            &&literal.vEdges(i).label()==-1&&literal.cells(i).label()==-1,"pool high-label chunk contains an unexpected live record");
        const auto before=stock(actual);
        auto* vertex=&actual.verts(label);auto* horizontal=&actual.hEdges(label);auto* vertical=&actual.vEdges(label);auto* cell=&actual.cells(label);
        require(record(*vertex,0)==record(Vert{},0)&&record(*horizontal,0)==record(Edge{},0)
            &&record(*vertical,0)==record(Edge{},0)&&record(*cell,0)==record(Cell{},0),"pool exposed absent free slot was not a real default entity");
        vertex->setCoord({.25,.5});literal.verts(label).setCoord({.25,.5});
        horizontal->relatedPreIndex={7,8};literal.hEdges(label).relatedPreIndex={7,8};
        vertical->relatedPreCoord={{.5,.75}};literal.vEdges(label).relatedPreCoord={{.5,.75}};
        cell->setPreIndex(23);literal.cells(label).setPreIndex(23);
        require(&actual.verts(label)==vertex&& &actual.hEdges(label)==horizontal&& &actual.vEdges(label)==vertical&& &actual.cells(label)==cell,
            "pool published entity reference moved after access");
        const auto after=stock(actual);for(std::size_t k=0;k<after.size();++k)
            require(after[k].constructedSlots==before[k].constructedSlots+static_cast<std::size_t>(last-first)
                &&after[k].publishedChunks==before[k].publishedChunks+1,"pool fresh high-label chunk/short-tail construction count differs");
        ++highLabels;
    }
    const auto a=passiveCache(actual,scratch/"actual"),b=zeroLsCache(literal,scratch/"dense_literal");
    require(a==b,"pool actual mutable high-label free payload did not serialize like dense original storage");++writerPairs;
}
void concurrentFirstGetter(const QuadtreeMesh& source,const fs::path& scratch)
{
    auto candidate=source.clone();const auto before=stock(*candidate);
    const int label=candidate->vertLabelCapacity()-1;
    require(label==4096&&source.verts(label).label()==-1,"pool concurrent short-block fixture is not actually free");
    std::array<const Vert*,2> pointers{{nullptr,nullptr}};
    std::array<std::exception_ptr,2> failures{};
    std::array<bool,2> passive{{false,false}};
    std::array<std::thread,2> workers;std::atomic<unsigned> ready{0};std::atomic<bool> start{false};
    struct Join {std::array<std::thread,2>& workers;std::atomic<bool>& start;
        ~Join(){start.store(true,std::memory_order_release);for(auto& t:workers)if(t.joinable())t.join();}} join{workers,start};
    // Each thread owns its FP receipt; no assertion or counter is shared.
    for(std::size_t i=0;i<workers.size();++i)workers[i]=std::thread([&,i](){try{
        ready.fetch_add(1,std::memory_order_release);while(!start.load(std::memory_order_acquire))std::this_thread::yield();
        const Floating entry;const auto* value=&candidate->verts(label);const Floating exit;
        bool same=std::memcmp(&entry.native,&exit.native,sizeof entry.native)==0&&entry.error==exit.error;
#if defined(_M_X64) || defined(__SSE2__)
        same=same&&entry.csr==exit.csr;
#endif
        pointers[i]=value;passive[i]=same;
    }catch(...){failures[i]=std::current_exception();}});
    while(ready.load(std::memory_order_acquire)!=workers.size())std::this_thread::yield();start.store(true,std::memory_order_release);
    for(auto& worker:workers)worker.join();
    for(std::size_t i=0;i<workers.size();++i){if(failures[i])std::rethrow_exception(failures[i]);
        require(passive[i],"pool concurrent getter changed worker native FP/raw CSR/errno");}
    require(pointers[0]&&pointers[0]==pointers[1]&&pointers[0]==&candidate->verts(label),
        "pool simultaneous first getters did not return the same stable real object");
    require(record(*pointers[0],0)==record(Vert{},0),"pool concurrently published entity was partially initialized");
    const auto after=stock(*candidate);
    require(after[0].constructedSlots==before[0].constructedSlots+1
        &&after[0].publishedChunks==before[0].publishedChunks+1,"pool last short chunk published more than once or wrong slot count");
    auto expected=before;expected[0]=after[0];sameStock(expected,after);
    passiveCache(*candidate,scratch);++concurrentPairs;
}
void cacheTransactions(const QuadtreeMesh& reference,const DictionaryMesh& input,const fs::path& scratch)
{
    const auto valid=zeroLsCache(reference,scratch/"valid");auto receiver=reference.clone();
    const auto original=passiveCache(*receiver,scratch/"receiver_entry");const auto originalStock=stock(*receiver);
    const auto reject=[&](std::vector<char> bad,const char* name){const auto folder=scratch/name;zeroLsWrite(folder,bad);
        bool rejected=false;try{receiver->readCache(folder.generic_string()+"/");}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"pool invalid native cache was accepted");sameStock(originalStock,stock(*receiver));
        require(bytes(folder/"mesh_cache.bin")==bad&&passiveCache(*receiver,scratch/(std::string(name)+"_state"))==original,
            "pool rejected cache partially committed storage/metadata or changed input file");++cacheRejections;};
    auto bad=valid;const int negative=-1;std::memcpy(bad.data()+8,&negative,sizeof negative);reject(bad,"negative_count");
    bad=valid;bad.pop_back();reject(bad,"truncated_tail");
    receiver->readCache((scratch/"valid").generic_string()+"/");dense(*receiver);
    require(passiveCache(*receiver,scratch/"restored")==valid,"pool dense-adopt QMC2 changed complete native payload");
    sameMesh(*receiver,reference,scratch/"full_adopt_oracle");
    auto twoInput=input;twoInput.levelSetCount=2;twoInput.minLevel_neg={0,0};twoInput.minLevel_pos={0,0};twoInput.maxLevel_local={3,3};
    QuadtreeMesh positiveReceiver(twoInput);positiveReceiver.readCache((scratch/"valid").generic_string()+"/");
    dense(positiveReceiver);require(positiveReceiver.levelSetCount()==0,"pool cache adoption retained receiver LS width");
    sameMesh(positiveReceiver,reference,scratch/"positive_receiver_adopt");
}
void initialSource(const QuadtreeMesh& mesh,const BackgroundMesh& original,int capacity)
{
    require(mesh.levelSetCount()==0&&original.levelSetCount()==0,"pool initial source changed its explicit zero-LS request");
    const auto before=stock(mesh);
    const std::array<int,4> count{{mesh.vertIndexNum(),mesh.hEdgeIndexNum(),mesh.vEdgeIndexNum(),mesh.cellIndexNum()}};
    require(count==std::array<int,4>{{original.vertIndexNum(),original.hEdgeIndexNum(),original.vEdgeIndexNum(),original.cellIndexNum()}},
        "pool initial indexed prefix differs from original BackgroundMesh");
    for(std::size_t k=0;k<before.size();++k){const auto chunks=(static_cast<std::size_t>(count[k])+255)/256;
        require(!before[k].denseMode&&before[k].logicalCapacity==static_cast<std::size_t>(capacity)
            &&before[k].publishedChunks==chunks&&before[k].constructedSlots==std::min(chunks*256,before[k].logicalCapacity)
            &&before[k].constructedSlots<before[k].logicalCapacity,"pool initial construction did not stop at chunks covering the original active prefix");}
    const auto prefix=[&](std::size_t k,const int* actualIndex,const int* originalIndex,auto value,auto oldValue){
        for(std::size_t label=0;label<before[k].constructedSlots;++label){const auto& t=value(static_cast<int>(label));using Entity=std::decay_t<decltype(t)>;
            if constexpr(std::is_same_v<Entity,Vert>)require(t.onCutEdge.empty()&&t.onCutCell.empty()&&t.nearestNodeCoord.empty()&&t.nearestNodeLabel.empty(),
                "pool initial constructed vertex allocated zero-LS vectors");
            else require(t.isCut.empty(),"pool initial constructed edge/cell allocated zero-LS vectors");
            if(label<static_cast<std::size_t>(count[k]))require(actualIndex[label]==originalIndex[label]
                &&record(t,0)==record(oldValue(static_cast<int>(label)),0),"pool initial live index/value differs from independent BackgroundMesh");
            else require(record(t,0)==record(Entity{},0),"pool initial chunk free slot is not the original default entity");}
    };
    prefix(0,mesh.vertIndexToLabel,original.vertIndexToLabel,[&](int i)->const auto&{return mesh.verts(i);},[&](int i)->const auto&{return original.verts(i);});
    prefix(1,mesh.hEdgeIndexToLabel,original.hEdgeIndexToLabel,[&](int i)->const auto&{return mesh.hEdges(i);},[&](int i)->const auto&{return original.hEdges(i);});
    prefix(2,mesh.vEdgeIndexToLabel,original.vEdgeIndexToLabel,[&](int i)->const auto&{return mesh.vEdges(i);},[&](int i)->const auto&{return original.vEdges(i);});
    prefix(3,mesh.cellIndexToLabel,original.cellIndexToLabel,[&](int i)->const auto&{return mesh.cells(i);},[&](int i)->const auto&{return original.cells(i);});
    sameStock(before,stock(mesh));++initialSources;
}
void scene(const char* name,Vector2D lower,Vector2D upper,const fs::path& scratch)
{
    DictionaryMesh input;input.minCoord=lower;input.maxCoord=upper;input.Nx=input.Ny=4;input.levelSetCount=0;
    input.maxCellNum=4097;input.maxInitialIterations=0;input.maxLevel_local={3};
    QuadtreeMesh source(input),literal(input);const BackgroundMesh original(input);initialSource(source,original,input.maxCellNum);
    const auto initial=passiveCache(source,scratch/"initial_source");
    literal.readCache((scratch/"initial_source").generic_string()+"/");dense(literal);
    require(passiveCache(literal,scratch/"initial_dense_reference")==initial,"pool initial lazy writer/dense-adopt complete QMC2 bytes differ");
    sameMesh(source,literal,scratch/"initial_full_oracle");++initialWriterPairs;
    require(settle(source,literal,{1,lower,upper},scratch/"dense_refine")>0,"pool source lacks actual mixed refinement");
    require(settle(source,literal,{2,lower,upper},scratch/"dense_move")>0,"pool source lacks actual retirement/reallocation");
    std::size_t localParents=0;
    for(int i=0;i<source.cellLabelCapacity();++i)if(source.cells(i).label()!=-1&&source.cells(i).son[0]>=0){
        require(std::find(source.cellIndexToLabel,source.cellIndexToLabel+source.cellIndexNum(),i)==source.cellIndexToLabel+source.cellIndexNum(),
            "pool source live ancestor is incorrectly an indexed leaf");++localParents;}
    require(localParents>0,"pool independent source has no live nonindexed parent records");parentRecords+=localParents;
    dirtyFree(source);dirtyFree(literal);const NativeMeshIdentity sourceIdentity(source);
    const auto sourceBytes=zeroLsCache(source,scratch/"source_dirty");resetFreeReference(literal);
    const auto expected=zeroLsCache(literal,scratch/"independent_dense_clone");
    require(sourceBytes!=expected,"pool original free-payload reset fixture is vacuous");
    const Floating entry;auto candidate=source.clone();const Floating exit;sameFloating(entry,exit);sparse(*candidate);
    require(passiveCache(*candidate,scratch/"lazy_clone")==expected,"pool lazy clone differs from ORIGINAL dense clone free/live rule");++writerPairs;
    cloneEnvironment(source,expected,scratch/"clone_environment");concurrentFirstGetter(source,scratch/"concurrent");
    // Observe sparse storage BEFORE sameMesh intentionally visits all slots.
    const auto before=stock(*candidate);const auto first=before[1].constructedSlots+before[2].constructedSlots+before[3].constructedSlots;
    firstMaintenance(*candidate,literal,scratch/"first_sparse_noop");
    // firstMaintenance() has intentionally materialized the full capacity for
    // unchanged all-record/QMC2 oracle. Fresh clone+writer above is the sparse observation.
    auto highCandidate=literal.clone();sparse(*highCandidate);
    mutableHighLabels(*highCandidate,literal,scratch/"high_labels");sameMesh(*highCandidate,literal,scratch/"high_records");
    // Retired free payload is reset once more before real subsequent allocation.
    auto evolving=highCandidate->clone();resetFreeReference(literal);
    require(passiveCache(*evolving,scratch/"repeat_clone")==zeroLsCache(literal,scratch/"repeat_dense"),
        "pool second original free reset disagrees with independent dense reference");
    require(settle(*evolving,literal,{0,lower,upper},scratch/"coarsen")>0,"pool real clone coarsening did not retire slots");
    require(settle(*evolving,literal,{3,lower,upper},scratch/"rerefine")>0,"pool real clone rerefinement did not allocate/recycle slots");
    cacheTransactions(literal,input,scratch/"cache");
    sourceIdentity.validate(source);require(zeroLsCache(source,scratch/"source_after")==sourceBytes,
        "pool clone/adapt/cache/getter operations mutated original source including free payload");
    std::cout<<"ENTITY_POOL_STORAGE_SCENE,"<<name<<",lazy_initial_slots,"<<before[0].constructedSlots<<'/'<<before[1].constructedSlots<<'/'
        <<before[2].constructedSlots<<'/'<<before[3].constructedSlots<<",logical_capacity,"<<before[0].logicalCapacity<<",entry_transfer_slots,"<<first<<'\n';++scenes;
}
void positive(const fs::path& scratch)
{
    for(const int count:{1,2}){DictionaryMesh input;input.minCoord={0,0};input.maxCoord={4,4};input.Nx=input.Ny=4;
        input.levelSetCount=count;input.maxCellNum=513;input.maxInitialIterations=0;input.minLevel_neg.assign(count,0);
        input.minLevel_pos.assign(count,0);input.maxLevel_local.assign(count,3);
        QuadtreeMesh source(input);dense(source);auto candidate=source.clone();dense(*candidate);
        checkLsWidth(source,count);checkLsWidth(*candidate,count);
        sameMesh(*candidate,source,scratch/std::to_string(count));++positiveCases;}
}
void run(const fs::path& scratch)
{
    scene("XY",{-2,-1},{2,3},scratch/"xy");scene("shifted_coordinates",{.125,-.25},{4.125,3.75},scratch/"shifted");positive(scratch/"positive");
    require(scenes==2&&parentRecords>0&&writerPairs==4&&firstClears==2&&highLabels==6&&concurrentPairs==2
        &&cacheRejections==4&&positiveCases==2&&environmentPairs>=2&&initialSources==2&&initialWriterPairs==2,"pool independent finite storage coverage incomplete");
    std::cout<<"ENTITY_POOL_STORAGE_PASS,scenes,"<<scenes<<",independent_dense_reference_checks,"<<denseReferences
        <<",nonindexed_live_parents,"<<parentRecords<<",passive_full_QMC2_pairs,"<<writerPairs<<",first_sparse_maintenance,"<<firstClears
        <<",high_label_mutable_records,"<<highLabels<<",concurrent_first_getters,"<<concurrentPairs<<",cache_transaction_rejections,"<<cacheRejections
        <<",positive_dense_LS_cases,"<<positiveCases<<",full_FP_errno_clone_pairs,"<<environmentPairs
        <<",actual_initial_lazy_sources,"<<initialSources<<",initial_passive_writer_dense_adopt_pairs,"<<initialWriterPairs
        <<",allocation_fault_injection_NA,private_unused_tail_NA,no_speed_claim\n";
}
}

int main(int argc,char** argv)
{
    if(argc==2&&std::string(argv[1])=="--entity-pool-storage-only"){
        const auto scratch=fs::current_path()/("mesh_entity_pool_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        int result=0;try{EntityPoolStorageContract::run(scratch);}catch(const std::exception& e){std::cerr<<"ENTITY_POOL_STORAGE_FAIL "<<e.what()<<'\n';result=1;}
        std::error_code ignored;fs::remove_all(scratch,ignored);return result;
    }
    if(argc==2&&std::string(argv[1])=="--native-tree-workspace-only"){
        const auto scratch=fs::current_path()/("mesh_native_tree_workspace_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        int result=0;try{NativeTreeWorkspaceContract::run(scratch);}catch(const std::exception& e){std::cerr<<"NATIVE_TREE_WORKSPACE_FAIL "<<e.what()<<'\n';result=1;}
        std::error_code ignored;fs::remove_all(scratch,ignored);return result;
    }
    if(argc==2&&std::string(argv[1])=="--zero-ls-only"){
        const auto scratch=fs::current_path()/("mesh_zero_ls_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        int result=0;try{zeroLsCases(scratch);}catch(const std::exception& e){std::cerr<<"ZERO_LS_FAIL "<<e.what()<<'\n';result=1;}
        std::error_code ignored;fs::remove_all(scratch,ignored);return result;
    }
    if(argc==2&&std::string(argv[1])=="--refinement-prerequisites-only"){
        try{refinementPrerequisiteCases();return 0;}catch(const std::exception& e){std::cerr<<"REFINEMENT_PREREQUISITES_FAIL "<<e.what()<<'\n';return 1;}
    }
    if(argc==2&&std::string(argv[1])=="--thin-mesh-only"){
        try{thinMeshTransferCases();return 0;}catch(const std::exception& e){std::cerr<<"THIN_MESH_FAIL "<<e.what()<<'\n';return 1;}
    }
    if(argc==2&&std::string(argv[1])=="--native-target-only"){
        try{nativeTargetTransfer();return 0;}catch(const std::exception& e){std::cerr<<"NATIVE_TARGET_ONLY_FAIL "<<e.what()<<'\n';return 1;}
    }
    if(argc!=1){std::cerr<<"usage: mesh_adaptation_contract [--native-target-only|--thin-mesh-only]\n";return 2;}
    const auto scratch = fs::current_path() / ("mesh_adaptation_" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    int rc = 0;
    try {
        verify({"unit_square_level1",{0,0},{1,1},4,4,1,{.01,.01},{.99,.49},40}, scratch/"unit");
        verify({"translated_rectangle_level2",{-1.5,.5},{2.5,3.5},5,3,2,{-.7,.5},{.9,2.},0}, scratch/"translated");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        rc = 1;
    }
    std::error_code ignored;
    fs::remove_all(scratch,ignored);
    return rc;
}
