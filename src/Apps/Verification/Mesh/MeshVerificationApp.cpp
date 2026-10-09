#include "Verification/Mesh/MeshVerificationApp.h"
#include "Core/IO/StrictDictionary.h"
#include "Core/IO/TecplotMeshOutput.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>

namespace {
namespace fs = std::filesystem;
void require(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}
Vector2D pair(const DictionaryInput::Node& node)
{
    auto tokens = node.values;
    tokens.erase(std::remove(tokens.begin(), tokens.end(), ","), tokens.end());
    if (tokens.size() != 4 || tokens.front() != "(" || tokens.back() != ")")
        node.fail("two finite code-unit values (x, y)");
    auto value = node;
    value.values = {tokens[1]}; const double x = value.number();
    value.values = {tokens[2]}; return {x, value.number()};
}
DictionaryMesh readInput()
{
    for (const auto& entry : fs::directory_iterator("system")) {
        const auto name = entry.path().filename().string();
        require(entry.is_regular_file() && (name == "run.dict" || name == "geometry.dict" || name == "mesh.dict"),
            "MeshVerification consumes only run.dict, geometry.dict and mesh.dict");
    }
    DictionaryInput::Document run("system/run.dict"), geometry("system/geometry.dict"), mesh("system/mesh.dict");
    run.root().allow({"Run"});
    run.root().block("Run").allow({"app"});
    run.root().block("Run").get("app").choice({"MeshVerification"});
    geometry.root().allow({"Domain"});
    const auto& domain = geometry.root().block("Domain"); domain.allow({"min", "max"});
    const auto lo = pair(domain.get("min")), hi = pair(domain.get("max"));
    require(hi.x > lo.x && hi.y > lo.y && (hi.x-lo.x)*(hi.y-lo.y) > 0. &&
        std::isfinite((hi.x-lo.x)*(hi.y-lo.y)), "invalid finite mesh domain");
    const auto& root = mesh.root(); root.allow({"adaptiveMesh", "BackgroundMesh", "QuadtreeMesh", "GeometryRefinement"});
    root.get("adaptiveMesh").choice({"true"});
    root.block("BackgroundMesh").allow({"divide"});
    const auto divide = pair(root.block("BackgroundMesh").get("divide"));
    require(divide.x >= 2 && divide.y >= 2 && divide.x <= 64 && divide.y <= 64 &&
        divide.x == std::floor(divide.x) && divide.y == std::floor(divide.y), "MeshVerification requires integer Nx/Ny in [2,64]");
    const auto& quad = root.block("QuadtreeMesh"); quad.allow({"maxCellNum", "maxInitialIterations"});
    const int capacity = quad.get("maxCellNum").integer();
    require(quad.get("maxInitialIterations").integer() == 0, "MeshVerification schedules adaptation explicitly; maxInitialIterations must be 0");
    const auto& refinement = root.block("GeometryRefinement"); refinement.allow({"enabled", "refineRegion"});
    refinement.get("enabled").choice({"true"});
    refinement.block("refineRegion").allow({"entry"});
    const auto& box = refinement.block("refineRegion").block("entry");
    box.allow({"name", "shape", "xmin", "xmax", "ymin", "ymax", "maxLevel"});
    const auto name = box.get("name").word(); box.get("shape").choice({"box"});
    const Vector2D boxLo{box.get("xmin").number(), box.get("ymin").number()};
    const Vector2D boxHi{box.get("xmax").number(), box.get("ymax").number()};
    const int level = box.get("maxLevel").integer();
    require(boxLo.x >= lo.x && boxHi.x <= hi.x && boxLo.y >= lo.y && boxHi.y <= hi.y &&
        boxHi.x > boxLo.x && boxHi.y > boxLo.y, "refinement box must have positive area inside the domain");
    require(level >= 1 && level <= 3, "MeshVerification supports maxLevel in [1,3]");
    const int nx = static_cast<int>(divide.x) << level, ny = static_cast<int>(divide.y) << level;
    require(nx*ny <= 4096 && capacity >= 2*(nx+1)*(ny+1) && capacity <= 65536,
        "MeshVerification requires <=4096 finest uniform leaves and 2*(Nx*2^level+1)*(Ny*2^level+1) <= maxCellNum <= 65536");
    // Keep the established reader as the source of the actual mesh/rules.
    // This entry's capability filter rejects every unused physical setting.
    DictionaryMesh result; result.readDict();
    require(result.Nx == divide.x && result.Ny == divide.y && result.maxCellNum == capacity && result.maxInitialIterations == 0 &&
        result.minCoord.x == lo.x && result.minCoord.y == lo.y && result.maxCoord.x == hi.x && result.maxCoord.y == hi.y &&
        result.geometryRefinement.enabled && result.geometryRefinement.regionRules.size() == 1 && result.maxLevel_local.front() == level,
        "mesh input must use the existing reader's multiline block/record syntax");
    const auto& rule = result.geometryRefinement.regionRules.front();
    require(rule.name == name && rule.shape == DictionaryMesh::GeometryRefinementConfig::RegionShape::kBox && rule.maxLevel == level &&
        rule.min.x == boxLo.x && rule.min.y == boxLo.y && rule.max.x == boxHi.x && rule.max.y == boxHi.y,
        "region values were not consumed by the existing mesh reader");
    return result;
}
double reference(const QuadtreeMesh& mesh, Vector2D point)
{
    const double x = (point.x-mesh.minCoord.x)/(mesh.maxCoord.x-mesh.minCoord.x);
    const double y = (point.y-mesh.minCoord.y)/(mesh.maxCoord.y-mesh.minCoord.y);
    return 1. + 2.*x - 3.*y + 4.*x*y;
}
struct Probe {
    std::vector<VertField<double>> nodal;
    CellField<double> constant;
    explicit Probe(QuadtreeMesh& mesh) : nodal{VertField<double>(mesh,0.)}, constant(mesh,2.) {
        for (int i=0;i<mesh.vertIndexNum();++i) nodal[0].value[i] = reference(mesh,mesh.verts(mesh.vertIndexToLabel[i]).coord());
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
    throw std::runtime_error("mesh verification targets did not settle within the refinement-level bound");
}
struct Metrics { double area=0., nodalError=0., cellError=0., closureError=0.; int coarseFine=0, level=0; };
Metrics measure(const QuadtreeMesh& mesh, const Probe& probe)
{
    require(probe.nodal[0].size() == mesh.vertIndexNum() && probe.constant.size() == mesh.cellIndexNum(), "stale transferred probe fields");
    Metrics result;
    for (int i=0;i<mesh.vertIndexNum();++i) {
        const auto& vertex = mesh.verts(mesh.vertIndexToLabel[i]);
        require(vertex.index() == i && probe.nodal[0].label(i) == vertex.label(), "invalid vertex/field identity");
        const double error = std::abs(probe.nodal[0].value[i]-reference(mesh,vertex.coord()));
        require(std::isfinite(error), "nonfinite transferred nodal probe"); result.nodalError = std::max(result.nodalError,error);
    }
    FaceSegmentGraph graph; graph.build(mesh);
    for (int i=0;i<mesh.cellIndexNum();++i) {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
        require(cell.index() == i && cell.son[0] == -1 && probe.constant.label(i) == cell.label(), "invalid leaf/field identity");
        result.area += cell.area(); result.level = std::max(result.level,cell.level());
        const double error = std::abs(probe.constant.value[i]-2.);
        require(std::isfinite(error), "nonfinite transferred cell probe"); result.cellError = std::max(result.cellError,error);
        Vector2D closure{0.,0.}; const auto range = graph.segmentRangeForCellIndex(i);
        for (auto j=range.first;j<range.second;++j) {
            const auto& segment = graph.segments()[j];
            require(segment.K == cell.label() && segment.Af > 0. && segment.df > 0., "invalid face segment geometry");
            closure = closure + segment.n*segment.Af;
            if (segment.L >= 0) {
                const int difference = std::abs(cell.level()-mesh.cells(segment.L).level());
                require(difference <= 1, "mesh violates 2:1 face-neighbor balance"); result.coarseFine += difference;
            }
        }
        result.closureError = std::max(result.closureError,std::hypot(closure.x,closure.y)/(cell.size().x+cell.size().y));
    }
    const double expectedArea = (mesh.maxCoord.x-mesh.minCoord.x)*(mesh.maxCoord.y-mesh.minCoord.y);
    require(std::abs(result.area/expectedArea-1.) <= 1e-12 && result.nodalError <= 1e-11 &&
        result.cellError <= 1e-11 && result.closureError <= 1e-12, "mesh area, field transfer or face-closure reference failed");
    return result;
}
void publish(const QuadtreeMesh& mesh, const Probe& probe, int stage, int sweeps, std::ofstream& history)
{
    const auto metrics = measure(mesh,probe);
    auto zone = TecplotMeshOutput::nodalZone(mesh,{"bilinear_probe"},probe.nodal);
    zone.title = "Core mesh verification stage " + std::to_string(stage);
    zone.auxiliary = {{"coordinates","code"},{"stage",std::to_string(stage)},{"physical_time","none"}};
    zone.fields.push_back({"constant_probe",TecplotWriter::Location::CellCentered,{}});
    zone.fields.push_back({"level",TecplotWriter::Location::CellCentered,{}});
    for (int i=0;i<mesh.cellIndexNum();++i) {
        zone.fields[1].values.push_back(probe.constant.value[i]);
        zone.fields[2].values.push_back(mesh.cells(mesh.cellIndexToLabel[i]).level());
    }
    TecplotWriter::write(fs::path("output/fields")/("mesh_stage"+std::to_string(stage)+".dat"),zone);
    history << stage << ',' << sweeps << ',' << mesh.cellIndexNum() << ',' << mesh.vertIndexNum() << ',' << metrics.level << ','
        << metrics.coarseFine << ',' << metrics.area << ',' << metrics.nodalError << ',' << metrics.cellError << ',' << metrics.closureError << '\n';
    history.flush(); require(bool(history), "failed writing mesh verification history");
}
std::vector<char> bytes(const fs::path& path)
{
    std::ifstream file(path,std::ios::binary); require(bool(file), "missing mesh cache");
    return {std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>()};
}
}

int MeshVerificationApp::run()
{
    const char* stage = "input";
    try {
        const auto input = readInput();
        if (checkOnly_) { std::cout << "[MeshVerification] Read-only input/capability check passed.\n"; return 0; }
        require(!fs::exists("output") && !fs::exists("cache"), "MeshVerification requires a new run directory; existing output/cache is preserved");
        stage = "initial mesh";
        QuadtreeMesh mesh(input); Probe probe(mesh); QuadtreeMeshRefiner refiner(input);
        require(refiner.attachGeometry(nullptr,nullptr), "region-only refinement did not initialize");
        fs::create_directories("output/fields"); fs::create_directories("output/series"); fs::create_directories("output/meta");
        std::ofstream history("output/series/mesh_verification.csv"); history << std::setprecision(17);
        history << "stage,sweeps,cells,vertices,max_level,oriented_coarse_fine_segments,area,nodal_error,cell_error,face_closure_error\n";
        publish(mesh,probe,0,0,history);
        stage = "region refinement"; const int firstSweeps = adapt(mesh,probe,refiner,false);
        require(mesh.cellIndexNum() > input.Nx*input.Ny, "verification did not exercise actual refinement");
        publish(mesh,probe,1,firstSweeps,history);
        stage = "cache roundtrip";
        fs::create_directories("cache/refined"); fs::create_directories("cache/roundtrip");
        mesh.writeCache("cache/refined/");
        QuadtreeMesh restored(input); restored.readCache("cache/refined/");
        restored.writeCache("cache/roundtrip/");
        require(bytes("cache/refined/mesh_cache.bin") == bytes("cache/roundtrip/mesh_cache.bin"), "QMC2 roundtrip changed bytes");
        // QMC2 stores mesh only. Copy live probe fields with the verified same
        // native indexing; do not claim a physical-field checkpoint restore.
        Probe continued = probe; measure(restored,continued);
        stage = "coarsening restored mesh"; const int coarseSweeps = adapt(restored,continued,refiner,true);
        require(restored.cellIndexNum() == input.Nx*input.Ny, "restored mesh did not return to the background leaves");
        publish(restored,continued,2,coarseSweeps,history);
        stage = "refinement after coarsening"; const int secondSweeps = adapt(restored,continued,refiner,false);
        publish(restored,continued,3,secondSweeps,history);
        const auto original = TecplotMeshOutput::nodalZone(mesh,{"bilinear_probe"},probe.nodal);
        const auto repeated = TecplotMeshOutput::nodalZone(restored,{"bilinear_probe"},continued.nodal);
        require(original.connectivity == repeated.connectivity && original.fields[0].values == repeated.fields[0].values &&
            original.points.size() == repeated.points.size(), "refine/coarsen/refine changed nodal topology or values");
        for (std::size_t i=0;i<original.points.size();++i) require(original.points[i].x == repeated.points[i].x &&
            original.points[i].y == repeated.points[i].y, "repeated refinement changed coordinates");
        std::ofstream metadata("output/meta/mesh_verification.txt");
        const auto& box = input.geometryRefinement.regionRules.front();
        metadata << std::setprecision(17) << "app=MeshVerification\ncoordinates=code\nphysical_time=none\naccepted_stages=4\ncache=QMC2_mesh_only\ncache_roundtrip_exact=1\n"
            << "domain_min=" << input.minCoord.x << ',' << input.minCoord.y << "\ndomain_max=" << input.maxCoord.x << ',' << input.maxCoord.y
            << "\ndivide=" << input.Nx << ',' << input.Ny << "\nmaxCellNum=" << input.maxCellNum << "\nregion_name=" << box.name
            << "\nregion_min=" << box.min.x << ',' << box.min.y << "\nregion_max=" << box.max.x << ',' << box.max.y << "\nmaxLevel=" << box.maxLevel << '\n'
            << "nodal_reference=1+2*x-3*y+4*x*y (domain-normalized coordinates)\ncell_reference=2\n"
            << "field_tolerance=1e-11\nrelative_area_and_face_closure_tolerance=1e-12\n"
            << "fields=output/fields/mesh_stage[0-3].dat\nhistory=output/series/mesh_verification.csv\n";
        metadata.flush(); require(bool(metadata), "failed writing mesh verification metadata");
        std::cout << "[MeshVerification] Four stages passed: initial, refined, restored/coarsened, refined again; QMC2 roundtrip exact.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[MeshVerification] failed at " << stage << ": " << error.what() << '\n';
        return std::string(stage) == "input" ? 2 : 3;
    }
}
