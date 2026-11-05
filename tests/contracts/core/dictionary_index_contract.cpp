#include "Core/IO/Dictionary.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/DictionaryBlockText.h"
#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include <functional>
#include <map>
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Units/UnitsContext.h"
#include <cmath>
#include <limits>
#include <set>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

struct Reader : Dictionary {
    using Dictionary::parseIndex;
    using Dictionary::parseIndex_brace;
};
template<class Error, class F>
void readerRejects(F action, const std::filesystem::path& source, const char* reason)
{
    bool caught=false;
    try { action(); }
    catch (const Error& error) {
        const std::string message=error.what();
        caught=message.find(source.string())!=std::string::npos && message.find(reason)!=std::string::npos;
    }
    if (!caught) throw std::runtime_error("dictionary failure lost its exception type, source or reason");
}
void readerFailures()
{
    namespace fs=std::filesystem;
    const auto cwd=fs::current_path();
    const auto folder=cwd/("dictionary_failures_"+
        std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    fs::create_directories(folder);
    const auto geometry=folder/"shape.dict", mesh=folder/"grid.dict";
    const auto absent=folder/"absent.dict";
    const auto write=[](const fs::path& path,const std::string& text) {
        std::ofstream stream(path); stream<<text;
        if (!stream) throw std::runtime_error("cannot write dictionary fixture");
    };
    const std::string validGeometry="Domain\n{\n min (0.2, -0.5)\n max (1.2, 0.5)\n}\n";
    const std::string validMesh="adaptiveMesh false\nBackgroundMesh\n{\n divide (8, 10)\n}\n";
    const auto readMesh=[&] { DictionaryMesh value; value.readDict(mesh.string(),geometry.string()); };
    readerRejects<std::runtime_error>([&] { DictionaryGeometry value; value.readDict(absent.string()); },
        absent,"cannot open geometry dictionary");
    readerRejects<std::runtime_error>([&] { DictionaryMesh value; value.readDict(absent.string(),geometry.string()); },
        absent,"cannot open mesh dictionary");
    write(mesh,validMesh);
    readerRejects<std::runtime_error>(readMesh,geometry,"cannot open geometry dictionary");
    write(geometry,validGeometry);
    write(mesh,"BackgroundMesh\n{\n divide (0, 10)\n}\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"invalid grid division");
    write(mesh,"adaptiveMesh false\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"missing divide");
    write(mesh,validMesh);
    write(geometry,"Domain\n{\n min (0.2, -0.5)\n}\n");
    readerRejects<std::invalid_argument>(readMesh,geometry,"missing Domain/min/max");
    write(geometry,"Domain\n{\n min (0.2, -0.5)\n max (0.1, 0.5)\n}\n");
    readerRejects<std::invalid_argument>(readMesh,geometry,"invalid domain bounds");
    write(geometry,validGeometry);
    const std::string region="GeometryRefinement\n{\n refineRegion\n{\n entry\n{\n";
    write(mesh,validMesh+region+" name focus\n shape circle\n}\n}\n}\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"entry 'focus' uses unsupported shape 'circle'");
    write(mesh,validMesh+region+" shape box\n}\n}\n}\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"<unnamed> is missing required field 'name'");
    write(mesh,validMesh+region+" name focus\n");
    readerRejects<std::invalid_argument>(readMesh,mesh,"unterminated GeometryRefinement.refineRegion entry");
    // A new record remains usable in the same process after every failure.
    write(mesh,validMesh);
    DictionaryGeometry parsedGeometry; parsedGeometry.readDict(geometry.string());
    DictionaryMesh parsedMesh; parsedMesh.readDict(mesh.string(),geometry.string());
    if (!parsedGeometry.hasDomain || parsedGeometry.domainMin.x!=.2 || parsedMesh.Nx!=8 || parsedMesh.Ny!=10
        || parsedMesh.minCoord.y!=-.5 || parsedMesh.maxCoord.x!=1.2 || fs::current_path()!=cwd
        || fs::exists(folder/"output") || fs::exists(folder/"cache"))
        throw std::runtime_error("dictionary failure changed process state or later valid input");
}
void extensions()
{
    struct Restore {
        std::filesystem::path original=std::filesystem::current_path();
        ~Restore() { std::filesystem::current_path(original); }
    } restore;
    const auto folder=restore.original/("dictionary_extensions_"+
        std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(folder/"system"); std::filesystem::current_path(folder);
    {
        std::ofstream out("system/geometry.dict");
        out<<"ExtraPhysics\n{\n enabled off\n name \"literal {\"\n Domain\n{\n min (99, 99)\n max (100, 100)\n enabled true\n}\n}\n"
              "Domain\n{\n min (0, 0)\n max (1, 1)\n}\n";
    }
    DictionaryGeometry geometry; geometry.readDict();
    if (!geometry.hasDomain || geometry.domainMin.x!=0. || geometry.domainMax.x!=1.
        || !geometry.extensionBlock("ExtraPhysics") || geometry.extensionEnabled("ExtraPhysics")
        || geometry.extensionEnabled("Absent"))
        throw std::runtime_error("opaque extension nesting/capability/source isolation");
}
void blockScope()
{
    DictionaryBlockText reader;
    using Scope=DictionaryBlockText::Scope;
    const auto check=[](bool condition,const char* message) {
        if(!condition)throw std::runtime_error(message);
    };
    std::string body="stale";
    check(reader.extractBlockText("\xEF\xBB\xBFReactions { bom 1 }","Reactions",body,Scope::TopLevel)
        && body==" bom 1 ","UTF-8 BOM must not hide the first block");
    const std::string text=R"(Output {
 Diagnostics { Reactions { enabled true } }
 label "literal } Reactions {"
}
reactionNames (Reactions { ignored })
name Reactions
Other { value 9 }
// Reactions { ignored }
Reactions # ignored { }
{
 physical true
 nested { value "escaped \" } {" }
}
)";
    check(reader.extractBlockText(text,"Reactions",body,Scope::TopLevel) &&
        body.find("physical true")!=std::string::npos && body.find("escaped")!=std::string::npos,
        "root block must skip diagnostic, list, scalar, quoted and commented names/braces");
    check(reader.extractBlockText(text,"Reactions",body) && body==" enabled true ",
        "nested solver block extraction remains available explicitly or by default");
    for(const auto& input:{std::string("Reactions false\nOther { x 1 }"),
        std::string("Output { Reactions { enabled true } }"),std::string("ReactionsExtra { x 1 }"),
        std::string("Reactions { missing"),std::string("\"Reactions\" { x 1 }")})
        check(!reader.extractBlockText(input,"Reactions",body,Scope::TopLevel) && body.empty(),
              "false block match or stale body on failed extraction");
    check(reader.extractFirstBlockText("gas { second 2 }\nGas { first 1 }\nGas { later 3 }",
        {"Gas","gas"},body,Scope::TopLevel) && body==" first 1 ","caller alias priority or first occurrence changed");
    check(!reader.extractFirstBlockText("Gas {}",{},body) && body.empty(),"empty aliases retain stale body");
    check(reader.extractBlockText("A {}; Reactions { x 1 }","Reactions",body,Scope::TopLevel)
        && body==" x 1 ","adjacent top-level blocks");
}

namespace
{
void csgCheck(bool condition,const char* reason)
{ if(!condition) throw std::runtime_error(reason); }
template<class F> void csgReject(F action,const char* reason)
{
    bool rejected=false;
    try { action(); } catch(const std::exception&) { rejected=true; }
    csgCheck(rejected,reason);
}
DictionaryMesh csgMeshInput(Vector2D upper,int nx,int ny,int level=1)
{
    DictionaryMesh input;
    input.minCoord={0,0};input.maxCoord=upper;input.Nx=nx;input.Ny=ny;
    input.maxCellNum=std::max(256,8*nx*ny);input.maxInitialIterations=0;input.levelSetCount=0;
    input.maxLevel_local={level};return input;
}
void csgRectangle(DictionaryGeometry& geometry,Vector2D lo,Vector2D hi)
{
    const int start=static_cast<int>(geometry.lines.size());
    const Vector2D points[]={{lo.x,lo.y},{hi.x,lo.y},{hi.x,hi.y},{lo.x,hi.y}};
    std::vector<SignInt> path;
    for(int i=0;i<4;++i)
    {
        Line line;line.points[0]=points[i];line.points[1]=points[(i+1)%4];
        line.setType(LineType::straight);geometry.lines.push_back(line);path.push_back({start+i,false});
    }
    geometry.levelSetName.push_back("rectangle"+std::to_string(geometry.levelSetName.size()));
    geometry.innerPoints.push_back({0.5*(lo.x+hi.x),0.5*(lo.y+hi.y)});
    geometry.moving.push_back(false);geometry.linesOnInterface.push_back(path);
}
void csgExpression(DictionaryGeometry& geometry,const std::vector<LevelSetOperation>& operations,
    const std::vector<SignInt>& a,const std::vector<SignInt>& b)
{
    geometry.phaseName.push_back("phase"+std::to_string(geometry.phaseName.size()));
    geometry.levelSetOperation.push_back(operations);geometry.levelSetIndexA.push_back(a);
    geometry.levelSetIndexB.push_back(b);geometry.boundaryLines.emplace_back();
}
void csgBooleanValues()
{
    DictionaryGeometry geometry;csgRectangle(geometry,{.5,.5},{2.5,2.5});
    csgRectangle(geometry,{1.5,1.5},{3.5,3.5});
    csgExpression(geometry,{LevelSetOperation::min},{{0,false}},{{1,false}});
    csgExpression(geometry,{LevelSetOperation::max},{{0,false}},{{1,false}});
    csgExpression(geometry,{LevelSetOperation::max},{{0,false}},{{1,true}});
    csgExpression(geometry,{LevelSetOperation::non},{{0,true}},{{0,true}});
    csgExpression(geometry,{LevelSetOperation::min,LevelSetOperation::max,LevelSetOperation::max},
        {{0,false},{0,false},{2,false}},{{1,false},{1,true},{3,true}});
    auto input=csgMeshInput({4,4},8,8);QuadtreeMesh mesh(input);
    const FittedDomainGeometry::InitialPhaseSampler sampler(mesh,geometry);
    // The constructor's ImmersedGeometry seed and its dictionary copy have
    // already died. Mutating the source must not mutate the owned sampler.
    geometry.lines[0].points[0]={99,99};
    csgCheck(sampler.levelSetValue(0,{.5,1.0})==0.0,"initial CSG exact zero was displaced");
    for(double x:{.25,.75,1.75,2.75,3.75}) for(double y:{.25,.75,1.75,2.75,3.75})
    {
        const Vector2D point{x,y};const double a=sampler.levelSetValue(0,point),b=sampler.levelSetValue(1,point);
        const bool inA=x>.5&&x<2.5&&y>.5&&y<2.5;
        const bool inB=x>1.5&&x<3.5&&y>1.5&&y<3.5;
        csgCheck((a<0)==inA&&(b<0)==inB,"owned primitive sign disagrees with rectangle oracle");
        const bool expected[]={inA||inB,inA&&inB,inA&&!inB,!inA,(inA||inB)&&!(inA&&!inB)};
        const double values[]={std::min(a,b),std::max(a,b),std::max(a,-b),-a,
            std::max(std::min(a,b),-std::max(a,-b))};
        for(int phase=1;phase<=5;++phase)
            csgCheck(sampler.phaseValue(phase,point)==values[phase-1]
                &&(sampler.phaseValue(phase,point)<0)==expected[phase-1],"CSG union/intersection/complement/intermediate oracle");
    }
    csgReject([&]{sampler.phaseValue(0,{1,1});},"zero physical phase ID accepted");
    csgReject([&]{sampler.phaseValue(6,{1,1});},"out-of-range physical phase ID accepted");
    csgReject([&]{sampler.levelSetValue(2,{1,1});},"out-of-range level set accepted");
    csgReject([&]{sampler.phaseValue(1,{std::numeric_limits<double>::infinity(),1});},"non-finite CSG point accepted");
    geometry.lines[0].points[0]={.5,.5};
    geometry.levelSetIndexA[0][0]={2,false};
    csgReject([&]{FittedDomainGeometry::InitialPhaseSampler invalid(mesh,geometry);},"forward CSG intermediate accepted");
    geometry.levelSetIndexA[0][0]={0,false};geometry.levelSetOperation[0][0]=LevelSetOperation::non;
    csgReject([&]{FittedDomainGeometry::InitialPhaseSampler invalid(mesh,geometry);},"non operation with different references accepted");
    geometry.levelSetOperation[0][0]=static_cast<LevelSetOperation>(999);
    csgReject([&]{FittedDomainGeometry::InitialPhaseSampler invalid(mesh,geometry);},"unknown CSG operation accepted");
}
DictionaryGeometry csgNozzle()
{
    auto root=std::filesystem::path(__FILE__);
    for(int i=0;i<4;++i)root=root.parent_path();
    UnitsScales scales;scales.L0_m=1e-6;
    const auto units=UnitsContext::MakeSI(scales);
    DictionaryGeometry geometry;geometry.readDict((root/"cases/electrospray/system/nozzle.csg").string(),&units);
    csgCheck(geometry.domainMax.x==5&&geometry.domainMax.y==16,"original nozzle units changed");
    return geometry;
}
void csgNozzleSelection(const QuadtreeMesh& mesh,const DictionaryGeometry& geometry,
    const FittedDomainGeometry::InitialPhaseSampler& sampler)
{
    const auto selection=FittedDomainGeometry::selectFixedSolids(mesh,geometry,sampler,{1});
    double xy=0,rz=0;std::set<int> solidBoundaryIds;
    for(int i=0;i<mesh.cellIndexNum();++i)
    {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);const auto point=cell.centerCoord();
        // TEST ONLY analytic reference; production uses the original CSG.
        const bool electrode=point.x>1.0&&point.y<3.0;
        csgCheck(selection.excludedPhaseByCell[i]==(electrode?1:0),"original nozzle selected-solid mask");
        if(!electrode){xy+=cell.area();rz+=2*3.14159265358979323846*point.x*cell.area();}
    }
    csgCheck(std::abs(xy-68)<2e-10&&std::abs(rz-328*3.14159265358979323846)<2e-9,
        "original nozzle XY/RZ selected-domain measure");
    FaceSegmentGraph graph;graph.build(mesh);
    csgCheck(selection.boundaryIdsBySegment.size()==graph.segments().size(),"native segment mapping size");
    for(std::size_t i=0;i<graph.segments().size();++i)
    {
        const auto& face=graph.segments()[i];const int k=mesh.cells(face.K).index();
        const int l=face.L<0?-1:mesh.cells(face.L).index();
        const bool exposed=!selection.excludedPhaseByCell[k]&&(l<0||selection.excludedPhaseByCell[l]);
        csgCheck((selection.boundaryIdsBySegment[i]>0)==exposed,"native exposed boundary coverage");
        if(!exposed)continue;
        int expected=0;
        if(l>=0){expected=face.isHorizontal?2:3;solidBoundaryIds.insert(expected);}
        else if(face.boundary==Boundary::south)expected=5;
        else if(face.boundary==Boundary::east)expected=8;
        else if(face.boundary==Boundary::north)expected=9;
        else expected=face.faceCenter.y<3?12:face.faceCenter.y<4?11:10;
        csgCheck(selection.boundaryIdsBySegment[i]==expected,"original nozzle native boundary ID");
    }
    csgCheck(solidBoundaryIds==std::set<int>({2,3}),"original fixed electrode boundaries were not both selected");
}
void csgCartesianDistances(const QuadtreeMesh& mesh,const DictionaryGeometry& geometry,
    const FittedDomainGeometry::InitialPhaseSampler& sampler)
{
    // Original physical wet wall, including the exact native7 failing site.
    for(double y:{0.,.025,.075,.1,.3,1.1,2.975,3.})
        csgCheck(sampler.levelSetValue(1,{1,y})==0&&sampler.phaseValue(2,{1,y})==0,
            "original Cartesian wet wall acquired an artificial distance/phase");
    const auto checkSides=[&](const FittedDomainGeometry::InitialPhaseSampler& query){
        for(double y:{.025,.075,1.1,2.975}){
            const double left=std::nextafter(1.,0.),right=std::nextafter(1.,2.);
            csgCheck(query.phaseValue(2,{left,y})==-(1.-left)&&query.phaseValue(2,{right,y})==right-1.,
                "one-ULP wet-wall displacement was made zero or lost its true CSG sign");}
    };
    checkSides(sampler);
    auto reversed=geometry;reversed.lines[2].reverse();
    for(auto& path:reversed.linesOnInterface)for(auto& ref:path)if(ref.value==2)ref.isNegative=!ref.isNegative;
    for(auto& path:reversed.signClosureLines)for(auto& ref:path)if(ref.value==2)ref.isNegative=!ref.isNegative;
    auto probeInput=csgMeshInput({5,16},40,80);QuadtreeMesh probe(probeInput);
    const FittedDomainGeometry::InitialPhaseSampler reverseSampler(probe,reversed);
    for(double y:{.025,.075,.3,1.1,2.975})
        csgCheck(reverseSampler.levelSetValue(1,{1,y})==0&&reverseSampler.phaseValue(2,{1,y})==0,
            "reversed Cartesian wet wall lost its exact zero");
    checkSides(reverseSampler);
    DictionaryGeometry horizontal;csgRectangle(horizontal,{0,1},{3,4});
    csgExpression(horizontal,{LevelSetOperation::non},{{0,false}},{{0,false}});
    auto horizontalInput=csgMeshInput({4,4},8,8);QuadtreeMesh horizontalMesh(horizontalInput);
    const FittedDomainGeometry::InitialPhaseSampler horizontalSampler(horizontalMesh,horizontal);
    auto horizontalReverse=horizontal;horizontalReverse.lines[0].reverse();horizontalReverse.linesOnInterface[0][0].isNegative=true;
    const FittedDomainGeometry::InitialPhaseSampler horizontalReverseSampler(horizontalMesh,horizontalReverse);
    for(double x:{0.,.025,.075,.3,1.1,2.975,3.}){
        csgCheck(horizontalSampler.phaseValue(1,{x,1})==0&&horizontalReverseSampler.phaseValue(1,{x,1})==0,
            "horizontal/reversed segment acquired an artificial distance/phase");
        if(x>0&&x<3){const double below=std::nextafter(1.,0.),above=std::nextafter(1.,2.);
            csgCheck(horizontalSampler.phaseValue(1,{x,below})==1.-below&&horizontalSampler.phaseValue(1,{x,above})==-(above-1.)&&
                horizontalReverseSampler.phaseValue(1,{x,below})==1.-below&&horizontalReverseSampler.phaseValue(1,{x,above})==-(above-1.),
                "one-ULP horizontal displacement was made zero or lost its CSG sign");}}
    csgCheck(std::abs(sampler.phaseValue(2,{1,-.025}))==.025&&
        std::abs(horizontalSampler.phaseValue(1,{-.025,1}))==.025,
        "Cartesian closest-point endpoint clamping changed physical distance");
    (void)mesh;
}
void csgPhysicalLiquidGasDistance(QuadtreeMesh& mesh,const DictionaryGeometry& definition,
    const FittedDomainGeometry::InitialPhaseSampler& sampler)
{
    using Query=std::function<ImmersedGeometry<QuadtreeMesh>::InitialCurveSample(Vector2D)>;
    std::vector<Query> queries;
    { // Both the temporary dictionary and geometry die before every query.
        auto geometry=definition;ImmersedGeometry<QuadtreeMesh> seed(geometry,mesh);
        std::map<int,std::set<int>> physical;
        for(const auto& boundary:definition.boundaries){
            if(boundary.ignored||boundary.kind==GeometricBoundary::domain)continue;
            if(!((boundary.phase==1&&boundary.opposedPhase==2)||(boundary.phase==2&&boundary.opposedPhase==1)))continue;
            physical[boundary.levelSetIndex].insert(boundary.lineIds.begin(),boundary.lineIds.end());
        }
        csgCheck(physical.size()==1&&physical.begin()->first==1&&physical.begin()->second==std::set<int>{3},
            "original physical phase pair does not isolate the L/G arc from its shared wet-wall level set");
        for(const auto& entry:physical){std::vector<SignInt> lines;
            for(const int line:entry.second)lines.push_back({line,false});
            queries.push_back(seed.initialCurveSampler(mesh,lines,entry.first));
        }
    }
    const auto distance=[&](Vector2D point){double result=std::numeric_limits<double>::infinity();
        for(const auto& query:queries)result=std::min(result,query(point).unsignedDistance);return result;};
    // Independent TEST-only circular-arc oracle and the original exact endpoints.
    csgCheck(distance({1,3})==0&&distance({0,4})<2e-14,
        "physical L/G distance omitted the pin or true-axis endpoint");
    for(const Vector2D point:std::vector<Vector2D>{{.5,3.5},{.8,3.8},{.625,3+std::sqrt(1-.625*.625)}}){
        const double expected=std::abs(std::hypot(point.x,point.y-3)-1);
        csgCheck(std::abs(distance(point)-expected)<2e-14,"selected physical L/G arc distance differs from independent circle");
        if(distance(point)>1e-13)csgCheck(queries.front()(point).negative==(sampler.phaseValue(2,point)<0),
            "selected distance lines changed the original complete level-set parity/CSG sign");
    }
    csgCheck(sampler.phaseValue(2,{1,1.5})==0&&distance({1,1.5})==1.5&&
        std::abs(sampler.phaseValue(2,{.99,1.5}))<.011&&distance({.99,1.5})>1.4,
        "wet solid wall was incorrectly treated as the nearby L/G initial-band curve");
    csgCheck(distance({.01,2})>1&&sampler.phaseValue(2,{.01,2})<0,
        "sign-only closure was incorrectly used as physical L/G distance");
}
void csgFittedSolid()
{
    csgBooleanValues();auto geometry=csgNozzle();
    auto input=csgMeshInput({5,16},40,80,1);QuadtreeMesh mesh(input);
    const FittedDomainGeometry::InitialPhaseSampler sampler(mesh,geometry);
    csgCartesianDistances(mesh,geometry,sampler);
    csgPhysicalLiquidGasDistance(mesh,geometry,sampler);
    csgNozzleSelection(mesh,geometry,sampler);
    auto coarseInput=csgMeshInput({5,16},20,40,1);QuadtreeMesh cutting(coarseInput);
    const FittedDomainGeometry::InitialPhaseSampler cuttingSampler(cutting,geometry);
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(cutting,geometry,cuttingSampler,{1});},"original z=3 cutting base mesh accepted");
    DictionaryMesh::RefineBox box;box.min={0,0};box.max={5,16};box.level=1;
    cutting.applyRefineBoxes({box});csgNozzleSelection(cutting,geometry,cuttingSampler);
    // Actual refinement followed by actual coarsening, with no evolving LS.
    mesh.applyRefineBoxes({box});csgNozzleSelection(mesh,geometry,sampler);
    std::vector<VertField<double>> none;std::vector<VertField<double>*> vertices;
    std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
    std::vector<CellField<double>*> cells;CellField<double> gradient(mesh,0.0),target(mesh,0.0);
    const int before=mesh.cellIndexNum();mesh.adapt(none,gradient,vertices,horizontal,vertical,cells,&target);
    csgCheck(mesh.cellIndexNum()<before,"actual fitted mesh coarsening did not occur");
    csgNozzleSelection(mesh,geometry,sampler);
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,geometry,sampler,{0});},"invalid selected solid phase accepted");
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,geometry,sampler,{4});},"out-of-range selected solid phase accepted");
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,geometry,sampler,{1,1});},"duplicate solid phase accepted");
    auto changed=geometry;changed.lines[2].points[0].x+=.125;
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,changed,sampler,{1});},"changed native geometry reused stale CSG certificate");
    auto changedCsg=geometry;changedCsg.levelSetIndexA[0][0].isNegative=!changedCsg.levelSetIndexA[0][0].isNegative;
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,changedCsg,sampler,{1});},"changed CSG reused stale certificate");
    auto changedClosure=geometry;changedClosure.signClosureLines[0][0].isNegative=!changedClosure.signClosureLines[0][0].isNegative;
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,changedClosure,sampler,{1});},"changed sign closure reused stale certificate");
    auto missing=geometry;
    missing.boundaries.erase(std::remove_if(missing.boundaries.begin(),missing.boundaries.end(),
        [](const auto& b){return b.id==2;}),missing.boundaries.end());
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,missing,sampler,{1});},"incomplete native solid coverage accepted");
    auto ambiguous=geometry;auto duplicate=geometry.boundaries.back();
    for(const auto& boundary:geometry.boundaries)if(boundary.id==3)duplicate=boundary;
    duplicate.id=99;ambiguous.boundaries.push_back(duplicate);
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,ambiguous,sampler,{1});},"ambiguous native solid coverage accepted");
    auto moving=geometry;
    for(auto& boundary:moving.boundaries)if(boundary.id==3)boundary.kind=GeometricBoundary::moving;
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,moving,sampler,{1});},"moving selected-solid boundary accepted");
    // A complete closed solid lies away from every center and corner of the
    // one-cell mesh: whole-leaf geometry, rather than samples, must reject it.
    DictionaryGeometry hidden;csgRectangle(hidden,{.1,.1},{.2,.2});
    csgExpression(hidden,{LevelSetOperation::non},{{0,false}},{{0,false}});
    csgExpression(hidden,{LevelSetOperation::non},{{0,true}},{{0,true}});
    auto tinyInput=csgMeshInput({1,1},1,1);QuadtreeMesh tiny(tinyInput);
    const FittedDomainGeometry::InitialPhaseSampler hiddenSampler(tiny,hidden);
    csgCheck(hiddenSampler.phaseAt({.5,.5})==2,"hidden solid fixture center");
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(tiny,hidden,hiddenSampler,{1});},"hidden closed solid passed sample-only certificate");
    // Curved solid selection remains outside this fixed whole-leaf capability.
    csgReject([&]{FittedDomainGeometry::selectFixedSolids(mesh,geometry,sampler,{2});},"curved selected liquid accepted as fitted solid");
    std::cout<<"Initial CSG Boolean algebra/ownership and fitted solid XY/RZ native boundary contracts passed\n";
}
}

int main(int argc,char** argv) try {
    if(argc==2&&std::string(argv[1])=="--csg-fitted-solid-only"){csgFittedSolid();return 0;}
    readerFailures();
    extensions();
    blockScope();
    Reader reader;
    for (bool brace : {false,true}) {
        const auto parse=[&](const std::string& text, std::vector<SignInt>& result) {
            return brace ? reader.parseIndex_brace("{"+text+"}",result)
                         : reader.parseIndex("("+text+")",result);
        };
        for (const auto text : {"bad","1; 2","1, -","(1,2)","1.5","1,2 suffix"}) {
            std::vector<SignInt> values;
            if (parse(text,values)) throw std::runtime_error("malformed indices accepted");
        }
        std::vector<SignInt> values;
        if (!parse("1, -12  3",values) || values.size()!=3 || values[0].value!=1
            || values[0].isNegative || values[1].value!=12 || !values[1].isNegative
            || values[2].value!=3 || values[2].isNegative)
            throw std::runtime_error("valid signed index list changed");
    }
    std::cout<<"Malformed index lists return instead of hanging; signed lists preserved\n";
    return 0;
} catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
