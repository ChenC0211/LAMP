#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Modules/Hydrodynamics/IO/HydroGeometryLoader.h"
#include "Modules/Hydrodynamics/IO/DictionaryHydrodynamics.h"
#include "Modules/Hydrodynamics/HydroBoundaryVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceBoundary.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/Geometry/LevelSetGaussTrace.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedMeshTransfer.h"
#include "Modules/Hydrodynamics/HydroFittedRefinement.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <cmath>
#include <iostream>
#include <type_traits>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <limits>
#include <optional>

void hydroSplineTimeContract(bool fullInterface);
void hydroCompatiblePressureContract();
void hydroCompatibleMomentsContract();
void hydroMomentViscosityContract();
void hydroMomentRZViscosityContract();
void hydroMomentRZTransportContract();
void hydroMomentRZReconstructionContract();
void hydroMomentRZTimeContract();
void hydroMomentRZManufacturedContract();
void hydroMomentRZCoordinatesContract();
void hydroMomentTimeContract();
void hydroMomentPortsContract();
void hydroMomentPortSplineContract();
void hydroParametricSupplyContract();
void hydroMomentSupplyContract();
void hydroMomentComponentsContract();
void hydroParametricComponentsTransferContract();
void hydroParametricNeckContract();
void hydroMomentNeckContract();
void hydroMomentNeckReadinessContract();
void hydroMomentTransportContract();
void hydroMomentNSManufacturedContract();
void hydroMomentPhysicalInputContract();
void hydroMomentMaterialContract();
void hydroMomentMaterialFlowContract();
void hydroMomentSplineContract();
void hydroMomentAmrSplineContract();
void hydroMomentAmrTwoPhaseContract();
void hydroMomentMeshTransferContract();
void hydroMomentPortMeshTransferContract();
void hydroMomentNativeRefinementContract();
void hydroMomentNativeMaterialContract();
void hydroMomentCartesianMaterialGramContract();
void hydroMomentNativeStateContract();
void hydroMomentPhysicalComplementContract();
void hydroMomentBlockSolveContract();
void hydroMomentNativeFluxContract();
void hydroMomentMassCurrentContract();
void hydroMomentNativeMomentumContract();
void hydroMomentNativeCoupledContract();
void hydroMomentNativeRemeshContract();
void hydroMomentNativeAccuracyContract();
void hydroMomentNativeDensityTransportContract();
void hydroParametricWeightedPhaseContract();
void hydroMomentLocalResolutionContract();
void hydroMomentSplineInterfaceContract();
void hydroMomentSplineFluxRootsContract();
void hydroMomentPositiveFluxContract();
void hydroMomentTwoPhaseContract();
void hydroMomentCapillaryContract();
void hydroMomentMaterialSpaceContract();
void hydroMomentRZSplineContract();
void hydroMomentRZInterfaceContract(bool motion=true);
void hydroMomentRZTwoPhaseContract();
void hydroMomentRZFlowReplay(const char* path);
void hydroMomentFreeSlipPhaseProbe(const char* path,const std::string& task="velocity");
void hydroMomentCapillaryMetricContract();
void hydroMomentRZKineticContract();
void hydroParametricMassContract();
void hydroParametricMeasureContract();
void hydroParametricMaterialSupportContract();
void hydroParametricCurvatureContract();
void hydroParametricNonlinearCurvatureContract();
void hydroParametricAmrContract();
void hydroParametricCapillaryContract(int steps);
void hydroParametricMidpointContract(int steps);
void hydroParametricFluxReplay(const char* path);
void hydroEnergyViscosityContract();
namespace {
void require(bool valid, const char* reason) { if (!valid) throw std::runtime_error(reason); }
#define NO_HYDRO_RUNTIME_FIELD(field) \
    template<class T,class=void> struct has_runtime_##field : std::false_type {}; \
    template<class T> struct has_runtime_##field<T,std::void_t<decltype(std::declval<T>().field)>> : std::true_type {}; \
    static_assert(!has_runtime_##field<HydrodynamicsConfig>::value,"runtime fields must not reenter physical config");
NO_HYDRO_RUNTIME_FIELD(steps)
NO_HYDRO_RUNTIME_FIELD(dt)
NO_HYDRO_RUNTIME_FIELD(outputName)
NO_HYDRO_RUNTIME_FIELD(writeFields)
NO_HYDRO_RUNTIME_FIELD(writeInitialFields)
NO_HYDRO_RUNTIME_FIELD(writeEvery)
NO_HYDRO_RUNTIME_FIELD(logEvery)
#undef NO_HYDRO_RUNTIME_FIELD

template<class T,class=void> struct hasInterfaceSettings : std::false_type {};
template<class T> struct hasInterfaceSettings<T,std::void_t<decltype(std::declval<T>().advectInterface)>> : std::true_type {};
template<class T,class=void> struct hasInterfaceBoundary : std::false_type {};
template<class T> struct hasInterfaceBoundary<T,std::void_t<decltype(std::declval<T>().interfaceBoundary)>> : std::true_type {};
template<class T,class=void> struct hasRedistance : std::false_type {};
template<class T> struct hasRedistance<T,std::void_t<decltype(std::declval<T>().reinitializeEvery)>> : std::true_type {};
static_assert(!hasInterfaceSettings<HydrodynamicsConfig>::value && !hasInterfaceBoundary<HydrodynamicsConfig>::value
    && !hasRedistance<HydrodynamicsConfig>::value,"Interface input must not reenter Hydro physics");

void interfaceInputContract()
{
    namespace fs=std::filesystem;
    const auto path=fs::current_path()/("nodal_input_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count())+".dict");
    const auto write=[&](const std::string& text) {
        std::ofstream stream(path);stream<<text;stream.close();require(bool(stream),"interface input fixture write failed");
    };
    const std::string body="steps 3; dt 0.001; coordinateSystem planar2D; materialModel twoPhase; "
        "interfaceLevelSet liquid; densityNegative 1; densityPositive 2; viscosityNegative 0.01; viscosityPositive 0.02; "
        "interfaceThickness 0.08; interfaceForceGeometry hermite; interfaceForceLayer regularized; ";
    DictionaryHydrodynamics input;
    write("Hydrodynamics { "+body+"}");
    require(input.readDict(path.string()) && !input.interfaceConfig().advect && input.interfaceConfig().reinitializeEvery==0
        && input.interfaceConfig().boundary.model.empty(),"legacy Interface defaults changed");
    write("Hydrodynamics { "+body+"advectInterface true; reinitializeEvery 2; }");
    require(input.readDict(path.string()) && input.interfaceConfig().advect && input.interfaceConfig().reinitializeEvery==2
        && input.config().densityPositive==2.,"transport and material settings must be separately retained");
    const std::string contact="Hydrodynamics { "+body+"advectInterface true; reinitializeEvery 0; "
        "boundaryPatches { supply { side south; type normalVelocity; start 0; end 1; value -0.01; } } "
        "interfaceBoundary { model prescribedScalar; fixedContacts ((0.25, 0)); "
        "inflow { supply { offset 0.5; slope -2; origin 0.5; } } } }";
    write(contact);require(input.readDict(path.string()),"prescribed scalar input failed");
    const auto& boundary=input.interfaceConfig().boundary;
    require(boundary.model=="prescribedScalar" && boundary.fixedContacts.size()==1
        && boundary.fixedContacts[0].x==.25 && boundary.inflow.size()==1
        && boundary.inflow[0].patch=="supply","Interface must own original contacts and inlet data");
    auto conditions=HydroInterfaceBoundarySupport::scalarConditions(boundary,input.config().boundaryPatches);
    require(conditions.fixedZeroVertices.size()==1 && conditions.inflow.size()==1,
        "Hydro patch binding lost scalar constraints");
    // The callback owns its scalar values; a later successful read must not change it.
    write("Hydrodynamics { "+body+"}");
    require(input.readDict(path.string()) && !input.interfaceConfig().advect && input.interfaceConfig().boundary.inflow.empty()
        && input.interfaceConfig().boundary.fixedContacts.empty(),"reread must reset all Interface settings");
    require(!input.readDict(path.string()+".missing") && input.config().densityPositive==2.,
        "missing file must retain the original reader semantics");
    require(conditions.inflow[0].value({.25,0},0.)==1. && conditions.inflow[0].value({.75,0},10.)==0.,
        "resolved affine inlet must preserve coordinates, values and owning lifetime");
    write(contact);require(input.readDict(path.string()),"contact reread failed");
    auto patches=input.config().boundaryPatches;patches[0].value=0.;
    bool rejected=false;
    try { HydroInterfaceBoundarySupport::scalarConditions(input.interfaceConfig().boundary,patches); }
    catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"scalar binding must still reject an inlet without inward velocity");
}

void runInputContract()
{
    namespace fs=std::filesystem;
    const auto path=fs::current_path()/("hydro_run_input_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count())+".dict");
    const auto write=[&](const std::string& body) {
        std::ofstream stream(path);stream<<"Hydrodynamics { "<<body<<" viscosity 0.01;";
        if(body.find("density ")==std::string::npos)stream<<" density 1;";
        stream<<" }\n";stream.close();
        require(bool(stream),"run input fixture write failed");
    };
    DictionaryHydrodynamics input;
    write("steps 3; dt 0.125; coordinateSystem planar2D;");
    require(input.readDict(path.string()),"standalone Hydro reader could not open explicit path");
    const auto& initial=input.runInput();
    require(initial.steps==3 && initial.dt==.125 && initial.outputName=="hydro" && initial.writeFields
        && !initial.writeInitialFields && initial.writeEvery==1 && initial.logEvery==1,
        "legacy run defaults changed while splitting physics");
    write("steps 4; dt 0.02; coordinateSystem planar2D; density 4.5; outputName custom_42; writeFields false; writeInitialFields true; writeEvery 7; logEvery 2;");
    require(input.readDict(path.string()),"custom run input failed");
    const auto& custom=input.runInput();
    require(custom.steps==4 && custom.dt==.02 && custom.outputName=="custom_42" && !custom.writeFields
        && custom.writeInitialFields && custom.writeEvery==7 && custom.logEvery==2
        && input.config().density==4.5,"run settings must remain separate from physical values");
    write("steps 1; dt 0.001; coordinateSystem planar2D;");
    require(input.readDict(path.string()) && input.runInput().outputName=="hydro" && input.config().density==1.,
        "successful reread must reset both configuration groups");
    require(!input.readDict(path.string()+".missing") && input.runInput().steps==1,
        "missing dictionary retains original false-return semantics");
    const auto reject=[&](const std::string& body,const char* expected) {
        write(body);std::string message;
        try { input.readDict(path.string()); } catch(const std::runtime_error& e) { message=e.what(); }
        require(message.find(expected)!=std::string::npos,"run input validation changed its original reason");
    };
    reject("steps 0; dt 0.1; coordinateSystem planar2D;","steps must be a positive integer");
    reject("steps 2147483647; dt 0.1; coordinateSystem planar2D;","steps must be a positive integer");
    reject("steps 3; dt 0; coordinateSystem planar2D;","dt must be positive with a finite end time");
    reject("steps 3; dt 1e308; coordinateSystem planar2D;","dt must be positive with a finite end time");
    reject("steps 3; dt 0.1; coordinateSystem planar2D; writeEvery 0;","writeEvery and logEvery must be positive");
    reject("steps 3; dt 0.1; coordinateSystem planar2D; logEvery 0;","writeEvery and logEvery must be positive");
    reject("steps 3; dt 0.1; coordinateSystem planar2D; outputName bad/name;","outputName must be a nonempty filename stem");
    write("steps 3; dt 0.1; coordinateSystem planar2D; linearResidualTolerance 2e-5;");
    require(input.readDict(path.string()),"legacy tolerance input failed");
    const auto old=input.config().conservation;
    require(old.momentumBalance==2e-5 && old.initialContinuity==2e-5 && old.outletAdmissibility==2e-5
        && old.boundaryFluxBalance==2e-5,"legacy reader did not expand the original acceptance tolerance");
    reject("steps 3; dt 0.1; coordinateSystem planar2D; Conservation { momentumBalance 1e-8; }",
        "Conservation requires separate Hydro inputs");
    const std::string flow="coordinateSystem planar2D; linearResidualTolerance 2e-5; ";
    const std::string ports="boundaryPatches { outlet { side north; type pressureOutlet; start 0; end 1; value 0; } } ";
    const auto readSeparate=[&](const std::string& text) { write(flow+text);return readNodalHydrodynamicsConfig(path.string(),{}); };
    const auto defaults=readSeparate("");
    require(defaults.momentConvection=="central","Compatible-Moment convection must default to central");
    require(defaults.linearResidualTolerance==2e-5 && defaults.conservation.momentumBalance==1e-8
        && defaults.conservation.initialContinuity==1e-8 && defaults.conservation.outletAdmissibility==1e-8
        && defaults.conservation.boundaryFluxBalance==1e-8,"formal conservation defaults inherited linear accuracy");
    const auto explicitValues=readSeparate(ports+"Conservation { momentumBalance 2e-9; initialContinuity 3e-9; outletAdmissibility 4e-9; boundaryFluxBalance 5e-9; }");
    require(explicitValues.linearResidualTolerance==2e-5 && explicitValues.conservation.momentumBalance==2e-9
        && explicitValues.conservation.initialContinuity==3e-9 && explicitValues.conservation.outletAdmissibility==4e-9
        && explicitValues.conservation.boundaryFluxBalance==5e-9,"formal conservation settings lost their separate consumers");
    for(const auto& bad: {"momentumBalance 0;", "momentumBalance -1;", "momentumBalance 1;",
        "momentumBalance 1e-8; momentumBalance 1e-7;", "unknown 1e-8;", "initialContinuity 1e-8;",
        "outletAdmissibility 1e-8;", "boundaryFluxBalance 1e-8;"}) {
        bool rejected=false;
        try { readSeparate(std::string("Conservation { ")+bad+" }"); } catch(const std::exception&) { rejected=true; }
        require(rejected,"invalid or unused formal conservation control was accepted");
    }
}
void domainInputContract()
{
    namespace fs=std::filesystem;
    const auto cwd=fs::current_path();
    const auto directory=cwd/("hydro_domain_input_"+std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    fs::create_directories(directory);
    const auto geometryPath=(directory/"shape.input").generic_string();
    const auto meshPath=(directory/"grid.input").generic_string();
    const auto write=[](const std::string& path,const std::string& text) {
        std::ofstream stream(path);stream<<text;stream.close();
        require(bool(stream),"input fixture write failed");
    };
    const std::string domain="Domain { min (0, 0); max (2, 3); dimension 2; }\n";
    const std::string mesh="adaptiveMesh false\nBackgroundMesh { divide (8, 12); }\nQuadtreeMesh { maxCellNum 256; maxInitialIterations 0; }\n";
    write(geometryPath,domain);write(meshPath,mesh);
    HydrodynamicsConfig config;DictionaryGeometry geometry;DictionaryMesh grid;
    require(HydroGeometryInput::load(config,{},geometry,grid,geometryPath,meshPath)==-1,
        "constant material must have no interface");
    require(geometry.hasDomain && geometry.levelSetName.empty() && grid.Nx==8 && grid.Ny==12
        && grid.maxCellNum==256 && grid.maxInitialIterations==0 && grid.minCoord.x==0.
        && grid.maxCoord.x==2. && grid.maxCoord.y==3.,"geometry and mesh preflight must share bounds");
    const std::string contour="Point { (0.5, 0.5) (1.5, 0.5) (1, 2) }\n"
        "Line { (1, 2) { type straight; } (2, 3) { type straight; } (3, 1) { type straight; } }\n"
        "LevelSet { liquid { lines (1, 2, 3); innerPoint (1, 1); moving true; } }\n";
    write(geometryPath,domain+contour);config.materialModel="twoPhase";config.interfaceLevelSet="liquid";
    DictionaryGeometry liquid;DictionaryMesh liquidGrid;
    require(HydroGeometryInput::load(config,{},liquid,liquidGrid,geometryPath,meshPath)==0,
        "named material interface must retain its index");
    require(liquid.levelSetName==std::vector<std::string>{"liquid"} && liquid.lines.size()==3
        && liquid.linesOnInterface[0].size()==3 && liquid.moving[0] && liquidGrid.Nx==8,
        "module must publish the original contour and mesh records");
    const auto rejected=[&](const std::string& shape,const std::string& gridText,const std::string& marker,
                            bool geometryFirst=false) {
        write(geometryPath,shape);write(meshPath,gridText);
        DictionaryGeometry rejectedGeometry;DictionaryMesh rejectedMesh;rejectedMesh.Nx=123;
        std::string message;
        try { HydroGeometryInput::load(config,{},rejectedGeometry,rejectedMesh,geometryPath,meshPath); }
        catch(const std::runtime_error& error) { message=error.what(); }
        require(message.find(marker)!=std::string::npos,"module input rejection lost the original reason");
        require(message.find(geometryFirst?geometryPath:meshPath)!=std::string::npos,
            "module input rejection must name the supplied file");
        if(geometryFirst)require(rejectedMesh.Nx==123,"geometry must reject before touching mesh input");
    };
    config=HydrodynamicsConfig{};
    rejected("Domain { min (0, 0); max (0, 3); }","bad mesh", "max greater than min",true);
    rejected(domain+contour,mesh,"constant materialModel requires Domain-only geometry",true);
    rejected(domain,"adaptiveMesh true","adaptiveMesh must be false");
    rejected(domain,"adaptiveMesh false\nBackgroundMesh { divide (3, 8); }\nQuadtreeMesh { maxCellNum 100; }",
        "divide dimensions must be integers in [4, 256]");
    rejected(domain,mesh+"BackgroundMesh {}","duplicate key/block 'BackgroundMesh'");
    rejected(domain,"adaptiveMesh false\nBackgroundMesh { divide (8, 8); }\nQuadtreeMesh { maxCellNum 10; }",
        "maxCellNum must cover");
    rejected(domain,"adaptiveMesh false\nBackgroundMesh { divide (8, 8); }\nQuadtreeMesh { maxCellNum 100; maxInitialIterations 1; }",
        "initialization AMR is unsupported");
    config.materialModel="twoPhase";config.interfaceLevelSet="liquid";config.coordinateSystem="axisymmetricRZ";
    config.interfaceForceGeometry="conservative";config.interfaceForceLayer="inventory";
    InterfaceNodalConfig conservative;conservative.transport="conservative";conservative.advect=true;
    const std::string columnText="AxisymmetricColumn { levelSet liquid; radius .6; relativeAmplitude .001; }";
    write(geometryPath,domain+columnText);write(meshPath,mesh);
    std::optional<HydroGeometryInput::AxisymmetricColumn> column;
    DictionaryGeometry columnGeometry;DictionaryMesh columnMesh;
    require(HydroGeometryInput::load(config,conservative,columnGeometry,columnMesh,geometryPath,meshPath,&column)==-1
        && column && column->radius==.6 && column->relativeAmplitude==.001 && columnGeometry.lines.empty()
        && columnGeometry.levelSetName.empty(),"analytic column created a second auxiliary contour");
    const auto badColumn=[&](std::string text,HydrodynamicsConfig flow,InterfaceNodalConfig interface,bool consumer=true) {
        write(geometryPath,domain+text);DictionaryGeometry shape;DictionaryMesh grid;
        std::optional<HydroGeometryInput::AxisymmetricColumn> result;bool failed=false;
        try{HydroGeometryInput::load(flow,interface,shape,grid,geometryPath,meshPath,consumer?&result:nullptr);}
        catch(const std::exception&){failed=true;}
        require(failed,"unsupported analytic column input was accepted");
    };
    badColumn(columnText,config,conservative,false);
    auto wrongFlow=config;wrongFlow.coordinateSystem="planar2D";badColumn(columnText,wrongFlow,conservative);
    wrongFlow=config;wrongFlow.boundaryPatches={{"inlet",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,2,0}};
    badColumn(columnText,wrongFlow,conservative);
    auto wrongInterface=conservative;wrongInterface.transport="nodal";badColumn(columnText,config,wrongInterface);
    badColumn("AxisymmetricColumn { levelSet liquid; radius .6; relativeAmplitude 1; }",config,conservative);
    badColumn("AxisymmetricColumn { levelSet liquid; radius 2; relativeAmplitude 0; }",config,conservative);
    badColumn(columnText+contour,config,conservative);
    require(fs::current_path()==cwd && !fs::exists(directory/"output") && !fs::exists(directory/"cache"),
        "input service must neither change working directory nor publish physical output");
}
template<class Action> void rejects(Action action) {
    bool rejected = false;
    try { action(); } catch (const std::exception&) { rejected = true; }
    require(rejected, "invalid Hydro query was accepted");
}
template<class T, class = void> struct exposesPressure : std::false_type {};
template<class T> struct exposesPressure<T, std::void_t<decltype(std::declval<T&>().p)>> : std::true_type {};
template<class T, class = void> struct exposesVelocity : std::false_type {};
template<class T> struct exposesVelocity<T, std::void_t<decltype(std::declval<T&>().ux)>> : std::true_type {};
template<class T, class = void> struct exposesHistory : std::false_type {};
template<class T> struct exposesHistory<T, std::void_t<decltype(std::declval<T&>().uxPre)>> : std::true_type {};
static_assert(!exposesPressure<Hydrodynamics>::value && !exposesVelocity<Hydrodynamics>::value
    && !exposesHistory<Hydrodynamics>::value, "Hydro arrays must remain private");
static_assert(!std::is_reference_v<decltype(std::declval<HydroStateView>().velocityX[0])>,
    "read-only state must not expose writable references");
static_assert(!std::is_reference_v<decltype(std::declval<HydroStateView>().vertexDensity[0])>,
    "read-only material output must not expose writable references");
// Fixed-dt kernels must remain reachable only through the managed state/time
// contract. Probe every retired overload, while preserving the distinct old
// nodal-coupling API; these are public-interface checks, not solver replicas.
#define HYDRO_CALL_PROBE(name) \
    template<class T, class... Args> auto public_##name(int) \
        -> decltype(std::declval<T&>().name(std::declval<Args>()...), std::true_type{}); \
    template<class, class...> std::false_type public_##name(...);
HYDRO_CALL_PROBE(solve)
HYDRO_CALL_PROBE(solveRZ)
HYDRO_CALL_PROBE(solveFittedRZ)
HYDRO_CALL_PROBE(updateProperties)
HYDRO_CALL_PROBE(initializeConstantProperties)
HYDRO_CALL_PROBE(assignFittedMaterials)
#undef HYDRO_CALL_PROBE
using H=Hydrodynamics;using M=QuadtreeMesh&;using O=const HydroSolveOptions&;
using C=const CellField<double>&;using X=const VEdgeField<double>&;using Y=const HEdgeField<double>&;
static_assert(!decltype(public_solve<H,M,double,O>(0))::value
    && !decltype(public_solve<H,M,double,O,C,C>(0))::value
    && !decltype(public_solve<H,M,double,O,C,C,X,Y>(0))::value);
static_assert(!decltype(public_solveRZ<H,M,double,O>(0))::value
    && !decltype(public_solveRZ<H,M,double,O,C,C>(0))::value
    && !decltype(public_solveRZ<H,M,double,O,C,C,X,Y>(0))::value);
static_assert(!decltype(public_solveFittedRZ<H,M,double,O,const HydroFittedDomain&>(0))::value);
static_assert(!decltype(public_initializeConstantProperties<H,M,double,double>(0))::value
    && !decltype(public_updateProperties<H,M,LevelSet&,double,double,double,double,double>(0))::value
    && !decltype(public_assignFittedMaterials<H,M,const HydroFittedDomain&,const std::vector<int>&,
        const HydroPhaseProperties&,const HydroPhaseProperties&>(0))::value);
static_assert(decltype(public_solve<H,M,double,LegacyNodalState&,ImmersedGeometry<QuadtreeMesh>&>(0))::value
    && decltype(public_solveRZ<H,M,double,LegacyNodalState&,ImmersedGeometry<QuadtreeMesh>&>(0))::value
    && decltype(public_updateProperties<H,M,LevelSet&>(0))::value);
void equal(Vector2D a, Vector2D b) {
    require(a.x == b.x && a.y == b.y, "velocity reconstruction changed");
}

template<class A, class B> void equalField(const A& a, const B& b) {
    require(a.size()==b.size(), "candidate field dimensions changed");
    for(int i=0;i<a.size();++i)
        require(a.label(i)==b.label(i) && a[i]==b[i], "candidate changed numerical state/history");
}
void equalState(const HydroStateView& a, const HydroStateView& b) {
    equalField(a.physicalPressure,b.physicalPressure);equalField(a.density,b.density);equalField(a.viscosity,b.viscosity);
    equalField(a.velocityX,b.velocityX);equalField(a.previousVelocityX,b.previousVelocityX);equalField(a.predictorX,b.predictorX);
    equalField(a.velocityY,b.velocityY);equalField(a.previousVelocityY,b.previousVelocityY);equalField(a.predictorY,b.predictorY);
    equalField(a.divergence,b.divergence);equalField(a.speed,b.speed);
    equalField(a.vertexDensity,b.vertexDensity);equalField(a.vertexViscosity,b.vertexViscosity);
    equalField(a.externalForceX,b.externalForceX);equalField(a.externalForceY,b.externalForceY);
}

void velocityInitializationContract(bool axisymmetric, bool refined) {
    DictionaryMesh input; input.Nx=8; input.Ny=9; input.maxCellNum=512; input.maxInitialIterations=0;
    input.minCoord={axisymmetric?0.:-2.,3.}; input.maxCoord={axisymmetric?3.:1.,5.};
    QuadtreeMesh mesh(input);
    if(refined) {
        DictionaryMesh::RefineBox box;
        box.min={input.minCoord.x+1.,3.6}; box.max={input.minCoord.x+2.,4.4}; box.level=1;
        mesh.applyRefineBoxes({box}); require(mesh.cellIndexNum()>72,"initial field test did not refine");
    }
    const auto geometry=HydroStepGeometry::complete(axisymmetric);
    // Prescribed physical pressure survives both velocity profiles; invalid
    // pressure is rejected before publishing either pressure or histories.
    CellField<double> initialPressure(mesh,0.0);
    for(int i=0;i<initialPressure.size();++i) initialPressure.value[i]=.25+.001*i;
    for(auto profile:{HydroInitialVelocity::Zero,HydroInitialVelocity::Smooth}) {
        Hydrodynamics withPressure(mesh);
        initialPressure.value[0]=std::numeric_limits<double>::infinity();
        rejects([&]{withPressure.initializeCompleteVelocity(mesh,geometry,profile,.125,&initialPressure);});
        require(!withPressure.state().initialStateSet,"invalid pressure published initial state");
        initialPressure.value[0]=.25;
        withPressure.initializeCompleteVelocity(mesh,geometry,profile,.125,&initialPressure);
        for(int i=0;i<initialPressure.size();++i)
            require(withPressure.state().physicalPressure[i]==initialPressure.value[i],"initial profile discarded physical pressure");
    }
    Hydrodynamics hydro(mesh), doubled(mesh), reference(mesh), zero(mesh), zeroAmplitude(mesh);
    for(double bad:{-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
        for(auto profile:{HydroInitialVelocity::Zero,HydroInitialVelocity::Smooth})
            rejects([&]{hydro.initializeCompleteVelocity(mesh,geometry,profile,bad);});
    rejects([&]{hydro.initializeCompleteVelocity(mesh,geometry,static_cast<HydroInitialVelocity>(99),.125);});
    auto invalidGeometry=geometry; invalidGeometry.kind=HydroDomainKind::FittedRZ;
    rejects([&]{hydro.initializeCompleteVelocity(mesh,invalidGeometry,HydroInitialVelocity::Zero,.125);});
    DictionaryMesh otherInput=input; otherInput.Nx=4; QuadtreeMesh other(otherInput);
    rejects([&]{hydro.initializeCompleteVelocity(other,geometry,HydroInitialVelocity::Smooth,.125);});
    require(!hydro.state().initialStateSet && !hydro.state().progress.stopped,
        "invalid initial profile published or stopped state");
    hydro.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Smooth,.125);
    reference.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Smooth,.125);
    doubled.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Smooth,.25);
    zero.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Zero,1.e100);
    zeroAmplitude.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Smooth,0.);
    equalState(zero.state(),zeroAmplitude.state());
    const auto state=hydro.state(), twice=doubled.state(), rest=zero.state();
    require(state.initialStateSet && state.progress.geometry==geometry.kind
        && state.progress.completedSteps==0 && state.progress.acceptedSteps==0
        && state.progress.endTime==0. && state.progress.acceptedEndTime==0.
        && !state.progress.pendingAcceptance && !state.progress.stopped && !state.materials,
        "prescribed initialization solved, accepted, or advanced a physical step");
    equalField(state.velocityX,state.previousVelocityX); equalField(state.velocityX,state.predictorX);
    equalField(state.velocityY,state.previousVelocityY); equalField(state.velocityY,state.predictorY);
    bool nonzeroX=false,nonzeroY=false,nonzeroDivergence=false,peak=false;
    for(int i=0;i<state.velocityX.size();++i) {
        const auto& edge=mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        require(twice.velocityX[i]==2.*state.velocityX[i] && rest.velocityX[i]==0.,"x amplitude scaling changed");
        if(edge.boundary()!=Boundary::non) require(state.velocityX[i]==0.,"initial x boundary is nonzero");
        nonzeroX |= state.velocityX[i]!=0.;
        const auto point=edge.midpointCoord();
        if(!refined && point.x==.5*(input.minCoord.x+input.maxCoord.x) && point.y==4.) {
            require(state.velocityX[i]==.125,"smooth profile lost its center peak on a translated rectangle"); peak=true;
        }
    }
    for(int i=0;i<state.velocityY.size();++i) {
        const auto& edge=mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        require(twice.velocityY[i]==2.*state.velocityY[i] && rest.velocityY[i]==0.,"y amplitude scaling changed");
        if(edge.boundary()!=Boundary::non) require(state.velocityY[i]==0.,"initial y boundary is nonzero");
        nonzeroY |= state.velocityY[i]!=0.;
    }
    for(int i=0;i<state.physicalPressure.size();++i) {
        require(state.physicalPressure[i]==0. && std::isfinite(state.divergence[i]) && std::isfinite(state.speed[i])
            && rest.physicalPressure[i]==0. && rest.divergence[i]==0. && rest.speed[i]==0.,
            "initial pressure/diagnostics are invalid");
        nonzeroDivergence |= state.divergence[i]!=0.;
    }
    require(nonzeroX && nonzeroY && nonzeroDivergence && (refined || peak),
        "smooth seed must be nontrivial and remain unprojected");
    rejects([&]{hydro.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Zero,0.);});
    rejects([&]{hydro.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Smooth,.5);});
    rejects([&]{hydro.initializeState(mesh,geometry);});
    equalState(hydro.state(),reference.state());
    require(!hydro.state().progress.stopped,"repeated initialization changed failure state");
}

void pressureMetricsContract() {
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    auto& cell=mesh.cells(mesh.findCellLabel({.31,.35}));
    const double pi=std::acos(-1.0);
    const auto near=[](double actual,double expected) {
        require(std::isfinite(actual) && std::abs(actual-expected)<1e-13,
            "native pressure geometry changed face coefficients or flux rate");
    };
    const auto xy=HydroNativePressureDetail::metrics<false>(mesh,cell);
    near(xy.west(),.8);near(xy.east(),.8);near(xy.axial(),1.25);
    near(xy.fluxRate(.2,2.,1.,3.,1.),1.75);
    const auto rz=HydroNativePressureDetail::metrics<true>(mesh,cell);
    near(rz.west(),.4*pi);near(rz.east(),.6*pi);near(rz.axial(),.78125*pi);
    near(rz.fluxRate(.2,2.,1.,3.,1.),1.28125*pi);
    auto& axisCell=mesh.cells(mesh.findCellLabel({.01,.35}));
    const auto axis=HydroNativePressureDetail::metrics<true>(mesh,axisCell);
    require(axis.west()==0.,"axis must have zero west face area");
    near(axis.east(),.2*pi);near(axis.axial(),.15625*pi);
}

void momentumMetricsContract() {
    using namespace HydroNativeMomentumDetail;
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=512;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    DictionaryMesh::RefineBox box;box.min={.375,.3};box.max={.625,.7};box.level=1;mesh.applyRefineBoxes({box});
    const double pi=std::acos(-1.0);int covered[2][3]={};
    const auto near=[](double actual,double expected) {
        require(std::isfinite(actual) && std::abs(actual-expected)<1e-12,
            "native momentum metric differs from rectangular/annular geometry");
    };
    for(bool radial:{false,true})for(int i=0;i<(radial?mesh.vEdgeIndexNum():mesh.hEdgeIndexNum());++i) {
        const auto& edge=radial?mesh.vEdges(mesh.vEdgeIndexToLabel[i]):mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if(edge.masterCell<0 || edge.slaveCell<0)continue;
        const auto& master=mesh.cells(edge.masterCell);const auto& slave=mesh.cells(edge.slaveCell);
        const int relation=master.level()==slave.level()?0:master.level()<slave.level()?1:2;
        ++covered[radial][relation];const auto& cell=relation==1?slave:master;
        const auto h=cell.size();const double r=edge.midpointCoord().x;
        const double cap=pi*((r+.5*h.x)*(r+.5*h.x)-(r-.5*h.x)*(r-.5*h.x));
        const double east=2*pi*(r+.5*h.x)*h.y,west=2*pi*(r-.5*h.x)*h.y;
        near(volume<false>(edge,cell),h.x*h.y);near(volume<true>(edge,cell),cap*h.y);
        const auto verifyFlux=[&](const auto& xy,const auto& rz) {
            near(xy.dx,h.x);near(xy.dy,h.y);near(rz.dx,h.x);near(rz.dy,h.y);
            near(xy.westArea(),h.y);near(xy.eastArea(),h.y);near(xy.axialArea(),h.x);
            near(rz.westArea(),west);near(rz.eastArea(),east);near(rz.axialArea(),cap);
            near(xy.diffusionNet(2.,1.,3.,1.),h.y+2*h.x);
            near(rz.diffusionNet(2.,1.,3.,1.),2*east-west+2*cap);
            near(xy.diffusionNet(.7,.7,.7,.7),0.);
            near(rz.diffusionNet(.7,.7,.7,.7),.7*(east-west));
        };
        if(radial) {
            const auto verify=[&](const auto& xy,const auto& rz) {
                near(xy.shear(),h.x/h.y);near(xy.east(),h.y/h.x);near(xy.west(),h.y/h.x);
                near(xy.cross(-.7),-.7);near(xy.normalViscosity(.11,.17),.28*h.y/h.x);
                near(rz.shear(),cap/h.y);near(rz.east(),east/h.x);near(rz.west(),west/h.x);
                near(rz.cross(-.7),-.7*cap/h.x);near(rz.normalViscosity(.11,.17),(.11*west+.17*east)/h.x);
            };
            if(relation==0)verify(radialMetrics<false,FaceRefinement::SameLevel>(edge,master,slave),radialMetrics<true,FaceRefinement::SameLevel>(edge,master,slave));
            else if(relation==1)verify(radialMetrics<false,FaceRefinement::MasterCoarse>(edge,master,slave),radialMetrics<true,FaceRefinement::MasterCoarse>(edge,master,slave));
            else verify(radialMetrics<false,FaceRefinement::SlaveCoarse>(edge,master,slave),radialMetrics<true,FaceRefinement::SlaveCoarse>(edge,master,slave));
            if(relation==0)verifyFlux(radialFluxMetrics<false,FaceRefinement::SameLevel>(edge,master,slave),radialFluxMetrics<true,FaceRefinement::SameLevel>(edge,master,slave));
            else if(relation==1)verifyFlux(radialFluxMetrics<false,FaceRefinement::MasterCoarse>(edge,master,slave),radialFluxMetrics<true,FaceRefinement::MasterCoarse>(edge,master,slave));
            else verifyFlux(radialFluxMetrics<false,FaceRefinement::SlaveCoarse>(edge,master,slave),radialFluxMetrics<true,FaceRefinement::SlaveCoarse>(edge,master,slave));
        } else {
            const auto& start=mesh.verts(edge.edgeToVert[0]);const auto& end=mesh.verts(edge.edgeToVert[1]);
            const auto xy=axialMetrics<false>(edge,cell,start,end);const auto rz=axialMetrics<true>(edge,cell,start,end);
            verifyFlux(axialFluxMetrics<false>(edge,cell,start,end),axialFluxMetrics<true>(edge,cell,start,end));
            near(xy.shearWest(),h.y/h.x);near(xy.shearEast(),h.y/h.x);near(xy.normal(),h.x/h.y);
            near(xy.crossWest(-.7),-.7);near(xy.crossEast(.3),.3);
            near(rz.shearWest(),west/h.x);near(rz.shearEast(),east/h.x);near(rz.normal(),cap/h.y);
            near(rz.crossWest(-.7),-.7*west/h.y);near(rz.crossEast(.3),.3*east/h.y);
        }
    }
    for(const auto& orientation:covered)for(int count:orientation)require(count>0,"missing native momentum face refinement relation");
}

void nativeProjectionContract(bool axisymmetric,bool refined,bool varying,const std::filesystem::path& output,
    bool inviscid=false) {
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=512;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    if(refined) {
        DictionaryMesh::RefineBox box;box.min={.375,.3};box.max={.625,.7};box.level=1;
        mesh.applyRefineBoxes({box});require(mesh.cellIndexNum()>80,"projection test did not refine");
    }
    VEdgeField<double> vx(mesh,0.);HEdgeField<double> vy(mesh,0.);
    const double pi=std::acos(-1.0);
    for(int i=0;i<vx.size();++i) {
        const auto& edge=mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        if(edge.masterCell<0 || edge.slaveCell<0)continue;
        const auto p=edge.midpointCoord();
        vx.value[i]=.05*std::sin(pi*p.x)*std::sin(2*pi*p.y);
    }
    for(int i=0;i<vy.size();++i) {
        const auto& edge=mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        if(edge.masterCell<0 || edge.slaveCell<0)continue;
        const auto p=edge.midpointCoord();
        vy.value[i]=-.03*std::sin(2*pi*p.x)*std::sin(pi*p.y);
    }
    Hydrodynamics hydro(mesh);const auto geometry=HydroStepGeometry::complete(axisymmetric);
    hydro.initializeState(mesh,geometry,{nullptr,&vx,&vy});
    CellField<double> density(mesh,0.),viscosity(mesh,0.);
    VertField<double> vertexDensity(mesh,0.),vertexViscosity(mesh,0.);
    if(varying) {
        for(int i=0;i<mesh.cellIndexNum();++i) {
            const auto p=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
            density.value[i]=1.2+.4*p.x+.3*p.y;viscosity.value[i]=inviscid?0.:.04+.01*p.x+.005*p.y;
        }
        for(int i=0;i<mesh.vertIndexNum();++i) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
            vertexDensity.value[i]=1.2+.4*p.x+.3*p.y;vertexViscosity.value[i]=inviscid?0.:.04+.01*p.x+.005*p.y;
        }
        hydro.updateMaterials(mesh,{geometry,0.,HydroMaterialEvolution::Frozen,
            HydroNativeMaterial{&density,&viscosity,&vertexDensity,&vertexViscosity}});
    } else hydro.updateMaterials(mesh,{geometry,0.,HydroMaterialEvolution::Frozen,HydroConstantMaterial{1.2,inviscid?0.:.04}});
    const std::string name=std::string(axisymmetric?"rz":"xy")+(refined?"_refined":"_uniform")
        +(varying?"_variable":"")+(inviscid?"_inviscid":"");
    std::ofstream table;
    if(!output.empty()) {
        std::filesystem::create_directories(output);table.open(output/(name+".csv"));
        require(bool(table),"cannot open projection reference");
        table<<std::scientific<<std::setprecision(17)<<"step,field,index,label,value\n";
    }
    HydroSolveOptions options;
    for(int step=1;step<=2;++step) {
        const auto oldX=hydro.state().velocityX.copy();const auto oldY=hydro.state().velocityY.copy();
        const auto report=hydro.advance(mesh,{geometry,{step,(step-1)*.001,.001},options,
            {nullptr,nullptr,nullptr,nullptr,(step-1)*.001}});
        if(!report.success)throw std::runtime_error(name+": "+report.stage+": "+report.reason);
        require(report.momentumUsedSuperLU && report.pressureUsedSuperLU && report.maxSpeed>1e-4
            && report.momentumRelativeResidual<1e-8 && report.pressureRelativeResidual<1e-8
            && report.divergenceAfterProjection<1e-8,"nonzero native projection failed");
        std::cout<<std::setprecision(17)<<"Native balance: "<<name<<" step="<<step
            <<" coupled="<<report.coupledCorrectionUsed<<" linear="<<report.momentumRelativeResidual
            <<" full="<<report.momentumBalanceRelativeResidual
            <<" forceDensity="<<report.momentumBalanceMaxForceDensity<<'\n';
        require(report.coupledCorrectionUsed==!refined
            && std::isfinite(report.momentumBalanceRelativeResidual)
            && std::isfinite(report.momentumBalanceMaxForceDensity),"missing or misclassified full momentum residual");
        // With mu=0 the mass-weighted projection exactly closes the predictor
        // equation, also across virtual-child faces and with variable density.
        // Viscous splitting is reported separately; linear accuracy alone is
        // not evidence that the final velocity satisfies full CN momentum.
        if(inviscid)require(report.momentumBalanceRelativeResidual<1e-10
            && report.momentumBalanceMaxForceDensity<1e-8,"inviscid projection has a momentum splitting defect");
        else if(refined)require(report.momentumBalanceRelativeResidual>10*options.linearResidualTolerance
            && report.momentumBalanceMaxForceDensity>1e-4,"viscous splitting defect was hidden by a linear residual");
        const auto state=hydro.state();
        for(int i=0;i<oldX.size();++i)require(state.previousVelocityX[i]==oldX.value[i],"projection lost x history");
        for(int i=0;i<oldY.size();++i)require(state.previousVelocityY[i]==oldY.value[i],"projection lost y history");
        for(int i=0;i<vx.size();++i) {
            const auto& edge=mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
            if(edge.masterCell<0 || edge.slaveCell<0)require(state.velocityX[i]==0.,"projection changed closed x boundary");
        }
        for(int i=0;i<vy.size();++i) {
            const auto& edge=mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
            if(edge.masterCell<0 || edge.slaveCell<0)require(state.velocityY[i]==0.,"projection changed closed y boundary");
        }
        if(table.is_open()) {
            const auto field=[&](const char* key,const auto& values) {
                for(int i=0;i<values.size();++i)table<<step<<','<<key<<','<<i<<','<<values.label(i)<<','<<values[i]<<'\n';
            };
            field("pressure",state.physicalPressure);field("ux",state.velocityX);field("uy",state.velocityY);
            field("ux_previous",state.previousVelocityX);field("uy_previous",state.previousVelocityY);
            field("ux_predictor",state.predictorX);field("uy_predictor",state.predictorY);
            field("divergence",state.divergence);field("speed",state.speed);
        }
    }
    if(table.is_open()) {table.flush();require(bool(table),"cannot write projection reference");}
    std::cout<<"Native projection: "<<name<<"; nonzero two-step velocity, continuity, boundaries and histories passed\n";
}

void boundaryProjectionAlgebra() {
    using Projection=HydroBoundaryProjection;
    static_assert(std::is_same_v<decltype(std::declval<const Projection&>().equations()),const HydroPressureSystem&>,
        "projection equations must be read-only");
    for(double scale:{1.,.01}) {
        std::vector<Projection::Face> faces{{0,-1,2.,0.,-1.,.5,true},
            {1,0,2.,.25*scale,0.,0.,false},{1,-1,2.,.5*scale,1.,0.,false}};
        const std::vector<double> predicted{.5,.75,.25};
        const Projection projection(2,faces,[&](std::size_t i){return predicted[i];});
        faces[0].area=999.;faces[1].high=99;
        const auto& equation=projection.equations();
        require(projection.face(0).area==2. && projection.face(1).high==1,"projection borrowed mutable face data");
        require(equation.matrix.coeff(0,0)==scale && equation.matrix.coeff(0,1)==-scale
            && equation.matrix.coeff(1,0)==-scale && equation.matrix.coeff(1,1)==3.*scale
            && equation.rhs[0]==-.5 && equation.rhs[1]==1.,"two-cell boundary Schur equation differs from analytic reference");
        Eigen::VectorXd multiplier(2);multiplier<<-.25/scale,.25/scale;
        std::vector<double> corrected=predicted;corrected[0]=-99.;
        const auto read=[&](std::size_t i){return corrected[i];};
        const auto write=[&](std::size_t i,double value){corrected[i]=value;};
        projection.apply(multiplier,read,write,true,"nonfinite test projection");
        require(corrected[0]==-99. && corrected[1]==.5 && corrected[2]==.5,
            "in-place initial correction must preserve fixed faces and apply inverse mass");
        corrected=predicted;corrected[0]=-99.;
        projection.apply(multiplier,read,write,false,"nonfinite test projection");
        require(corrected==std::vector<double>({.5,.5,.5}),"step correction must impose the fixed boundary constraint");
        const auto flux=projection.fluxes(read);
        require(flux.integrated.isZero(0.) && flux.absolute[0]==2. && flux.absolute[1]==2.,
            "boundary and interior incidence must close each physical cell flux");
        const auto saved=corrected;Eigen::VectorXd wrong(1);wrong.setZero();
        rejects([&]{projection.apply(wrong,read,write,false,"nonfinite test projection");});
        require(corrected==saved,"invalid multiplier size wrote a velocity");
        multiplier[1]=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{projection.apply(multiplier,read,write,true,"nonfinite test projection");});
        require(corrected==saved,"nonfinite first free-face correction was published");
        rejects([&]{Projection invalid(2,faces,[&](std::size_t i){return predicted[i];});});
    }
}

void boundaryProjectionContract(bool axisymmetric,bool varying,const std::filesystem::path& output) {
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);Hydrodynamics hydro(mesh);
    const auto geometry=HydroStepGeometry::complete(axisymmetric);
    CellField<double> pressure(mesh,.03),density(mesh,1.2),viscosity(mesh,.04);
    VertField<double> vertexDensity(mesh,1.2),vertexViscosity(mesh,.04);
    VEdgeField<double> vx(mesh,0.);HEdgeField<double> vy(mesh,0.);
    if(varying) {
        for(int i=0;i<density.size();++i) {
            const auto p=mesh.cells(mesh.cellIndexToLabel[i]).centerCoord();
            density.value[i]+=.4*p.x+.3*p.y;
        }
        for(int i=0;i<vertexDensity.size();++i) {
            const auto p=mesh.verts(mesh.vertIndexToLabel[i]).coord();
            vertexDensity.value[i]+=.4*p.x+.3*p.y;
        }
    }
    hydro.initializeState(mesh,geometry,{&pressure,&vx,&vy});
    hydro.updateMaterials(mesh,{geometry,0.,HydroMaterialEvolution::Frozen,
        HydroNativeMaterial{&density,&viscosity,&vertexDensity,&vertexViscosity}});
    HydroSolveOptions options;
    options.boundaryPatches={{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0.,1.,-.05},
        {"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0.,1.,0.}};
    const auto initialized=hydro.initializeCompatibleVelocity(mesh,options,axisymmetric);
    require(initialized.success && initialized.usedSuperLU && initialized.maxSpeed>.04
        && initialized.divergenceBeforeProjection>0. && initialized.divergenceAfterProjection<1e-10,
        "actual zero-time boundary projection failed");
    auto state=hydro.state();
    require(state.progress.completedSteps==0 && state.progress.acceptedSteps==0 && state.progress.endTime==0.,
        "zero-time projection consumed a physical step");
    for(int i=0;i<pressure.size();++i)require(state.physicalPressure[i]==.03,"initial multiplier became physical pressure");
    equalField(state.velocityX,state.previousVelocityX);equalField(state.velocityX,state.predictorX);
    equalField(state.velocityY,state.previousVelocityY);equalField(state.velocityY,state.predictorY);
    const auto savedX=state.velocityX.copy();const auto savedY=state.velocityY.copy();
    const auto savedDiv=state.divergence.copy();const auto savedSpeed=state.speed.copy();
    auto invalid=options;invalid.boundaryPatches[0].value=.05;
    const auto rejected=hydro.initializeCompatibleVelocity(mesh,invalid,axisymmetric);
    require(!rejected.success && rejected.usedSuperLU && rejected.stage=="projection-validation"
        && rejected.reason.find("backflow")!=std::string::npos,"backflow must be rejected after the actual initial solve");
    state=hydro.state();
    for(int i=0;i<savedX.size();++i)require(state.velocityX[i]==savedX.value[i]
        && state.previousVelocityX[i]==savedX.value[i] && state.predictorX[i]==savedX.value[i],"failed initialization changed x layers");
    for(int i=0;i<savedY.size();++i)require(state.velocityY[i]==savedY.value[i]
        && state.previousVelocityY[i]==savedY.value[i] && state.predictorY[i]==savedY.value[i],"failed initialization changed y layers");
    for(int i=0;i<pressure.size();++i)require(state.physicalPressure[i]==.03
        && state.divergence[i]==savedDiv.value[i] && state.speed[i]==savedSpeed.value[i],"failed initialization published diagnostics/pressure");
    require(!state.progress.stopped && state.progress.completedSteps==0,"initial rejection locked the NS history");
    const std::string name=std::string("boundary_")+(axisymmetric?"rz":"xy")+(varying?"_variable":"_constant");
    std::ofstream table;
    if(!output.empty()) {
        std::filesystem::create_directories(output);table.open(output/(name+".csv"));
        require(bool(table),"cannot open boundary projection reference");
        table<<std::scientific<<std::setprecision(17)<<"step,field,index,label,value\n";
    }
    for(int step=0;step<=2;++step) {
        if(step) {
            const auto report=hydro.advance(mesh,{geometry,{step,(step-1)*.001,.001},options,
                {nullptr,nullptr,nullptr,nullptr,(step-1)*.001}});
            if(!report.success)throw std::runtime_error(name+": "+report.stage+": "+report.reason);
            require(report.pressureUsedSuperLU && report.coupledUsedSuperLU && report.maxSpeed>.04
                && report.divergenceAfterProjection<1e-10,"boundary step projection/coupled closure failed");
        }
        state=hydro.state();
        if(table.is_open()) {
            const auto field=[&](const char* key,const auto& values) {
                for(int i=0;i<values.size();++i)table<<step<<','<<key<<','<<i<<','<<values.label(i)<<','<<values[i]<<'\n';
            };
            field("pressure",state.physicalPressure);field("ux",state.velocityX);field("uy",state.velocityY);
            field("ux_previous",state.previousVelocityX);field("uy_previous",state.previousVelocityY);
            field("ux_predictor",state.predictorX);field("uy_predictor",state.predictorY);
            field("divergence",state.divergence);field("speed",state.speed);
        }
    }
    if(table.is_open()){table.flush();require(bool(table),"cannot write boundary projection reference");}
    std::cout<<"Boundary projection: "<<name<<"; initial state, numerical rejection and two steps passed\n";
}

void fittedBoundaryProtection() {
    DictionaryMesh input;input.Nx=input.Ny=16;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;
    const auto active=[](const QuadtreeMesh&,const Cell& c){const auto p=c.centerCoord();return p.x<.5 || p.y>.5;};
    graph.build(mesh,active);
    std::vector<HydroFittedDomain::Boundary> boundaries(graph.segments().size());
    for(std::size_t i=0;i<boundaries.size();++i) {
        const auto& face=graph.segments()[i];
        if(face.L>=0 && active(mesh,mesh.cells(face.L)))continue;
        boundaries[i].type=face.faceCenter.x==0. ? HydroFittedDomain::BoundaryType::Symmetry
            : face.faceCenter.y==1. ? HydroFittedDomain::BoundaryType::PressureOutlet
            : HydroFittedDomain::BoundaryType::NoSlipWall;
    }
    const auto domain=HydroFittedDomain::build(mesh,graph,boundaries,true);
    const auto indexAt=[&](Vector2D p){return mesh.cells(mesh.findCellLabel(p)).index();};
    std::vector<std::uint8_t> mask(mesh.cellIndexNum(),0);
    mask[indexAt({.25,.75})]=1;
    protectHydroFittedBoundarySupport(mesh,domain,mask);
    require(mask[indexAt({.46875,.28125})] && mask[indexAt({.53125,.28125})],
        "fitted wall protection must retain both fluid and solid neighbors");
    require(!mask[indexAt({.75,.25})] && mask[indexAt({.25,.75})],
        "boundary protection froze deep solid or lost another module's mark");
    const auto saved=mask;
    DictionaryMesh otherInput=input;otherInput.Nx=8;QuadtreeMesh other(otherInput);
    rejects([&]{protectHydroFittedBoundarySupport(other,domain,mask);});
    require(mask==saved,"foreign mesh modified support");
    mask.pop_back();rejects([&]{protectHydroFittedBoundarySupport(mesh,domain,mask);});
}

void candidateContract(bool axisymmetric, bool fitted, bool refined=false) {
    DictionaryMesh input;input.Nx=input.Ny=8;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    if(refined) {
        DictionaryMesh::RefineBox box;box.min={.375,.375};box.max={.625,.625};box.level=1;
        mesh.applyRefineBoxes({box});require(mesh.cellIndexNum()>64,"candidate test did not refine");
    }
    FaceSegmentGraph graph;std::optional<HydroFittedDomain> domain;
    if(fitted) {
        graph.build(mesh);std::vector<HydroFittedDomain::Boundary> boundary(graph.segments().size());
        for(std::size_t i=0;i<boundary.size();++i)if(graph.segments()[i].L<0) {
            const auto side=graph.segments()[i].side;
            if(side==FaceSide::South || side==FaceSide::North) {
                boundary[i].type=HydroFittedDomain::BoundaryType::PressureReservoir;
                boundary[i].pressure=side==FaceSide::South?.2:0.;
                boundary[i].incomingDensity=1.;boundary[i].incomingViscosity=.01;
            } else boundary[i].type=side==FaceSide::West?HydroFittedDomain::BoundaryType::Symmetry
                :HydroFittedDomain::BoundaryType::NoSlipWall;
        }
        domain.emplace(HydroFittedDomain::build(mesh,graph,boundary,true));
        std::vector<std::uint8_t> mask(mesh.cellIndexNum(),0);
        protectHydroFittedBoundarySupport(mesh,*domain,mask);
        for(int i=0;i<mesh.cellIndexNum();++i) {
            const auto& c=mesh.cells(mesh.cellIndexToLabel[i]);
            const auto p=c.centerCoord(),h=c.size();
            const bool boundaryCell=p.x==h.x/2 || p.y==h.y/2 || p.x==1-h.x/2 || p.y==1-h.y/2;
            require(bool(mask[i])==boundaryCell,"physical boundary protection missed a wall/port or froze an interior cell");
        }
    }
    const auto geometry=fitted?HydroStepGeometry::fittedRZ(*domain):HydroStepGeometry::complete(axisymmetric);
    Hydrodynamics hydro(mesh),reference(mesh);
    HydroSolveOptions options;options.accelerationY=.03;
    const auto material=[&](Hydrodynamics& h,double time) {
        h.updateMaterials(mesh,{geometry,time,HydroMaterialEvolution::TimeDependent,HydroConstantMaterial{1.,.01}});
    };
    for(auto* h:{&hydro,&reference}) {h->initializeState(mesh,geometry);material(*h,0.);}
    if(fitted) {
        auto wrong=options;wrong.completeConservation.emplace();
        const auto rejected=hydro.advance(mesh,{geometry,{1,0.,.001},wrong,{}});
        require(!rejected.success && rejected.stage=="step-input" && !hydro.state().progress.stopped,
            "fitted path silently consumed complete-domain conservation settings");
        Hydrodynamics uninitialized(mesh);
        rejects([&]{HydroFittedMeshTransfer transfer(mesh,uninitialized);});
        rejects([&]{HydroFittedMeshTransfer transfer(mesh,hydro);});
        equalState(hydro.state(),reference.state());
    }
    for(int step=1;step<=2;++step) {
        HydroStepInput request{geometry,{step,(step-1)*.001,.001},options,{nullptr,nullptr,nullptr,nullptr,(step-1)*.001}};
        {auto wrong=request;wrong.time.policy=HydroTimePolicy::ConservativeOneStep;
            const auto rejected=hydro.advance(mesh,wrong);
            require(!rejected.success && !hydro.state().progress.stopped && hydro.state().progress.completedSteps==step-1,
                "ordinary/fitted AB2 Hydro accepted conservative variable-time policy");}
        const auto restored=[&] {
            equalState(hydro.state(),reference.state());
            const auto a=hydro.state(),b=reference.state();
            require(a.progress.completedSteps==b.progress.completedSteps && a.progress.acceptedSteps==b.progress.acceptedSteps
                && a.progress.endTime==b.progress.endTime && a.progress.acceptedEndTime==b.progress.acceptedEndTime
                && !a.progress.stopped && !a.progress.pendingAcceptance && a.materials && b.materials
                && a.materials->revision==b.materials->revision && a.materials->sourceTime==b.materials->sourceTime,
                "transaction did not restore physical progress/material provenance");
        };
        {
            auto transaction=hydro.beginStepTransaction(mesh);
            rejects([&]{hydro.beginStepTransaction(mesh);});rejects([&]{transaction.commit();});
            if(fitted) rejects([&]{HydroFittedMeshTransfer transfer(mesh,hydro);});
            auto trial=hydro.prepareStep(mesh,request);
            require(trial.report().success && trial.pending(),"transaction fixture did not perform the physical solve");
            hydro.updateMaterials(mesh,{geometry,step*.001,HydroMaterialEvolution::TimeDependent,HydroConstantMaterial{3.,.2}});
            require(hydro.state().density[0]==3.,"trial target material did not change");
            transaction.rollback();require(!transaction.pending() && !trial.pending(),"rollback retained pending handles");
            rejects([&]{trial.accept(mesh);});rejects([&]{transaction.commit();});
        }
        restored();
        {
            auto transaction=hydro.beginStepTransaction(mesh);
            auto trial=hydro.prepareStep(mesh,request);
            require(trial.report().success,"rollback on scope exit did not reach a real candidate");
            rejects([&]{hydro.updateMaterials(mesh,{geometry,step*.001,HydroMaterialEvolution::TimeDependent,HydroConstantMaterial{-1.,.2}});});
            // Trial abandonment first stops the local solve; the enclosing
            // transaction must then restore it, including all history fields.
        }
        restored();
        {
            auto transaction=hydro.beginStepTransaction(mesh);
            auto failing=options;failing.linearResidualTolerance=-1.;
            const auto failed=hydro.prepareStep(mesh,{geometry,request.time,failing,request.force});
            require(!failed.report().success && hydro.state().progress.stopped,"kernel failure fixture did not stop");
        }
        restored();
        auto transaction=hydro.beginStepTransaction(mesh);
        auto candidate=hydro.prepareStep(mesh,request);
        const auto direct=reference.advance(mesh,request);
        if(!candidate.report().success)throw std::runtime_error(candidate.report().reason);
        require(direct.success && candidate.pending() && !candidate.accepted(),"missing actual candidate solve");
        require(candidate.report().momentumUsedSuperLU && candidate.report().pressureUsedSuperLU
            && candidate.report().coupledUsedSuperLU,"candidate bypassed SuperLU");
        equalState(hydro.state(),reference.state());
        auto p=hydro.state().progress;
        require(p.pendingAcceptance && !p.stopped && p.completedSteps==step && p.acceptedSteps==step-1,
            "local solve was accepted before coupling");
        require(p.acceptedEndTime==(step-1)*.001 && p.endTime==step*.001,"candidate time identity changed");
        // Immediate acceptance does not refresh time-dependent coefficients.
        // AMR must reject that stale source before a legitimate material update.
        if(fitted) {
            rejects([&]{HydroFittedMeshTransfer transfer(mesh,reference);});
            equalState(hydro.state(),reference.state());
            require(!reference.state().progress.stopped,"rejected capture poisoned accepted Hydro state");
        }
        material(hydro,step*.001);material(reference,step*.001);
        HydroStepInput next{geometry,{step+1,step*.001,.001},options,{nullptr,nullptr,nullptr,nullptr,step*.001}};
        require(!hydro.advance(mesh,next).success,"advance bypassed pending acceptance");
        { auto duplicate=hydro.prepareStep(mesh,next);require(!duplicate.report().success && !duplicate.pending(),"duplicate candidate accepted"); }
        require(candidate.pending() && !hydro.state().progress.stopped,"rejected entry invalidated the real candidate");
        rejects([&]{hydro.writeCache(mesh,"unused/");});rejects([&]{hydro.readCache(mesh,"unused/");});
        if(fitted)rejects([&]{HydroFittedMeshTransfer transfer(mesh,hydro);});
        QuadtreeMesh foreign(input);rejects([&]{candidate.accept(foreign);});
        auto moved=std::move(candidate);
        rejects([&]{candidate.accept(mesh);});rejects([&]{candidate.report();});
        require(!candidate.pending() && moved.pending(),"candidate move lost ownership");
        moved.accept(mesh);require(moved.accepted() && !moved.pending(),"candidate did not accept");
        transaction.commit();require(!transaction.pending(),"successful physical transaction was not released");
        rejects([&]{moved.accept(mesh);});
        p=hydro.state().progress;
        require(!p.pendingAcceptance && !p.stopped && p.acceptedSteps==step && p.acceptedEndTime==step*.001,
            "accepted progress did not match the numerical substep");
        equalState(hydro.state(),reference.state());
        if(fitted) {HydroFittedMeshTransfer transfer(mesh,hydro);require(transfer.pressure.size()==mesh.cellIndexNum(),"accepted capture failed");}
    }
    if(fitted) {
        HydroFittedMeshTransfer transfer(mesh,hydro);
        for(auto* f:{&transfer.velocityR,&transfer.previousVelocityR})
            for(int i=0;i<f->size();++i)f->value[i]=0.;
        for(auto* f:{&transfer.velocityZ,&transfer.previousVelocityZ})
            for(int i=0;i<f->size();++i)f->value[i]=0.;
        const HydroFittedNumerics numerics(1e-10);
        auto compatible=transfer.prepare(mesh,*domain,numerics);
        require(!compatible.report.usedSuperLU && compatible.report.current.projectionSolves==0
            && compatible.report.previous.projectionSolves==0,"compatible history performed a pressure solve");
        int pulse=-1;
        for(int i=0;i<mesh.vEdgeIndexNum();++i) {
            const auto& edge=mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
            if(edge.masterCell>=0 && edge.slaveCell>=0){pulse=i;break;}
        }
        require(pulse>=0,"no interior transfer probe face");
        transfer.velocityR.value[pulse]=.1;transfer.previousVelocityR.value[pulse]=-.03;
        auto projected=transfer.prepare(mesh,*domain,numerics);
        require(projected.report.usedSuperLU && projected.report.current.projectionSolves>0
            && projected.report.previous.projectionSolves>0,"both history levels must require projection");
        require(projected.report.current.localContinuityResidual<=numerics.localContinuity
            && projected.report.previous.localContinuityResidual<=numerics.localContinuity,
            "transferred histories failed original continuity gate");
        require(transfer.velocityR.value[pulse]==.1 && transfer.previousVelocityR.value[pulse]==-.03,
            "history projection changed staging fields");
        const auto state=projected.solver->state();
        require(state.progress.acceptedSteps==2 && state.progress.endTime==.002,
            "algebraic history projection advanced physical time");
        equalState(hydro.state(),reference.state());
        std::cout<<"Transfer projection: "<<(refined?"refined":"uniform")<<"; compatible solves=0; current="
            <<projected.report.current.projectionSolves<<" previous="<<projected.report.previous.projectionSolves<<'\n';
    }
    std::unique_ptr<HydroFittedVelocitySampler> fittedSample;
    std::unique_ptr<StaggeredVelocitySampler> completeSample;
    Vector2D saved;
    {
        auto abandoned=hydro.prepareStep(mesh,{geometry,{3,.002,.001},options,{nullptr,nullptr,nullptr,nullptr,.002}});
        require(abandoned.report().success && abandoned.pending(),"abandon test did not solve");
        if(fitted) {
            fittedSample=std::make_unique<HydroFittedVelocitySampler>(hydro.fittedVelocity(mesh,*domain,
                HydroVelocityTimeLevel::Current,HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible));
            saved=fittedSample->sample({.5,.5});
        } else {completeSample=std::make_unique<StaggeredVelocitySampler>(hydro.averagedNativeVelocity(mesh));saved=completeSample->sample({.5,.5});}
    }
    const auto stopped=hydro.state().progress;
    require(stopped.stopped && !stopped.pendingAcceptance && stopped.completedSteps==3 && stopped.acceptedSteps==2,
        "abandonment failed to stop an unaccepted local step");
    const auto retry=hydro.advance(mesh,{geometry,{4,.003,.001},options,{nullptr,nullptr,nullptr,nullptr,.003}});
    require(!retry.success && retry.requiresReinitialization,"abandoned Hydro allowed continuation");
    rejects([&]{material(hydro,.003);});
    if(fitted) {
        rejects([&]{hydro.fittedVelocity(mesh,*domain,HydroVelocityTimeLevel::Current,
            HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible);});
        rejects([&]{HydroFittedMeshTransfer transfer(mesh,hydro);});
        equal(fittedSample->sample({.5,.5}),saved);
    } else {rejects([&]{hydro.averagedNativeVelocity(mesh);});equal(completeSample->sample({.5,.5}),saved);}
    std::cout<<"Hydro candidate: "<<(fitted?(refined?"native fitted RZ":"fitted RZ"):(axisymmetric?"complete RZ":"complete XY"))
        <<"; unchanged fields/history, deferred acceptance and fail-stop passed\n";
}

void candidateLifetime() {
    static_assert(!std::is_copy_constructible_v<HydroStepCandidate> && std::is_nothrow_move_constructible_v<HydroStepCandidate>);
    DictionaryMesh input;input.Nx=input.Ny=8;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);const auto geometry=HydroStepGeometry::complete(false);HydroSolveOptions options;
    const auto make=[&] {
        auto h=std::make_unique<Hydrodynamics>(mesh);h->initializeState(mesh,geometry);
        h->updateMaterials(mesh,{geometry,0.,HydroMaterialEvolution::Frozen,HydroConstantMaterial{1.,.01}});return h;
    };
    const HydroStepInput request{geometry,{1,0.,.001},options,{}};
    auto a=make(),b=make();auto ca=a->prepareStep(mesh,request),cb=b->prepareStep(mesh,request);
    require(ca.pending() && cb.pending(),"lifetime test did not solve");
    ca=std::move(cb);
    require(a->state().progress.stopped && b->state().progress.pendingAcceptance && !cb.pending(),"move assignment did not abandon old owner");
    ca.accept(mesh);require(b->state().progress.acceptedSteps==1,"moved candidate accepted wrong owner");
    auto dying=make();auto orphan=dying->prepareStep(mesh,request);dying.reset();
    require(!orphan.pending() && !orphan.accepted() && orphan.report().success,"owner destruction invalidated retained report");
    rejects([&]{orphan.accept(mesh);});
    auto fresh=make();
    {auto bad=fresh->prepareStep(mesh,{geometry,{2,0.,.001},options,{}});require(!bad.report().success && !bad.pending(),"bad sequence accepted");}
    require(!fresh->state().progress.stopped,"preflight rejection stopped valid state");
    {auto stale=fresh->prepareStep(mesh,request);auto& v=mesh.verts(mesh.vertIndexToLabel[0]);const auto point=v.coord();
     v.setCoord(point+Vector2D{.001,0.});rejects([&]{stale.accept(mesh);});v.setCoord(point);}
    require(fresh->state().progress.stopped && fresh->state().progress.acceptedSteps==0,"invalidated geometry published candidate");
    auto failure=make();auto failing=options;failing.linearResidualTolerance=-1.;
    auto rejected=failure->prepareStep(mesh,{geometry,{1,0.,.001},failing,{}});
    require(!rejected.report().success && rejected.report().requiresReinitialization && !rejected.pending()
        && failure->state().progress.stopped,"numerical failure lost fail-stop semantics");
    rejects([&]{rejected.accept(mesh);});
    std::cout<<"Hydro candidate lifetime: move assignment, owner destruction, stale mesh and failure checks passed\n";
    static_assert(!std::is_copy_constructible_v<HydroStepTransaction> && std::is_nothrow_move_constructible_v<HydroStepTransaction>);
    auto left=make(),right=make();auto first=left->beginStepTransaction(mesh),second=right->beginStepTransaction(mesh);
    auto trial=left->prepareStep(mesh,request);require(trial.pending(),"transaction move fixture");
    first=std::move(second);
    require(!left->state().progress.stopped && left->state().progress.completedSteps==0 && !trial.pending()
        && first.pending() && !second.pending(),"transaction move assignment failed to restore old owner");
    first.rollback();
    auto temporary=make();auto detached=temporary->beginStepTransaction(mesh);temporary.reset();
    require(!detached.pending(),"destroyed owner left live transaction");rejects([&]{detached.rollback();});
    auto edited=make();
    {auto guard=edited->beginStepTransaction(mesh);auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);
        const auto p=vertex.coord();vertex.setCoord(p+Vector2D{.001,0});rejects([&]{guard.rollback();});vertex.setCoord(p);}
    require(!edited->state().progress.stopped,"restored mesh prevented rollback");
    auto modified=make();const auto origin=mesh.verts(mesh.vertIndexToLabel[0]).coord();
    {auto guard=modified->beginStepTransaction(mesh);mesh.verts(mesh.vertIndexToLabel[0]).setCoord(origin+Vector2D{.001,0});}
    mesh.verts(mesh.vertIndexToLabel[0]).setCoord(origin);
    require(modified->state().progress.stopped,"changed mesh silently received an old transaction snapshot");
    std::cout<<"Hydro transaction: physical rollback/retry, material history, ownership and changed-mesh isolation passed\n";
}

void conservationContract()
{
    DictionaryMesh input;input.Nx=8;input.Ny=10;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);const auto geometry=HydroStepGeometry::complete(false);
    const auto solve=[&](const HydroSolveOptions& options) {
        Hydrodynamics hydro(mesh);
        hydro.initializeCompleteVelocity(mesh,geometry,HydroInitialVelocity::Smooth,.05);
        hydro.updateMaterials(mesh,{geometry,0.,HydroMaterialEvolution::Frozen,HydroConstantMaterial{1.2,.04}});
        return hydro.advance(mesh,{geometry,{1,0.,.001},options,{}});
    };
    HydroSolveOptions options;options.linearResidualTolerance=2e-5;
    const auto baseline=solve(options);
    require(baseline.success && baseline.momentumBalanceRelativeResidual>0.,"conservation fixture requires an actual nonzero floating residual");
    options.completeConservation.emplace();
    options.completeConservation->momentumBalance=.5*baseline.momentumBalanceRelativeResidual;
    const auto physical=solve(options);
    require(!physical.success && physical.stage=="coupled-balance-check" && physical.requiresReinitialization
        && physical.momentumBalanceRelativeResidual==baseline.momentumBalanceRelativeResidual
        && physical.momentumRelativeResidual<options.linearResidualTolerance,
        "physical momentum gate did not act independently of successful linear solves");
    options.completeConservation->momentumBalance=1e-8;options.linearResidualTolerance=1e-30;
    const auto linear=solve(options);
    require(!linear.success && !linear.coupledCorrectionUsed && linear.reason.find("relative linear residual")!=std::string::npos,
        "physical tolerance replaced linear solve accuracy");
    options.linearResidualTolerance=2e-5;options.completeConservation->momentumBalance=-1.;
    const auto invalid=solve(options);
    require(!invalid.success && invalid.stateChange==HydroStateChange::Unchanged,"invalid conservation control changed physical fields");
    std::cout<<"Complete Hydro conservation: legacy expansion, independent physical/linear gates and invalid controls passed\n";
}

void velocityLatticeContract()
{
    DictionaryMesh input;input.Nx=128;input.Ny=168;input.maxCellNum=43008;input.maxInitialIterations=0;
    input.minCoord={0,0};input.maxCoord={4,5.2359877559829888};
    QuadtreeMesh mesh(input);
    VEdgeField<double> vx(mesh,0.);HEdgeField<double> vy(mesh,0.);
    const auto exact=[](Vector2D p){return Vector2D{.3+.2*p.x-.1*p.y,-.2+.4*p.x+.3*p.y};};
    for(int i=0;i<vx.size();++i)vx.value[i]=exact(mesh.vEdges(mesh.vEdgeIndexToLabel[i]).midpointCoord()).x;
    for(int i=0;i<vy.size();++i)vy.value[i]=exact(mesh.hEdges(mesh.hEdgeIndexToLabel[i]).midpointCoord()).y;
    const StaggeredVelocitySampler sample(mesh,vx,vy);
    double maxError=0;
    for(Vector2D p:{Vector2D{.1,.1},Vector2D{1.23,3.45},Vector2D{3.91,5.19}}) {
        const auto expected=exact(p),difference=sample.sample(p)-expected;
        const double error=std::hypot(difference.x,difference.y);
        maxError=std::max(maxError,error);
        require(error<=32*std::numeric_limits<double>::epsilon()*std::max(1.,std::hypot(expected.x,expected.y)),
            "native interpolation did not reproduce an affine velocity at rounding accuracy");
    }
    auto& edge=mesh.vEdges(mesh.vEdgeIndexToLabel[0]);const auto p=edge.midpointCoord();
    edge.setMidpointCoord(p+Vector2D{0,.001*mesh.cellSize.y});
    bool rejected=false;
    try {const StaggeredVelocitySampler invalid(mesh,vx,vy);}
    catch(const std::exception& error){rejected=std::string(error.what()).find("nonuniform native face locations")!=std::string::npos;}
    edge.setMidpointCoord(p);
    require(rejected,"off-lattice face was accepted or rejected for an unrelated reason");
    std::cout<<"Native velocity lattice: 128x168 non-dyadic domain affine field and off-lattice rejection passed; max error="<<maxError<<'\n';
}
void splineAccuracyContract()
{
    double previousValue=0,previousGradient=0;
    for(int n:{16,32}) {
        DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;
        input.maxCellNum=4*(n+1)*(n+1);input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);
        VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);const double pi=std::acos(-1.);
        const auto psi=[&](Vector2D p){return std::sin(pi*p.x)*std::sin(pi*p.y)/pi;};
        for(const auto face:domain.faces(mesh)) {
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            (face.horizontal?uy.value:ux.value)[edge.index()]=face.outwardSign?0:(face.horizontal?-1:1)*(psi(b)-psi(a))/face.area;
        }
        HydroSplineVelocity velocity(mesh,domain,ux,uy);double valueError=0,gradientError=0;
        // Two-point Gauss nodes coincide with leading spline second-derivative
        // error zeros; use four-point integration to measure the full error.
        constexpr double nodes[4]={-.8611363115940525752,-.3399810435848562648,.3399810435848562648,.8611363115940525752};
        constexpr double weights[4]={.3478548451374538574,.6521451548625461426,.6521451548625461426,.3478548451374538574};
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int a=0;a<4;++a)for(int b=0;b<4;++b) {
            const Vector2D p{(i+.5*(1+nodes[a]))/n,(j+.5*(1+nodes[b]))/n};const auto jet=velocity.firstJet(p);
            const double sx=std::sin(pi*p.x),cx=std::cos(pi*p.x),sy=std::sin(pi*p.y),cy=std::cos(pi*p.y);
            const auto v=jet.value-Vector2D{sx*cy,-cx*sy},dx=jet.dx-Vector2D{pi*cx*cy,pi*sx*sy},dy=jet.dy-Vector2D{-pi*sx*sy,-pi*cx*cy};
            const double weight=.25*weights[a]*weights[b];
            valueError+=weight*(v.x*v.x+v.y*v.y);gradientError+=weight*(dx.x*dx.x+dx.y*dx.y+dy.x*dy.x+dy.y*dy.y);
        }
        valueError=std::sqrt(valueError/(n*n));gradientError=std::sqrt(gradientError/(n*n));
        std::cout<<std::setprecision(17)<<"HYDRO_SPLINE_ACCURACY cells="<<n<<" velocity_rms="<<valueError<<" gradient_rms="<<gradientError
            <<" velocity_ratio="<<(previousValue?previousValue/valueError:0)<<" gradient_ratio="<<(previousGradient?previousGradient/gradientError:0)<<std::endl;
        if(previousValue)require(previousValue/valueError>6 && previousGradient/gradientError>3,"compatible smooth-field spline accuracy did not improve under refinement");
        previousValue=valueError;previousGradient=gradientError;
    }
}

void splineSolvedLayersContract(bool timeLinear=false)
{
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=8;
    input.maxCellNum=1024;input.maxInitialIterations=0;QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto native=NativeLeafDomain::build(mesh,graph,false);const auto domain=HydroStepGeometry::complete(false);
    auto hydro=std::make_unique<Hydrodynamics>(mesh);
    rejects([&]{hydro->splineVelocity(mesh,native,HydroVelocityTimeLevel::Current);});
    hydro->initializeCompleteVelocity(mesh,domain,HydroInitialVelocity::Smooth,.05);
    hydro->updateMaterials(mesh,{domain,0,HydroMaterialEvolution::Frozen,HydroConstantMaterial{1,.01}});
    rejects([&]{hydro->splineVelocity(mesh,native,HydroVelocityTimeLevel::Previous);});
    auto retained=hydro->splineVelocity(mesh,native,HydroVelocityTimeLevel::Current);const auto initial=retained.firstJet({.37,.61});
    auto state=InterfaceParametricMotion::initialize(mesh,native,InterfacePrescribedData::parametricEllipse({.5,.75},{.15,.15},16));
    InterfaceParametricProjection::Options geometryOptions;HydroSolveOptions solveOptions;double fluxError=0,layerChange=0,maximumFraction=0,maximumBalance=0;
    for(int k=1;k<=2;++k) {
        HydroForceDensity force;force.evaluationTime=(k-1)*.001;
        const auto report=hydro->advance(mesh,{domain,{k,force.evaluationTime,.001},solveOptions,force});
        require(report.success && report.momentumUsedSuperLU && report.pressureUsedSuperLU,"native C1 layer control did not complete its NS solve");
        const auto fields=hydro->state();
        for(auto layer:{HydroVelocityTimeLevel::Previous,HydroVelocityTimeLevel::Current}) {
            const auto snapshot=hydro->splineVelocity(mesh,native,layer);const bool old=layer==HydroVelocityTimeLevel::Previous;
            for(const auto face:native.faces(mesh)) {
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);const int i=edge.index();
                const double now=face.horizontal?fields.velocityY[i]:fields.velocityX[i],before=face.horizontal?fields.previousVelocityY[i]:fields.previousVelocityX[i];
                layerChange=std::max(layerChange,std::abs(now-before));const auto& f=snapshot.faceFlux(face.nativeRow);
                fluxError=std::max(fluxError,std::abs(static_cast<double>(f.integrate(0,0,1)+f.integrate(1,0,1))-face.area*(old?before:now)));
            }
        }
        const auto snapshot=hydro->splineVelocity(mesh,native,HydroVelocityTimeLevel::Current);
        const auto previous=hydro->splineVelocity(mesh,native,HydroVelocityTimeLevel::Previous);
        const auto flow=timeLinear?HydroInterfaceTransport::timeLinearParametricFlow(previous,snapshot,.001):HydroInterfaceTransport::frozenParametricFlow(snapshot);
        const auto flux=HydroInterfaceTransport::normalFlux(snapshot),beforeFlux=HydroInterfaceTransport::normalFlux(previous);
        const auto oldInventory=state.parametricGeometry().inventory();
        const InterfaceStepTime time{k,(k-1)*.001,k*.001,.001};
        auto next=timeLinear?InterfaceParametricMotion::prepareTimeLinear(mesh,native,state,time,flow,beforeFlux,flux,geometryOptions)
            :InterfaceParametricMotion::prepare(mesh,native,state,time,.001,flow,flux,geometryOptions);
        require(state.parametricGeometry().inventory()==oldInventory && state.progress().completedSteps==k-1,"time-layer transport modified its old interface state");
        maximumFraction=std::max(maximumFraction,next.state.parametricGeometry().report().maximumInventoryFractionError);
        maximumBalance=std::max(maximumBalance,std::abs(next.inventoryBalance));state=std::move(next.state);
    }
    hydro.reset();const auto after=retained.firstJet({.37,.61});
    require(initial.value.x==after.value.x && initial.value.y==after.value.y && layerChange>0 && fluxError<2e-16,
        "C1 snapshots lost native time-layer identity or owning lifetime");
    require(maximumFraction<=3e-12 && maximumBalance<2e-15,"solved native layers violated interface inventory budget");
    std::cout<<std::setprecision(17)<<(timeLinear?"HYDRO_SPLINE_TIME_SOLVED_LAYERS":"HYDRO_SPLINE_SOLVED_LAYERS")<<" ns_steps=2 interface_steps="<<state.progress().completedSteps
        <<" physical_time="<<state.progress().endTime<<" native_face_flux_error="<<fluxError<<" layer_change="<<layerChange
        <<" maximum_fraction_error="<<maximumFraction<<" inventory_balance="<<maximumBalance
        <<" retained_snapshot_unchanged=1 temporal_order_claim=0 capillary_coupling=0\n";
}

void splineVelocityContract()
{
    constexpr int n=8;DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;
    input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);
    const UniformLeafGrid grid(mesh);
    for(bool residual:{false,true}) {
        VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);
        const auto stream=[](Vector2D p){return 4*p.x*p.x*std::pow(1-p.x,2)*p.y*p.y*std::pow(1-p.y,2);};
        bool perturbed=false;
        for(const auto face:domain.faces(mesh)) {
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            double speed=(face.horizontal?-1:1)*(stream(b)-stream(a))/face.area;
            if(residual && !perturbed && !face.horizontal && face.outwardSign==0 && face.center.x>.3 && face.center.x<.7 && face.center.y>.3 && face.center.y<.7){speed+=.003;perturbed=true;}
            (face.horizontal?uy.value:ux.value)[edge.index()]=speed;
        }
        HydroSplineVelocity velocity(mesh,domain,ux,uy);
        const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);
        const auto flux=HydroInterfaceTransport::normalFlux(velocity);
        double fluxError=0,densityError=0,divergence=0,commutingError=0,jetJump=0,halfJetJump=0,mapDerivativeError=0;
        const auto jetDifference=[](const auto& a,const auto& b){return (a.dx-b.dx).mod()+(a.dy-b.dy).mod();};
        for(const auto face:domain.faces(mesh)) {
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);const auto chart=flux.faceChart(face.nativeRow);
            const double expected=face.area*(face.horizontal?uy.value[edge.index()]:ux.value[edge.index()]);
            const double integrated=static_cast<double>(flux.integratedNormalFlux(face.nativeRow,0,0,1)+flux.integratedNormalFlux(face.nativeRow,1,0,1));
            fluxError=std::max(fluxError,std::abs(integrated-expected));
            for(double s:{0.,.173,.5,.827,1.}) {
                const auto p=chart.begin+(chart.end-chart.begin)*s;const auto j=velocity.firstJet(p);const int half=s<=.5?0:1;
                const double density=static_cast<double>(flux.normalFluxDensity(face.nativeRow,half,false,2*s-half)/(face.area*.5));
                densityError=std::max(densityError,std::abs(density-(face.horizontal?j.value.y:j.value.x)));
                divergence=std::max(divergence,std::abs(j.dx.x+j.dy.y));
            }
            if(face.outwardSign==0) {
                const auto p=chart.begin+(chart.end-chart.begin)*.371;const Vector2D offset=face.horizontal?Vector2D{0,1e-7}:Vector2D{1e-7,0};
                jetJump=std::max(jetJump,jetDifference(velocity.firstJet(p-offset),velocity.firstJet(p+offset)));
                halfJetJump=std::max(halfJetJump,jetDifference(velocity.firstJet(p-offset*.5),velocity.firstJet(p+offset*.5)));
            }
        }
        constexpr double gauss=.28867513459481288225;
        for(int j=0;j<n;++j)for(int i=0;i<n;++i) {
            const auto& cell=mesh.cells(grid.cells()[j*n+i]);const auto lo=mesh.verts(cell.cellToVert[0]).coord();
            double integral=0;for(double a:{.5-gauss,.5+gauss})for(double b:{.5-gauss,.5+gauss}) {
                const auto jet=velocity.firstJet(lo+Vector2D{a/n,b/n});const double div=jet.dx.x+jet.dy.y;
                integral+=div/(4*n*n);divergence=std::max(divergence,std::abs(div));
            }
            const double native=(ux.value[mesh.vEdges(grid.vertical()[j*(n+1)+i+1]).index()]-ux.value[mesh.vEdges(grid.vertical()[j*(n+1)+i]).index()]
                +uy.value[mesh.hEdges(grid.horizontal()[(j+1)*n+i]).index()]-uy.value[mesh.hEdges(grid.horizontal()[j*n+i]).index()])/n;
            commutingError=std::max(commutingError,std::abs(integral-native));
        }
        int tangentQueries=0;
        for(Vector2D p:{Vector2D{.25,.5},Vector2D{.4375,.5},Vector2D{.5,.25},Vector2D{.5,.4375}})for(double dt:{-.04,.04}) {
            const auto map=flow(p,dt);constexpr double h=1e-6;
            const auto derivative=[&](Vector2D offset){return (flow(p-offset*2,dt).point-flow(p-offset,dt).point*8+
                flow(p+offset,dt).point*8-flow(p+offset*2,dt).point)*(1/(12*h));};
            mapDerivativeError=std::max({mapDerivativeError,(derivative({h,0})-map.dx).mod(),(derivative({0,h})-map.dy).mod()});++tangentQueries;
        }
        require(fluxError<2e-16 && densityError<2e-13 && commutingError<2e-14,"C1 velocity lost native flux or commuting divergence");
        require(residual?divergence>.01:divergence<2e-12,"C1 velocity changed its discrete divergence source");
        require(halfJetJump<=.51*jetJump+1e-13 && mapDerivativeError<2e-7,"C1 velocity or trajectory derivative is not continuous/consistent");
        std::cout<<std::setprecision(17)<<"HYDRO_SPLINE_VELOCITY residual="<<residual<<" face_flux_error="<<fluxError<<" face_density_error="<<densityError
            <<" maximum_divergence="<<divergence<<" commuting_divergence_error="<<commutingError<<" seam_jet_difference="<<jetJump
            <<" half_offset_difference="<<halfJetJump<<" tangent_queries="<<tangentQueries<<" map_derivative_fd_error="<<mapDerivativeError<<std::endl;
        if(residual)continue;
        auto state=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.5,.75},{.15,.15},16));
        InterfaceParametricProjection::Options options;double maximumBalance=0,maximumFractionError=0;
        for(int k=1;k<=2;++k) {
            const auto inventory=state.parametricGeometry().inventory();
            auto next=InterfaceParametricMotion::prepare(mesh,domain,state,{k,(k-1)*.04,k*.04,.04},k==1?.04:-.04,flow,flux,options);
            require(state.progress().completedSteps==k-1 && state.parametricGeometry().inventory()==inventory,"C1 velocity update changed the old interface");
            maximumBalance=std::max(maximumBalance,std::abs(next.inventoryBalance));
            maximumFractionError=std::max(maximumFractionError,next.state.parametricGeometry().report().maximumInventoryFractionError);state=std::move(next.state);
        }
        require(maximumBalance<2e-15 && maximumFractionError<=3e-12,"C1 native velocity exceeded the existing inventory budget");
        std::cout<<"HYDRO_SPLINE_STATE steps="<<state.progress().completedSteps<<" physical_time="<<state.progress().endTime
            <<" inventory_balance="<<maximumBalance<<" maximum_fraction_error="<<maximumFractionError<<" old_state_unchanged=1\n";
    }
    splineAccuracyContract();splineSolvedLayersContract();
}

void parametricBridgeContract()
{
    constexpr int n=8;constexpr double dt=.04;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;
    input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);
    const auto domain=NativeLeafDomain::build(mesh,graph,false);
    VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);
    const auto stream=[](Vector2D p){return 4*p.x*p.x*std::pow(1-p.x,2)*p.y*p.y*std::pow(1-p.y,2);};
    for(const auto face:domain.faces(mesh)) {
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double speed=(face.horizontal?-1:1)*(stream(b)-stream(a))/face.area;
        (face.horizontal?uy.value:ux.value)[edge.index()]=speed;
    }
    const auto velocity=HydroInterfaceTransport::completeVelocity(mesh,domain,ux,uy,{});
    const auto flow=HydroInterfaceTransport::frozenParametricFlow(mesh,domain,velocity);
    const auto flux=HydroInterfaceTransport::normalFlux(velocity);
    // Independent exhaustive ownership on this dyadic root grid. Compare the
    // entire map jet through native half-cell seams, in both time directions.
    std::vector<LevelSetGaussTrace::SmoothPatch> referencePatches;
    for(int index:domain.cellIndexByRow()) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);
        const auto lo=mesh.verts(cell.cellToVert[0]).coord();
        for(int i=0;i<2;++i)for(int j=0;j<2;++j)
            referencePatches.push_back({lo+Vector2D{double(i)/(2*n),double(j)/(2*n)},
                lo+Vector2D{double(i+1)/(2*n),double(j+1)/(2*n)}});
    }
    const auto patchAt=[&](Vector2D p,Vector2D direction) {
        for(const auto& patch:referencePatches) {
            if(!LevelSetGaussTrace::contains(patch,p))continue;
            if((p.x==patch.lower.x && direction.x<0) || (p.x==patch.upper.x && direction.x>0) ||
               (p.y==patch.lower.y && direction.y<0) || (p.y==patch.upper.y && direction.y>0))continue;
            return patch;
        }
        throw std::runtime_error("reference point left the native domain");
    };
    const auto query=[&](Vector2D p){const auto j=velocity.volumeFirstJet(p);
        return LevelSetGaussTrace::VelocityJet{j.value,j.dx,j.dy};};
    int comparisons=0,crossings=0,unresolved=0;
    for(double x:{.25,.25+1e-9,.4375,.5,.75-1e-9})for(double y:{.25,.25+1e-9,.4375,.5,.75-1e-9})
        for(double duration:{-dt,dt}) {
            std::optional<LevelSetGaussTrace::Map> reference;
            std::optional<InterfaceParametricTransport::MapJet> produced;
            std::string referenceError,productionError;
            try{reference=LevelSetGaussTrace::advance(query,patchAt,{x,y},duration);}
            catch(const std::runtime_error& error){referenceError=error.what();}
            try{produced=flow({x,y},duration);}
            catch(const std::runtime_error& error){productionError=error.what();}
            if(!referenceError.empty() || !productionError.empty()) {
                require(referenceError==productionError,"native adapter changed a seam rejection");
                const auto jet=query({x,y});const auto acceleration=jet.dx*jet.value.x+jet.dy*jet.value.y;
                std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_UNRESOLVED x="<<x<<" y="<<y
                    <<" duration="<<duration<<" velocity="<<jet.value.x<<","<<jet.value.y
                    <<" convective_acceleration="<<acceleration.x<<","<<acceleration.y<<" reason="<<referenceError<<'\n';
                ++unresolved;continue;
            }
            const auto& expected=*reference;const auto& actual=*produced;
            require(actual.point.x==expected.point.x && actual.point.y==expected.point.y &&
                actual.dx.x==expected.derivative(0,0) && actual.dx.y==expected.derivative(1,0) &&
                actual.dy.x==expected.derivative(0,1) && actual.dy.y==expected.derivative(1,1) &&
                actual.timeDerivative.x==expected.timeDerivative.x && actual.timeDerivative.y==expected.timeDerivative.y &&
                actual.displacement.x==expected.displacement.x && actual.displacement.y==expected.displacement.y &&
                actual.coordinateResolution==expected.coordinateResolution,"production native flow changed the reference first jet");
            ++comparisons;crossings+=expected.events;
        }
    require(crossings>0,"native flow comparison did not cross a smooth-patch seam");
    double fluxError=0;
    for(const auto face:domain.faces(mesh)) {
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const double native=face.area*(face.horizontal?uy.value[edge.index()]:ux.value[edge.index()]);
        const double integral=static_cast<double>(flux.integratedNormalFlux(face.nativeRow,0,0,1)+flux.integratedNormalFlux(face.nativeRow,1,0,1));
        fluxError=std::max(fluxError,std::abs(integral-native));
    }
    require(fluxError<2e-16,"parametric bridge changed native face fluxes");
    rejects([&]{HydroInterfaceTransport::frozenParametricFlow(mesh,domain,velocity,0);});
    const auto rz=NativeLeafDomain::build(mesh,graph,true);
    rejects([&]{HydroInterfaceTransport::frozenParametricFlow(mesh,rz,velocity);});
    QuadtreeMesh other(input);
    DictionaryMesh::RefineBox box;box.min={.25,.25};box.max={.5,.5};box.level=1;other.applyRefineBoxes({box});
    rejects([&]{HydroInterfaceTransport::frozenParametricFlow(other,domain,velocity);});
    auto state=InterfaceParametricMotion::initialize(mesh,domain,
        InterfacePrescribedData::parametricEllipse({.5,.75},{.15,.15},16));
    InterfaceParametricProjection::Options options;
    double maximumBalance=0,maximumFractionError=0;
    for(int k=1;k<=2;++k) {
        const auto inventory=state.parametricGeometry().inventory();
        const auto* geometry=&state.parametricGeometry().geometry();
        auto next=InterfaceParametricMotion::prepare(mesh,domain,state,{k,(k-1)*dt,k*dt,dt},
            k==1?dt:-dt,flow,flux,options);
        require(state.progress().completedSteps==k-1 && state.parametricGeometry().inventory()==inventory &&
            &state.parametricGeometry().geometry()==geometry,"native preparation changed the previous state");
        require(next.state.progress().completedSteps==k && next.state.progress().endTime==k*dt &&
            next.state.parametricGeometry().inventory()!=inventory,"native parametric update did not advance its state");
        maximumBalance=std::max(maximumBalance,std::abs(next.inventoryBalance));
        maximumFractionError=std::max(maximumFractionError,next.state.parametricGeometry().report().maximumInventoryFractionError);
        state=std::move(next.state);
    }
    require(maximumBalance<2e-15 && maximumFractionError<=3e-12,"native parametric state exceeded the existing inventory budget");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_BRIDGE comparisons="<<comparisons<<" seam_crossings="<<crossings<<" matching_rejections="<<unresolved
        <<" map_jet_identical=1 native_flux_error="<<fluxError<<" steps="<<state.progress().completedSteps
        <<" physical_time="<<state.progress().endTime<<" inventory_balance="<<maximumBalance
        <<" maximum_fraction_error="<<maximumFractionError<<" old_state_unchanged=1\n";
}
}

void hydroMidpointTimeContract();
int main(int argc,char** argv)
try {
    if(argc==2 && std::string(argv[1])=="--velocity-lattice"){velocityLatticeContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-bridge"){parametricBridgeContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--spline-velocity"){splineVelocityContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--spline-time"){hydroSplineTimeContract(true);splineSolvedLayersContract(true);return 0;}
    if(argc==2 && std::string(argv[1])=="--spline-time-components"){hydroSplineTimeContract(false);splineSolvedLayersContract(true);return 0;}
    if(argc==2 && std::string(argv[1])=="--compatible-pressure"){hydroCompatiblePressureContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--compatible-moments"){hydroCompatibleMomentsContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-viscosity"){hydroMomentViscosityContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-viscosity"){hydroMomentRZViscosityContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-transport"){hydroMomentRZTransportContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-reconstruction"){hydroMomentRZReconstructionContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-time"){hydroMomentRZTimeContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-manufactured"){hydroMomentRZManufacturedContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-coordinates"){hydroMomentRZCoordinatesContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-time"){hydroMomentTimeContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-ports"){hydroMomentPortsContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-port-spline"){hydroMomentPortSplineContract();return 0;}
        if(argc==2 && std::string(argv[1])=="--parametric-supply"){hydroParametricSupplyContract();return 0;}
        if(argc==2 && std::string(argv[1])=="--moment-supply"){hydroMomentSupplyContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-components"){hydroMomentComponentsContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-components-transfer"){hydroParametricComponentsTransferContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-neck"){hydroParametricNeckContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-neck"){hydroMomentNeckContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-neck-readiness"){hydroMomentNeckReadinessContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-transport"){hydroMomentTransportContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-ns-manufactured"){hydroMomentNSManufacturedContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-physical-input"){hydroMomentPhysicalInputContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-material"){hydroMomentMaterialContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-material-flow"){hydroMomentMaterialFlowContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-spline"){hydroMomentSplineContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-amr-spline"){hydroMomentAmrSplineContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-amr-two-phase"){hydroMomentAmrTwoPhaseContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-mesh-transfer"){hydroMomentMeshTransferContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-port-mesh-transfer"){hydroMomentPortMeshTransferContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-refinement"){hydroMomentNativeRefinementContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-material"){hydroMomentNativeMaterialContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-cartesian-gram"){hydroMomentCartesianMaterialGramContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-state"){hydroMomentNativeStateContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-physical-complement"){hydroMomentPhysicalComplementContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-block-solve"){hydroMomentBlockSolveContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-weighted-phase"){hydroParametricWeightedPhaseContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-flux"){hydroParametricWeightedPhaseContract();hydroMomentNativeFluxContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-mass-current"){hydroMomentMassCurrentContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-momentum"){hydroMomentNativeMomentumContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-coupled"){hydroMomentNativeCoupledContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-capillary-metric"){hydroMomentCapillaryMetricContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-remesh"){hydroMomentNativeRemeshContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-accuracy"){hydroMomentNativeAccuracyContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-native-density-transport"){hydroMomentNativeDensityTransportContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-local-resolution"){hydroMomentLocalResolutionContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-spline-interface"){hydroMomentSplineInterfaceContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-spline-flux-roots"){hydroMomentSplineFluxRootsContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-positive-flux"){hydroMomentPositiveFluxContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-two-phase"){hydroMomentTwoPhaseContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-capillary"){hydroMomentCapillaryContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-material-space"){hydroMomentMaterialSpaceContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-spline"){hydroMomentRZSplineContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-interface"){hydroMomentRZInterfaceContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-surface"){hydroMomentRZInterfaceContract(false);return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-two-phase"){hydroMomentRZTwoPhaseContract();return 0;}
    if(argc==3 && std::string(argv[1])=="--moment-rz-flow-replay"){hydroMomentRZFlowReplay(argv[2]);return 0;}
    if(argc>=3 && std::string(argv[1])=="--moment-phase-free-slip-probe"){
        for(int input=2;input<argc;++input)hydroMomentFreeSlipPhaseProbe(argv[input]);return 0;}
    if(argc>=3 && std::string(argv[1])=="--moment-phase-free-resolution"){
        for(int input=2;input<argc;++input)hydroMomentFreeSlipPhaseProbe(argv[input],"resolution");return 0;}
    if(argc>=3 && std::string(argv[1])=="--moment-phase-material-replay"){
        for(int input=2;input<argc;++input)hydroMomentFreeSlipPhaseProbe(argv[input],"material");return 0;}
    if(argc>=3 && std::string(argv[1])=="--moment-phase-initial-layers"){
        for(int input=2;input<argc;++input)hydroMomentFreeSlipPhaseProbe(argv[input],"initialLayers");return 0;}
    if(argc==2 && std::string(argv[1])=="--moment-rz-kinetic"){hydroMomentRZKineticContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--energy-viscosity"){hydroEnergyViscosityContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--midpoint-time"){hydroMidpointTimeContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-mass"){hydroParametricMassContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-measure"){hydroParametricMeasureContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-material-support"){hydroParametricMaterialSupportContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-curvature"){hydroParametricCurvatureContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-nonlinear-curvature"){hydroParametricNonlinearCurvatureContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-amr"){hydroParametricAmrContract();return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-capillary"){hydroParametricCapillaryContract(0);return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-capillary-motion"){hydroParametricCapillaryContract(24);return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-capillary-motion-half-step"){hydroParametricCapillaryContract(48);return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-midpoint"){hydroParametricMidpointContract(0);return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-midpoint-oscillation"){hydroParametricMidpointContract(200);return 0;}
    if(argc==2 && std::string(argv[1])=="--parametric-midpoint-oscillation-half-step"){hydroParametricMidpointContract(400);return 0;}
    if(argc>=3 && std::string(argv[1])=="--parametric-flux-replay"){
        for(int input=2;input<argc;++input)hydroParametricFluxReplay(argv[input]);return 0;}
    require(argc<=2,"expected optional native projection output directory");
    const std::filesystem::path output=argc==2?argv[1]:"";
    velocityLatticeContract();
    domainInputContract();
    runInputContract();
    conservationContract();
    interfaceInputContract();
    for(bool rz:{false,true}) for(bool refined:{false,true}) velocityInitializationContract(rz,refined);
    momentumMetricsContract();
    boundaryProjectionAlgebra();
    for(bool axisymmetric:{false,true})for(bool refined:{false,true})for(bool varying:{false,true})
        nativeProjectionContract(axisymmetric,refined,varying,output);
    for(bool axisymmetric:{false,true})for(bool varying:{false,true})
        nativeProjectionContract(axisymmetric,true,varying,{},true);
    for(bool axisymmetric:{false,true})for(bool varying:{false,true})
        boundaryProjectionContract(axisymmetric,varying,output);
    pressureMetricsContract();fittedBoundaryProtection();candidateContract(false,false);candidateContract(true,false);
    candidateContract(true,true);candidateContract(true,true,true);candidateLifetime();
    DictionaryMesh input;
    input.Nx = input.Ny = 8; input.maxCellNum = 256; input.maxInitialIterations = 0;
    QuadtreeMesh mesh(input);
    auto owner = std::make_unique<Hydrodynamics>(mesh);
    auto& hydro = *owner;
    const auto domain = HydroStepGeometry::complete(false);
    rejects([&] { hydro.vertexVelocity(mesh, HydroVelocityTimeLevel::Current, false); });
    VEdgeField<double> vx(mesh, 0.0);
    HEdgeField<double> vy(mesh, 0.0);
    const double pi = std::acos(-1.0);
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i) {
        const auto& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        const auto p = edge.midpointCoord();
        if (edge.boundary() == Boundary::non) vx.value[i] = .05*std::sin(pi*p.x)*std::sin(pi*p.y);
    }
    hydro.initializeState(mesh, domain, {nullptr, &vx, &vy});
    hydro.updateMaterials(mesh, {domain, 0., HydroMaterialEvolution::Frozen, HydroConstantMaterial{1., .01}});
    rejects([&] { hydro.averagedNativeVelocity(mesh); });
    rejects([&] { hydro.vertexVelocity(mesh, static_cast<HydroVelocityTimeLevel>(99), false); });
    const auto initial = hydro.vertexVelocity(mesh, HydroVelocityTimeLevel::Current, false);
    const auto initialReference = hydroVertexVelocity(mesh, vx, vy, false);
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
        require(initial.first->value[i] == initialReference.first->value[i]
            && initial.second->value[i] == initialReference.second->value[i], "initial vertex values changed");
    auto pressureCopy = hydro.state().physicalPressure.copy();
    auto materialCopy = hydro.state().vertexDensity.copy();
    pressureCopy.value[0] = 321.; materialCopy.value[0] = 123.;
    require(hydro.state().physicalPressure[0] == 0. && hydro.state().vertexDensity[0] == 1.,
        "output copy aliases Hydro storage");
    std::vector<double> lastX(vx.value.get(), vx.value.get() + vx.size());
    std::vector<double> lastY(vy.value.get(), vy.value.get() + vy.size());
    std::unique_ptr<StaggeredVelocitySampler> retained;
    std::function<Vector2D(Vector2D)> retainedBoundary;
    Vector2D retainedValue, retainedWallValue;
    HydroSolveOptions options;
    for (int step = 1; step <= 2; ++step) {
        HydroForceDensity force; force.evaluationTime = (step-1)*.001;
        const auto report = hydro.advance(mesh, {domain, {step, force.evaluationTime, .001}, options, force});
        if (!report.success || !report.momentumUsedSuperLU || !report.pressureUsedSuperLU ||
            !report.coupledUsedSuperLU || !std::isfinite(report.maxSpeed) || report.maxSpeed <= 0. ||
            report.momentumRelativeResidual > 1e-8 || report.pressureRelativeResidual > 1e-8 ||
            report.divergenceAfterProjection > 1e-8) {
            std::cerr << report.stage << ": " << report.reason << '\n'; return 1;
        }
        const auto state = hydro.state();
        require(state.progress.completedSteps == step && state.progress.endTime == step*.001,
            "state query lost its physical time layer");
        for (int i = 0; i < state.velocityX.size(); ++i) {
            require(state.previousVelocityX[i] == lastX[i], "previous x velocity is not the actual previous step");
            lastX[i] = state.velocityX[i];
        }
        for (int i = 0; i < state.velocityY.size(); ++i) {
            require(state.previousVelocityY[i] == lastY[i], "previous y velocity is not the actual previous step");
            lastY[i] = state.velocityY[i];
        }
        auto currentX = state.velocityX.copy(), previousX = state.previousVelocityX.copy();
        auto currentY = state.velocityY.copy(), previousY = state.previousVelocityY.copy();
        const StaggeredVelocitySampler reference(mesh, previousX, previousY, currentX, currentY);
        const auto averaged = hydro.averagedNativeVelocity(mesh);
        for (int j = 1; j < 8; ++j) for (int i = 1; i < 8; ++i)
            equal(averaged.sample({i/8., j/8.}), reference.sample({i/8., j/8.}));
        for (auto layer : {HydroVelocityTimeLevel::Current, HydroVelocityTimeLevel::Previous}) {
            const bool old = layer == HydroVelocityTimeLevel::Previous;
            const auto nodes = hydro.vertexVelocity(mesh, layer, false);
            const auto direct = hydroVertexVelocity(mesh, old ? previousX : currentX, old ? previousY : currentY, false);
            for (int i = 0; i < mesh.vertIndexNum(); ++i)
                require(nodes.first->value[i] == direct.first->value[i] && nodes.second->value[i] == direct.second->value[i],
                    "named native time layer changed during vertex conversion");
        }
        const auto boundary = hydro.boundaryAdvectionVelocity(mesh, HydroVelocityTimeLevel::Current, {});
        const auto boundaryReference = hydroBoundaryAdvectionVelocity(mesh, currentX, currentY, {});
        for (const auto point : {Vector2D{0., .25}, Vector2D{1., .75}, Vector2D{.5, .5}})
            equal(boundary(point), boundaryReference(point));
        if (step == 1) {
            retained = std::make_unique<StaggeredVelocitySampler>(hydro.averagedNativeVelocity(mesh));
            retainedBoundary = boundary;
            retainedValue = retained->sample({.5, .5});
            retainedWallValue = retainedBoundary({0., .25});
        } else {
            equal(retained->sample({.5, .5}), retainedValue);
            equal(retainedBoundary({0., .25}), retainedWallValue);
        }
    }
    const auto final = hydro.vertexVelocity(mesh, HydroVelocityTimeLevel::Current, false);
    const double savedNode = final.first->value[mesh.vertIndexNum()/2];
    // A foreign layout is rejected before any interpolation or field access.
    DictionaryMesh otherInput = input; otherInput.Nx = otherInput.Ny = 4;
    QuadtreeMesh other(otherInput);
    rejects([&] { hydro.vertexVelocity(other, HydroVelocityTimeLevel::Current, false); });
    // Numerical failure stops samplers; already owning results remain usable.
    HydroSolveOptions failing = options; failing.linearResidualTolerance = -1.;
    const auto failed = hydro.advance(mesh, {domain, {3, .002, .001}, failing, {nullptr, nullptr, nullptr, nullptr, .002}});
    require(!failed.success && hydro.state().progress.stopped, "failed kernel did not stop the Hydro service");
    rejects([&] { hydro.vertexVelocity(mesh, HydroVelocityTimeLevel::Current, false); });
    owner.reset();
    equal(retained->sample({.5, .5}), retainedValue);
    equal(retainedBoundary({0., .25}), retainedWallValue);
    require(final.first->value[mesh.vertIndexNum()/2] == savedNode, "owned nodes did not survive Hydro destruction");
    std::cout << "Standalone Hydro: two SuperLU steps, private fields, named histories and owning velocity queries passed\n";
    return 0;
} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
