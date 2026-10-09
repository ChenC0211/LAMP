#include "Electrospray/ElectrosprayApp.h"
#include "Runtime/RunControls.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
#include "Electrospray/Setup/ElectrosprayModuleInputs.h"
#include "Electrospray/Setup/ElectrosprayHydroDomain.h"
#include "Electrospray/IO/ElectrosprayHydroDiagnostics.h"
#include "Electrospray/IO/ElectrosprayForceDiagnostics.h"
#include "Electrospray/IO/ElectrosprayNumericsReport.h"
#include "Electrospray/IO/ElectrosprayVersionedInput.h"
#include "Electrospray/IO/ElectrosprayResolvedReport.h"
#include "Electrospray/IO/ElectrosprayRunSession.h"
#include "Electrospray/IO/ElectrosprayMeshDiagnostics.h"
#include "Electrospray/Setup/ElectrosprayCaseBindings.h"

#include <algorithm>
#include <optional>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Electrospray/Compatibility/DictionaryElectrospray.h"
#include "Electrospray/Compatibility/ElectrosprayPhysicalConfig.h"
#include "Electrospray/Compatibility/ElectrosprayInputAdapters.h"
#include "Electrospray/IO/ElectrosprayTecplotOutput.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Electrospray/IO/ElectrosprayElectricDiagnostics.h"
#include "Electrospray/ElectrosprayTimeLoop.h"
#include "Electrospray/Setup/ElectrosprayInterface.h"
#include "Electrospray/IO/ElectrosprayInitialInterfaceDiagnostics.h"
#include "Electrospray/Coupling/ElectrosprayCapillarity.h"
#include "Electrospray/Coupling/ElectrosprayForces.h"
#include "Electrospray/Coupling/ElectrosprayConducting.h"
#include "Electrospray/Coupling/ElectrosprayMeshAdaptation.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/OutputLayout.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Electricity/SurfaceChargeField.h"

namespace
{
bool needsRefinement(const QuadtreeMesh& mesh, const CellField<double>& target)
{
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        if (!std::isfinite(target.value[i]))
            throw std::runtime_error("non-finite geometry refinement target");
        if (std::lround(target.value[i]) != mesh.cells(mesh.cellIndexToLabel[i]).level())
            return true;
    }
    return false;
}

void initializeGeometry(QuadtreeMesh& mesh, ImmersedGeometry<QuadtreeMesh>& geometry,
    const DictionaryGeometry& definition, const DictionaryMesh& meshDefinition, int iterations, bool diagnosticFields)
{
    // No flow/electric state exists during initial refinement. Use the shared
    // target-level provider, not the old Hydro-gradient fallback in adapt().
    geometry.initializeLevelSets(mesh, 0.0, diagnosticFields);
    geometry.markPhases(mesh);
    // Initial velocity is zero; physical indicators become available only
    // after the first accepted NS step. Geometry/phase minima still apply now.
    auto initialDefinition = meshDefinition;
    initialDefinition.physicsRefinement.enabled = false;
    QuadtreeMeshRefiner refiner(initialDefinition);
    refiner.attachGeometry(&definition, &geometry);
    CellField<double> target(mesh, 0.0), indicator(mesh, 0.0);
    for (int iteration = 0; iteration < iterations; ++iteration)
    {
        refiner.buildTargetLevels(mesh, target, &indicator);
        if (!needsRefinement(mesh, target)) break;
        std::vector<VertField<double>> unusedLevelSets;
        std::vector<VertField<double>*> vertexFields;
        for (auto& levelSet : geometry.levelSets) vertexFields.push_back(&levelSet);
        std::vector<HEdgeField<double>*> horizontalFields;
        std::vector<VEdgeField<double>*> verticalFields;
        std::vector<CellField<double>*> cellFields{&target};
        mesh.adapt(unusedLevelSets, indicator, vertexFields, horizontalFields, verticalFields, cellFields, &target);
        geometry.initializeLevelSets(mesh, 0.0, diagnosticFields);
        geometry.markPhases(mesh);
        std::cout << "[电喷雾] 初始网格细化 " << iteration + 1
                  << "：单元数=" << mesh.cellIndexNum() << "，顶点数=" << mesh.vertIndexNum() << '\n';
    }
    refiner.buildTargetLevels(mesh, target, &indicator);
    if (needsRefinement(mesh, target))
        throw std::runtime_error("initialAdaptIterations did not reach the configured geometry refinement targets");
    // Keep the dictionary's exact initial curve; populate derivatives/history/cuts
    // without replacing it by a pseudo-time redistanced initial condition.
    geometry.refreshLevelSetDerivedData(mesh, 0.0, false);
    geometry.checkCutInfo(mesh, true);
    geometry.markPhases(mesh);
}

void writeGeometry(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>& geometry,
    const DictionaryGeometry& definition, const RunControls& run)
{
    std::vector<std::string> names;
    std::vector<VertField<double>> fields;
    for (std::size_t i = 0; i < geometry.levelSets.size(); ++i)
    {
        names.push_back("levelSet_" + definition.levelSetName[i]);
        fields.push_back(geometry.levelSets[i]);
    }
    VertField<double> phase(mesh, -1.0);
    std::vector<int> counts(geometry.phases.size(), 0);
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int material = mesh.verts(mesh.vertIndexToLabel[i]).phase();
        if (material < 0 || material >= static_cast<int>(counts.size()))
            throw std::runtime_error("unassigned initial phase at vertex " + std::to_string(i));
        phase.value[i] = material + 1; // Match the dictionary's 1-based phase IDs.
        ++counts[material];
        for (const auto& field : fields)
            if (!std::isfinite(field.value[i]))
                throw std::runtime_error("non-finite initial level set");
    }
    names.push_back("phase_id");
    fields.push_back(phase);
    for (std::size_t i = 0; i < counts.size(); ++i)
        std::cout << "[电喷雾] 初始相 " << i + 1 << " '" << definition.phaseName[i]
                  << "'：顶点数=" << counts[i] << '\n';
    if (run.writeFields && run.writeInitial && run.diagnosticFields)
    {
        OutputLayout::fields();
        mesh.outputFields("fields/" + run.outputPrefix, names, fields);
        const auto path = OutputLayout::fields() / (run.outputPrefix + ".dat");
        if (!std::filesystem::exists(path) || std::filesystem::file_size(path) == 0)
            throw std::runtime_error("failed to write initial spray geometry");
        std::cout << "[电喷雾] 初始几何已输出：" << path.string()
                  << "；此文件尚不含电场或速度解。\n";
    }
}
}

ElectrosprayApp::ElectrosprayApp(const RunConfig& cfg) : cfg_(cfg) {}

int ElectrosprayApp::run()
{
    std::string stage = "configuration";
    try
    {
        ElectrosprayWorkflow workflow;
        std::optional<ElectrosprayConfig> legacyRun;
        std::optional<ElectrosprayPhysicalConfig> legacyPhysical;
        ElectrosprayModuleInputs legacyModules;
        std::unique_ptr<ElectrosprayVersionedInput> versioned;
        std::unique_ptr<ElectrosprayRunSession> session;
        if (cfg_.schemaVersion == 2)
        {
            versioned = std::make_unique<ElectrosprayVersionedInput>(readElectrosprayVersionedInput());
            workflow.problem = versioned->problem;
            cfg_.units = versioned->units; cfg_.normalization = versioned->normalization;
            cfg_.writeSIOutputs = true;
        }
        else
        {
            DictionaryElectrospray input;
            if (!input.readDict()) throw std::runtime_error("Electrospray block is required in system/run.dict");
            legacyRun = input.config();
            workflow = adaptLegacyElectrosprayWorkflow(*legacyRun);
            legacyPhysical = ElectrosprayPhysicalConfig::read();
            legacyModules = adaptLegacyElectrosprayInputs(*legacyPhysical, *legacyRun, workflow);
        }
        const auto& inputs = versioned ? versioned->modules : legacyModules;
        const auto& run = inputs.controls;
        const auto& flowConfig = inputs.flow;
        const auto& electricConfig = inputs.electric;
        const auto& liquid = inputs.materials.phases.at(flowConfig.liquidPhase);
        const auto& gas = inputs.materials.phases.at(flowConfig.gasPhase);
        const double surfaceTensionSI = inputs.materials.surfaceTensionSI.at(std::minmax(flowConfig.liquidPhase, flowConfig.gasPhase));
        const double electricLinearTolerance = inputs.systems.at("ElectricPotential").relativeTolerance;
        (void)electrosprayProblemLabel(workflow.problem);
        if(workflow.problem!=ElectrosprayProblem::Conducting && electricConfig.interfaceDiscretization!="cellHarmonic")
            throw std::runtime_error("q=0 spray branches retain cellHarmonic; sharpHermiteP3 is an explicit coupled choice");
        if (!cfg_.units) throw std::runtime_error("physical spray requires a units context");
        if (workflow.problem == ElectrosprayProblem::Conducting)
            preflightElectrosprayConduction(inputs, workflow, *cfg_.units);
        if (workflow.restartLoop >= 0)
            throw std::runtime_error("legacy initialization caches are not compatible with the restored spray workflow");
        if (!workflow.initializeLevelSets)
            throw std::runtime_error("a fresh spray run requires dictionary level-set initialization");
        if (!(run.dt > 0.0) || !std::isfinite(run.dt) || !(run.endTime > 0.0)
            || !std::isfinite(run.endTime) || !(run.writeInterval > 0.0) || !std::isfinite(run.writeInterval))
            throw std::runtime_error("dt, finalTime and writeDt must be finite and positive");
        if (run.logEverySteps < 1) throw std::runtime_error("logEvery must be a positive integer");
        if (run.outputPrefix.empty() || run.outputPrefix.find_first_of("/\\:") != std::string::npos
            || run.outputPrefix == "." || run.outputPrefix == "..")
            throw std::runtime_error("outputName must be a nonempty file prefix, not a path");
        DictionaryGeometry definition;
        if (!versioned) writeElectrosprayNumericsReport(*legacyPhysical, *legacyRun, workflow, cfg_.solvers.get());
        definition.readDict(cfg_.units.get());
        if (!definition.hasDomain || definition.levelSetName.empty() || definition.phaseName.empty())
            throw std::runtime_error("spray requires v3 Domain, LevelSet and Phase definitions");
        DictionaryMesh meshDefinition;
        meshDefinition.readDict(cfg_.units.get(), true);
        meshDefinition.ensureLevelSetCount(static_cast<int>(definition.levelSetName.size()));
        if (meshDefinition.runtimeAdaptation.enabled
            && (!versioned || workflow.problem != ElectrosprayProblem::Conducting || electricConfig.model != "leakyDielectric"
                || electricConfig.interfaceDiscretization != "sharpHermiteP3"))
            throw std::runtime_error("运行时 AMR 当前要求新版字典、导电运动液面与 sharpHermiteP3 电场");
        const int iterations = workflow.initialAdaptIterations < 0
            ? meshDefinition.maxInitialIterations : workflow.initialAdaptIterations;
        if (iterations < 0) throw std::runtime_error("initial refinement iterations must be nonnegative");
        if (versioned)
        {
            validateElectrosprayCaseBindings(definition, meshDefinition, inputs.flow, inputs.interfaceConfig, inputs.electric);
            const char* electricModelName = electricConfig.model == "leakyDielectric" ? "漏电介质" : "理想介电";
            std::cout << "[电喷雾] 算例配置与几何编号检查通过；电学模型："
                      << electricModelName << "；线性求解器：SuperLU。\n";
            if (cfg_.checkCase)
            {
                std::cout << "[电喷雾] 只读预检查完成：未建网格、未求解、未写输出。"
                             "实际网格对齐、几何质量和数值稳定性仍由初始化及时间推进检查。\n";
                return 0;
            }
            stage = "startup";
            session = std::make_unique<ElectrosprayRunSession>(std::filesystem::current_path());
            if (!session->prepare(cfg_.startMode, ElectrosprayRunSession::terminalInput(), std::cin, std::cout))
                return 0;
        }
        if (versioned) writeElectrosprayResolvedReport(*versioned, &meshDefinition);
        else writeElectrosprayResolvedReport(*legacyPhysical, *legacyRun, *cfg_.units, cfg_.schemaVersion, nullptr, &meshDefinition);
        OutputLayout::fields();
        std::cout << "[电喷雾] 初始化字典几何：计算域(代码长度)=(" << definition.domainMin.x
                  << ',' << definition.domainMin.y << ")--(" << definition.domainMax.x << ','
                  << definition.domainMax.y << ")，水平集数=" << definition.levelSetName.size()
                  << "，相数=" << definition.phaseName.size() << "。\n";
        stage = "geometry-initialize";
        QuadtreeMesh mesh(meshDefinition);
        ImmersedGeometry<QuadtreeMesh> geometry(definition, mesh);
        initializeGeometry(mesh, geometry, definition, meshDefinition, iterations,
            run.writeFields && run.writeInitial && run.diagnosticFields);
        stage = "geometry-output";
        writeGeometry(mesh, geometry, definition, run);
        stage = "initial-physical-electrostatics";
        auto electric = initializeElectrosprayElectrostatics(mesh, geometry, definition, flowConfig, inputs.interfaceConfig, electricConfig, liquid, gas, electricLinearTolerance, *cfg_.units);
        reportElectrosprayElectricInitialization(mesh, electric, electricConfig.conductorPotentialSI, *cfg_.units);
        if (run.writeFields && run.writeInitial && run.diagnosticFields)
            writeElectrosprayElectrostatics(mesh, electric, *cfg_.units,
                (OutputLayout::fields() / run.outputPrefix).string(), true, run.detailedLog);
        std::cout << "[电喷雾] 压力边界：供液端=" << flowConfig.supply.pressureSI
                  << " Pa，出口=" << flowConfig.outlet.pressureSI << " Pa；供液速度由求解确定。\n";
        stage = "fitted-fluid-domain";
        const auto fluidDomain = buildElectrosprayHydroDomain(mesh, electric, definition, inputs.flow,
            inputs.materials.phases.at(inputs.flow.liquidPhase), inputs.materials.phases.at(inputs.flow.gasPhase),
            inputs.electric.conductorPhaseId, *cfg_.units);
        reportElectrosprayHydroDomain(fluidDomain);
        if (workflow.problem == ElectrosprayProblem::PressureDriven)
        {
            stage = "frozen-material-hydrodynamics";
            runElectrosprayFrozenMaterialHydro(mesh, fluidDomain, electric, inputs, workflow, *cfg_.units);
            return 0; // Only the explicitly selected pressure-only NS interval.
        }
        stage = "free-interface-initialize";
        auto freeInterface = initializeElectrosprayInterface(mesh, geometry, definition, fluidDomain, electric, flowConfig, inputs.interfaceConfig, electricConfig.conductorPhaseId);
        reportElectrosprayInterfaceInitialization(mesh, freeInterface, inputs.interfaceConfig);
        if (run.writeFields && run.writeInitial && run.diagnosticFields) writeElectrosprayInitialInterface(mesh, freeInterface, *cfg_.units,
            (OutputLayout::fields() / run.outputPrefix).string());
        stage = "active-interface-geometry";
        const auto initialGeometryOutput = makeElectrosprayInitialGeometryObserver(mesh, freeInterface, *cfg_.units,
            run.writeFields && run.writeInitial && run.diagnosticFields ? (OutputLayout::fields() / run.outputPrefix).string() : std::string{});
        buildElectrosprayActiveInterfaceGeometry(mesh, fluidDomain, freeInterface, initialGeometryOutput);
        if (meshDefinition.runtimeAdaptation.enabled && !freeInterface.surface.geometry()->usesNativeLeaves())
            throw std::runtime_error("运行时 AMR 要求初始网格具有原生粗细层级并准备完整 Q4 液面，不能从同层 Q3 数据静默转换");
        if (workflow.problem == ElectrosprayProblem::MovingDielectric)
        {
            stage = "moving-dielectric-hydrodynamics";
            runElectrosprayMovingDielectricHydro(mesh, fluidDomain, electric, inputs, workflow, *cfg_.units, freeInterface);
            return 0; // Explicit q=0 moving interface; conducting charge is not advanced.
        }
        if (workflow.problem == ElectrosprayProblem::FixedInterfaceCapillary || workflow.problem == ElectrosprayProblem::FixedInterfaceDielectric)
        {
            // Preserve the old normal-only geometry and q=0 force numerics in
            // their explicit frozen branches. Coupled mode below builds ONE
            // local geometry for both capillary and charged Maxwell loads.
            stage = "fitted-capillary-force";
            const auto capillary = assembleElectrosprayCapillarity(mesh, fluidDomain, freeInterface.surface, liquid.densitySI, surfaceTensionSI,
                flowConfig.capillaryHalfWidthCells, *cfg_.units, 0.0, HydroForceEvolution::Frozen);
            writeElectrosprayCapillaryDiagnostics(mesh, capillary, liquid.densitySI, flowConfig.capillaryHalfWidthCells, run, *cfg_.units);
            if (workflow.problem == ElectrosprayProblem::FixedInterfaceCapillary)
            {
                stage = "frozen-interface-capillary-hydrodynamics";
                runElectrosprayFrozenCapillaryHydro(mesh, fluidDomain, electric, inputs, workflow, *cfg_.units,
                    capillary.forceX, capillary.forceY, capillary.source);
                return 0; // Only the explicitly selected immutable-interface NS interval.
            }
            stage = "fitted-electric-force";
            const auto forces = assembleElectrosprayDielectricForces(mesh, fluidDomain, electric, capillary, liquid.densitySI, electricConfig, *cfg_.units);
            writeElectrosprayDielectricForceDiagnostics(mesh, fluidDomain, forces, capillary, liquid.densitySI, run, *cfg_.units);
            stage = "frozen-interface-ehd-hydrodynamics";
            const auto load = forces.total.density();
            runElectrosprayFrozenEhdHydro(mesh, fluidDomain, electric, inputs, workflow, *cfg_.units,
                *load.faceX, *load.faceY, load.source);
            return 0; // Explicit fixed initial E/interface; no charge or shape evolution.
        }
        // Do not resurrect k=3 / beta=1 manufactured Helmholtz, its SparseLU
        // solve, or ignored finalTime as a spray run. Use the real boundaries.
        stage = "conducting-interface-state";
        // Retain the SAME active leaf scalar used by the free-surface loads.
        // The initial electric solve above still explicitly requires q=0;
        // do not infer conducting time advancement from this initial state.
        auto state = std::make_unique<ElectrosprayConductingState>(mesh);
        {
            stage = "fitted-charged-electrostatics";
            std::cout << "[电喷雾] 求解带电界面电场：使用已有 P3 跳跃重构、同单元曲率导数和 SuperLU；"
                      << "取样只使用有效流体单元，固体电极采用真实壁面边界。\n" << std::flush;
            auto electrical=initializeElectrosprayConductingState(mesh,freeInterface.surface.geometry(),
                *freeInterface.surface.panels(),freeInterface.surface.pinnedCoordinates(),
                electricConfig,liquid,gas,electricLinearTolerance,*cfg_.units,electric);
            const auto surfaceCharge=electrical.charge();
            const auto& chargePanels=surfaceCharge->geometry().panels();
            std::cout << "[电喷雾] 导电液面状态初始化：" << chargePanels.size()
                      << " 个面元，与流体共用同一液面和回转面积；初始自由电荷为零，尚未推进电荷。\n" << std::flush;
            const auto chargeSource=electrical.source();
            const auto& conducting=electrical.electric();
            std::cout << "[电喷雾] 拟合域 P3 电势/电流已求得：相对残差=" << conducting.report.relativeResidual
                      << "，界面修正面数=" << conducting.report.correctedPhysicalFaces
                      << "；尚未推进导电时间步。\n" << std::flush;
            if (run.writeFields && run.writeInitial && run.diagnosticFields)
            {
                stage = "conducting-initial-output";
                writeElectrosprayConductingInitialOutput(mesh,electric,conducting,chargeSource,*cfg_.units,run.outputPrefix);
            }
            stage = "fitted-charged-maxwell-loads";
            std::cout << "[电喷雾] 装配完整 Maxwell 与毛细载荷：共用实际压力梯度、液面局部曲率及流体侧壁面标量梯度。\n"
                      << std::flush;
            auto conductingForces = buildElectrosprayConductingForces(mesh,fluidDomain,
                freeInterface.surface.geometry(),surfaceCharge->geometry(),chargeSource,conducting,liquid, gas, surfaceTensionSI, electricConfig,flowConfig.capillaryHalfWidthCells,*cfg_.units,0.0);
            reportElectrosprayConductingForces(conductingForces,*cfg_.units);
            if (run.writeFields && run.writeInitial && run.diagnosticFields)
            {
                stage = "conducting-initial-load-output";
                writeElectrosprayConductingInitialLoads(mesh,fluidDomain,conductingForces,liquid, gas,*cfg_.units,run.outputPrefix);
            }
            state->electrical = std::move(electrical);
            state->forces = std::move(conductingForces);
        }
        stage = "conducting-time-coupling";
        std::unique_ptr<QuadtreeMeshRefiner> runtimeRefiner;
        std::unique_ptr<ElectrosprayMeshAdaptation> adaptation;
        if (meshDefinition.runtimeAdaptation.enabled)
        {
            runtimeRefiner = std::make_unique<QuadtreeMeshRefiner>(meshDefinition);
            if (!runtimeRefiner->attachGeometry(&definition, &geometry))
                throw std::runtime_error("运行时 AMR 缺少有效的公共几何加密规则");
            adaptation = std::make_unique<ElectrosprayMeshAdaptation>(ElectrosprayMeshAdaptation{
                *runtimeRefiner, definition, meshDefinition.runtimeAdaptation.everySteps,
                meshDefinition.runtimeAdaptation.maxSweeps});
            adaptation->currentPhases = !meshDefinition.geometryRefinement.phaseRules.empty();
            for (const auto& name : meshDefinition.physicsRefinement.indicators)
            {
                RefinementIndicatorField field;
                field.name = name;
                const auto& settings = meshDefinition.physicsRefinement;
                const ElectrosprayHydroScales scales(liquid.densitySI, *cfg_.units);
                field.reference = name == "Velocity" ? settings.velocityThresholdSI / scales.velocity
                    : settings.velocityGradientThresholdSI * scales.time;
                adaptation->flowIndicators.push_back(std::move(field));
            }
            std::cout << "[电喷雾 AMR] 已启用动态网格：每 " << adaptation->everySteps
                << " 个已接受步检查一次，每次最多 " << adaptation->maxSweeps
                << " 轮原四叉树细粗化；完整液面、固壁和开口模板保留保护层。\n";
            if (!adaptation->flowIndicators.empty())
                ElectrosprayMeshDiagnostics::flowPolicy(meshDefinition.physicsRefinement.velocityThresholdSI,
                    meshDefinition.physicsRefinement.velocityGradientThresholdSI);
        }
        runElectrosprayConductingHydro(mesh,fluidDomain,electric,inputs,workflow,*cfg_.units,
            freeInterface,std::move(state),
            versioned && versioned->chargeRelaxation ? &*versioned->chargeRelaxation : nullptr, adaptation.get());
        return 0; // Only the requested interval, not a physical spray acceptance claim.
    }
    catch (const std::exception& error)
    {
        std::cerr << "[电喷雾] 计算失败；阶段=" << stage << "；原因：" << error.what() << '\n';
        return 3;
    }
}
