#include "Electrospray/ElectrosprayTimeLoop.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Runtime/RunControls.h"
#include "Electrospray/Setup/ElectrosprayWorkflow.h"
#include <algorithm>
#include "Electrospray/Setup/ElectrosprayModuleInputs.h"
#include "Electrospray/Coupling/ElectrosprayMeshAdaptation.h"
#include "Electrospray/IO/ElectrosprayMeshDiagnostics.h"
#include "Electrospray/IO/ElectrosprayHydroDiagnostics.h"
#include "Electrospray/IO/ElectrosprayForceDiagnostics.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Electrospray/IO/ElectrosprayElectricDiagnostics.h"
#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"
#include "Electrospray/IO/ElectrosprayInterfaceMotionDiagnostics.h"
#include "Electrospray/Coupling/ElectrosprayCapillarity.h"
#include "Electrospray/Coupling/ElectrosprayForces.h"
#include "Electrospray/Coupling/ElectrosprayConducting.h"
#include "Electrospray/Setup/ElectrosprayMaterials.h"
#include "Electrospray/IO/ElectrosprayTecplotOutput.h"
#include "Electrospray/IO/ElectrosprayVolumeDiagnostics.h"
#include "Electrospray/IO/ElectrosprayChargeDiagnostics.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Interface/Compatibility/InterfacePhaseFlux.h"
#include "Core/IO/OutputLayout.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"

#include <filesystem>
#include <optional>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>


namespace
{
using Domain = HydroFittedDomain;

void runHydro(QuadtreeMesh& mesh, const HydroFittedDomain& originalDomain,
    const ElectrosprayElectricState& originalBinding, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units,
    const VEdgeField<double>* forceX, const HEdgeField<double>* forceY,
    std::shared_ptr<const HydroForceSource> frozenSource = {},
    ElectrosprayInterfaceState* movingInterface = nullptr,
    std::unique_ptr<ElectrosprayConductingState> conducting = {},
    const Quantum::Electricity::SurfaceChargeRelaxationOptions* relaxation = nullptr,
    ElectrosprayMeshAdaptation* adaptation = nullptr,
    const ElectrosprayAcceptedObserver& acceptedObserver = {})
{
    const auto& run = inputs.controls;
    const auto& flowConfig = inputs.flow;
    const auto& electricConfig = inputs.electric;
    const auto& liquid = inputs.materials.phases.at(flowConfig.liquidPhase);
    const auto& gas = inputs.materials.phases.at(flowConfig.gasPhase);
    const double surfaceTensionSI = inputs.materials.surfaceTensionSI.at(std::minmax(flowConfig.liquidPhase, flowConfig.gasPhase));
    const double electricLinearTolerance = inputs.systems.at("ElectricPotential").relativeTolerance;
    // One NS/geometry/volume loop. Only the electric/charge candidate stage
    // differs between explicit q=0 and full conducting motion.
    const bool charged = static_cast<bool>(conducting);
    const bool moving = movingInterface != nullptr;
    const auto frozenInterface = frozenSource ? frozenSource->interfaceGeometry() : nullptr;
    if (charged && !moving) throw std::invalid_argument("conducting spray requires a moving interface");
    if (relaxation && !charged) throw std::invalid_argument("charge iteration controls require conducting motion");
    const bool capillary = forceX != nullptr || moving;
    const bool ehd = workflow.problem == ElectrosprayProblem::FixedInterfaceDielectric;
    const ElectrosprayHydroScales scale(liquid.densitySI, units);
    const auto checkedSteps = run.fixedStepCount();
    if (!checkedSteps)
        throw std::invalid_argument(moving ? "moving spray finalTime must be an integer number of fixed dt steps" : ehd
            ? "frozen-interface EHD Hydro finalTime must be an integer number of fixed dt steps" : capillary
            ? "frozen-interface capillary Hydro finalTime must be an integer number of fixed dt steps"
            : "frozen-material Hydro finalTime must be an integer number of fixed dt steps");
    const int stepCount = *checkedSteps;
    auto domain = std::make_unique<HydroFittedDomain>(originalDomain);
    auto initial = std::make_unique<ElectrosprayElectricState>(originalBinding);
    auto hydro = std::make_unique<Hydrodynamics>(mesh);
    std::unique_ptr<ElectrosprayInterfaceState> interfaceOwner;
    if (adaptation)
    {
        if (!charged || !movingInterface->surface.geometry()->usesNativeLeaves())
            throw std::invalid_argument("runtime AMR requires the native-Q4 conducting workflow");
        interfaceOwner = std::make_unique<ElectrosprayInterfaceState>(std::move(*movingInterface));
        movingInterface = interfaceOwner.get();
    }
    std::unique_ptr<ElectrosprayElectricState> movingElectric;
    std::unique_ptr<ElectrosprayDielectricForces> movingForces;
    ElectrosprayMaterialState conductingMaterials;
    RunControls quietForceRun = run;
    quietForceRun.writeFields = false; // Never overwrite saved initial-force files.
    if (charged)
    {
        const auto preflight = preflightElectrosprayConduction(inputs, workflow, units);
        if (!conducting->electrical.charge() || !conducting->electrical.electric().report.success
            || !conducting->electrical.electric().report.conductionCurrent)
            throw std::invalid_argument("conducting spray requires its solved initial charge/current state");
        conductingMaterials = prepareElectrosprayMaterials(mesh,*domain,*movingInterface->surface.geometry(),
            flowConfig, electricConfig, liquid, gas,units,*initial);
        forceX = conducting->forces.total.density().faceX; forceY = conducting->forces.total.density().faceY;
        std::cout << "[电喷雾] " << (relaxation ? "步内欠松弛迭代；电荷松弛数=" : "导电显式松弛数=") << preflight.relaxationCfl
            << "，上限沿用 0.25；界面体积和电荷均保留原守恒检查。\n";
    }
    else if (moving)
    {
        movingElectric = std::make_unique<ElectrosprayElectricState>(updateElectrosprayDielectricState(
            mesh, *domain, movingInterface->surface.geometry(), flowConfig, electricConfig, liquid, gas, electricLinearTolerance, units, *initial,0.0));
        int synchronizedCells = 0;
        for (int k : domain->cellIndexByRow()) synchronizedCells += movingElectric->phaseByCell[k] != initial->phaseByCell[k];
        std::cout << "[电喷雾] 初始液面与物性同步：更新相标记的单元=" << synchronizedCells
            << "；尚未推进物理时间。\n";
        auto initialFrameRun = quietForceRun;
        initialFrameRun.writeFields = run.writeFields && run.writeInitial && run.diagnosticFields;
        initialFrameRun.outputPrefix = run.outputPrefix + "_moving_dielectric_step0";
        if (run.writeFields && run.writeInitial && run.diagnosticFields) writeElectrosprayElectrostatics(mesh, *movingElectric, units,
            (OutputLayout::fields() / initialFrameRun.outputPrefix).string(), false, run.detailedLog);
        const auto cap = assembleElectrosprayCapillarity(mesh, *domain, movingInterface->surface, liquid.densitySI, surfaceTensionSI,
            flowConfig.capillaryHalfWidthCells, units, 0.0, HydroForceEvolution::TimeDependent);
        writeElectrosprayCapillaryDiagnostics(mesh, cap, liquid.densitySI, flowConfig.capillaryHalfWidthCells, initialFrameRun, units);
        movingForces = std::make_unique<ElectrosprayDielectricForces>(assembleElectrosprayDielectricForces(
            mesh, *domain, *movingElectric, cap, liquid.densitySI, electricConfig, units));
        writeElectrosprayDielectricForceDiagnostics(mesh, *domain, *movingForces, cap, liquid.densitySI, initialFrameRun, units);
        forceX = movingForces->total.density().faceX; forceY = movingForces->total.density().faceY;
    }
    // Never retain references into a candidate bundle across its commit.
    const auto materialPhases = [&]() -> const std::vector<int>& {
        return charged ? conductingMaterials.phaseByCell : moving ? movingElectric->phaseByCell : initial->phaseByCell;
    };
    updateElectrosprayHydroMaterials(mesh, *domain, *hydro, materialPhases(), liquid, gas, scale, 0.0,
        moving ? HydroMaterialEvolution::TimeDependent : HydroMaterialEvolution::Frozen);
    hydro->initializeState(mesh, HydroStepGeometry::fittedRZ(*domain));
    HydroSolveOptions options;
    options.fittedNumerics = electrosprayHydroNumerics(inputs);
    const std::string prefix = (OutputLayout::series() / run.outputPrefix).string();
    const std::string fieldBase = (OutputLayout::fields() / run.outputPrefix).string();
    const std::string fieldPrefix = charged ? fieldBase + "_conducting" : moving ? fieldBase + "_moving_dielectric" : ehd ? fieldBase + "_frozen_ehd"
        : capillary ? fieldBase + "_frozen_capillary" : fieldBase;
    std::unique_ptr<ElectrosprayTecplotOutput> tecplot;
    if (run.writeFields)
        tecplot = std::make_unique<ElectrosprayTecplotOutput>(std::filesystem::path(fieldPrefix).filename().string(),run.diagnosticFields);
    ElectrosprayHydroDiagnostics diagnostics(prefix, charged, moving, ehd, capillary);
    std::unique_ptr<ElectrosprayMeshDiagnostics> meshDiagnostics;
    if (adaptation) meshDiagnostics = std::make_unique<ElectrosprayMeshDiagnostics>(prefix, scale, units);
    std::unique_ptr<ElectrosprayVolumeDiagnostics> volumeDiagnostics;
    std::unique_ptr<ElectrosprayChargeDiagnostics> chargeDiagnostics;
    if (charged) chargeDiagnostics = std::make_unique<ElectrosprayChargeDiagnostics>(prefix,units,relaxation!=nullptr);
    long double initialVolume = 0.0L, expectedVolumeChange = 0.0L;
    long double acceptedVolume = 0.0L, acceptedVolumeDefect = 0.0L;
    double previousLiquidFlux = 0.0;
    if (moving)
    {
        std::cout << "[电喷雾] 准备初始逐单元液体体积与相通量，用于本步守恒修正。" << std::endl;
        const auto initialVelocity = hydro->fittedVelocity(mesh, *domain, HydroVelocityTimeLevel::Current,
            HydroFittedVelocitySampler::Reconstruction::RZFluxCompatible);
        movingInterface->surface.initializePhase(mesh, HydroInterfaceTransport::domainView(*domain), HydroInterfaceTransport::normalFlux(initialVelocity));
        initialVolume = (movingInterface->surface.geometry()->usesNativeLeaves() ? std::visit([](const auto& measure) { return measure.summedNegativeVolume; }, movingInterface->surface.phase().volumeMeasure)
            : InterfacePhaseGeometry::volume(mesh,domain->geometry(),*movingInterface->surface.geometry(),*movingInterface->surface.panels()))
            * scale.length * scale.length * scale.length;
        std::ostringstream volumeMessage;
        volumeMessage << std::scientific << std::setprecision(6)
            << "[电喷雾] 初始液体体积=" << initialVolume << " m^3；包含自由液面、湿润固壁和开口围成的体积。\n";
        std::cout << volumeMessage.str();
        diagnostics.initializeMotion(prefix);
        volumeDiagnostics = std::make_unique<ElectrosprayVolumeDiagnostics>(prefix, scale.length);
        std::cout << "[电喷雾] 局部体积修正：每次重算液面与相通量；迭代准则为相对输运收支 1e-3、单元容量归一化残差 1e-10。"
            "这是新增迭代准则，不替代原全局体积诊断或物理验收。\n";
        if (charged)
            std::cout << "[电喷雾] 模式：导电运动液面；推进 NS、液面、物性、液面电荷和 P3 电场。"
                << (relaxation ? "电荷与电场在步内欠松弛迭代，两层历史沿同一材料轨迹搬运；"
                               : "旧时刻电流只积分一次，材料轨迹包含界面对流；")
                << "壁面无滑移、唇缘位置钉扎。\n";
        else
            std::cout << "[电喷雾] 模式：运动介电液面；推进NS、液面、物性与静电场。"
                "q=0，尚未推进导电电荷；壁面无滑移、唇缘位置钉扎。\n";
    }
    else if (ehd)
        std::cout << "[电喷雾] 模式：固定介电液面；NS含毛细力和极化电力，液面、物性、电场冻结，未推进电荷。\n";
    else if (capillary)
        std::cout << "[电喷雾] 模式：固定毛细液面；仅NS和毛细力，无电力，液面与物性冻结。\n";
    else
        std::cout << "[电喷雾] 模式：固定物性的压力驱动流；无液面推进、毛细力或电力。\n";
    {
        std::ostringstream message;
        message << std::scientific << std::setprecision(4)
            << "[电喷雾] 开始推进：总步数=" << static_cast<int>(stepCount)
            << "，时间步=" << run.dt*scale.time << " s，终时=" << run.endTime*scale.time
            << " s；线性求解器=SuperLU，logEvery=" << run.logEverySteps << "。\n"
            << "[电喷雾] 流量约定：供液量以流入为正；侧面/顶部以流出为正，负值表示回流。\n";
        std::cout << message.str() << std::flush;
    }
    if (run.writesFieldsAtStep(0,stepCount))
    {
        if (run.diagnosticFields)
            writeElectrosprayHydroFields(mesh, *domain, *hydro, materialPhases(), scale, fieldPrefix + "_step0");
        if (charged)
            tecplot->writeConducting(mesh,*domain,*hydro,materialPhases(),*initial,conducting->electrical,liquid, gas, surfaceTensionSI,units,*movingInterface,conducting->forces,0,0.0);
        else
            tecplot->write(mesh, *domain, *hydro, moving ? *movingElectric : *initial, liquid, surfaceTensionSI, units, movingInterface,
                movingForces.get(), forceX, forceY, 0, 0.0, workflow.problem);
    }
    int acceptedSteps = 0, outputSteps = 0;
    for (int completed = 0; completed < stepCount; ++completed)
    {
        const int step = completed + 1;
        const char* stage = "mesh-adaptation";
        try
        {
            // AMR belongs BETWEEN accepted physical steps. Prepare every target
            // object before touching the accepted mesh; no extra physical dt or
            // charge-current integration is used. Spatial frames remain one/step.
            if (adaptation && step > 1 && (step - 1) % adaptation->everySteps == 0)
            {
                std::cout << "[电喷雾 AMR] 第 " << step - 1
                    << " 步之后检查移动加密带；不推进物理时间。\n" << std::flush;
                std::unique_ptr<ElectrosprayMeshCandidate> candidate;
                try
                {
                    candidate = prepareElectrosprayMeshAdaptation(*adaptation,
                        {mesh, *domain, *initial, *movingInterface, *hydro, *conducting, conductingMaterials}, inputs, units,
                        ElectrosprayMeshDiagnostics::progress);
                }
                catch (const std::exception& error)
                {
                    throw std::runtime_error("动态 AMR 候选失败，未替换第 " + std::to_string(step - 1)
                        + " 步的网格与耦合状态；原因：" + error.what());
                }
                if (candidate)
                {
                    candidate->commit(mesh,domain,initial,interfaceOwner,hydro,conducting,conductingMaterials);
                    movingInterface = interfaceOwner.get();
                    const auto& transfer = candidate->report();
                    previousLiquidFlux = static_cast<double>(transfer.liquidFluxAfter
                        * scale.velocity * scale.length * scale.length);
                    acceptedVolume = transfer.volumeAfter * scale.length * scale.length * scale.length;
                    // Keep the original initial volume and cumulative transport:
                    // AMR must not reset the existing conservation diagnostics.
                    acceptedVolumeDefect = acceptedVolume - initialVolume - expectedVolumeChange;
                    forceX = conducting->forces.total.density().faceX;
                    forceY = conducting->forces.total.density().faceY;
                    meshDiagnostics->write(step - 1, (step - 1) * run.dt, transfer);
                }
                else std::cout << "[电喷雾 AMR] 网格已满足当前加密与保护要求，本次无需改变。\n" << std::flush;
            }
            if (charged)
            {
                forceX = conducting->forces.total.density().faceX; forceY = conducting->forces.total.density().faceY;
                std::cout << "[电喷雾] 导电步 " << step << '/' << static_cast<int>(stepCount)
                    << "：先用已接受的毛细/Maxwell 载荷求解 NS。\n" << std::flush;
            }
            HydroForceDensity load;
            load.evaluationTime = (step - 1) * run.dt;
            if (capillary)
            {
                load.faceX = forceX; load.faceY = forceY;
                // Carry the producer's record; never manufacture a fresh source
                // from the requested step or reinterpret old AMR fields as new.
                load.source = charged ? conducting->forces.total.source() : moving ? movingForces->total.source() : frozenSource;
            }
            stage = "hydro-step";
            auto hydroCandidate = hydro->prepareStep(mesh, {HydroStepGeometry::fittedRZ(*domain),
                {step, (step - 1) * run.dt, run.dt}, options, load,
                moving ? movingInterface->surface.geometry().get() : frozenInterface.get()});
            const auto& report = hydroCandidate.report();
            if (!report.success)
            {
                std::cerr << std::setprecision(17) << "[电喷雾] NS步骤失败，未提交当前步：step=" << step
                    << ", acceptedSteps=" << step - 1 << ", SuperLU(momentum,pressure,coupled)="
                    << report.momentumUsedSuperLU << ',' << report.pressureUsedSuperLU << ',' << report.coupledUsedSuperLU
                    << ", momentumResidual=" << report.momentumRelativeResidual
                    << ", projectionResidual=" << report.pressureRelativeResidual
                    << ", coupledResidual=" << report.coupledRelativeResidual
                    << ", fullMomentumResidual=" << report.momentumBalanceRelativeResidual
                    << ", maxDiv_code=" << report.divergenceAfterProjection << '\n';
                throw std::runtime_error("frozen Hydro step " + std::to_string(step)
                    + " failed at " + report.stage + ": " + report.reason);
            }
            std::optional<ElectrosprayBoundaryFlowRatesSI> boundaryFlow;
            if (moving)
            {
                stage = "interface-candidate";
                InterfaceGeometryObserver observeGeometry;
                if (step == 1 || step % run.logEverySteps == 0 || step == static_cast<int>(stepCount))
                    observeGeometry = [&](InterfaceMotionStage stage, const LevelSetInterfacePanels& panels) {
                        writeElectrosprayInterfaceCurvature(step, stage, panels, units);
                    };
                const double time = step * run.dt;
                auto candidate = prepareElectrosprayInterfaceMotion(mesh, *domain, *movingInterface, *hydro,
                    {step, (step-1)*run.dt, time, run.dt}, flowConfig.capillaryHalfWidthCells,
                    inputs.systems.at("InterfaceVolumeCorrection").relativeTolerance, workflow.legacyIterationLog.value_or(run.detailedLog)
                        ? InterfaceVolumeProgress(ElectrosprayVolumeDiagnostics::writeProgress) : InterfaceVolumeProgress{},
                    observeGeometry);
                const bool writeStep = run.writesFieldsAtStep(step,stepCount);
                auto frameRun = quietForceRun;
                frameRun.writeFields = writeStep && run.diagnosticFields;
                frameRun.outputPrefix = run.outputPrefix + (charged ? "_conducting_step" : "_moving_dielectric_step") + std::to_string(step);
                std::unique_ptr<ElectrosprayElectricState> nextElectric;
                std::unique_ptr<ElectrosprayDielectricForces> nextForces;
                std::unique_ptr<HydroFittedCapillaryState> nextCapillary;
                std::unique_ptr<ElectrosprayConductingCandidate> nextConducting;
                ElectrosprayMaterialState nextMaterials;
                // Hold this step's actual OLD current independently of the
                // charge owner that will be replaced at commit.
                const auto sourceCurrent = charged ? conducting->electrical.electric().report.conductionCurrent : nullptr;
                stage = "electric-force-candidate";
                if (charged)
                {
                    std::cout << "[电喷雾] 导电步 " << step
                        << "：液面与局部体积修正完成，更新电荷并重算 P3 电场和完整受力。\n" << std::flush;
                    nextMaterials = prepareElectrosprayMaterials(mesh,*domain,*candidate.state.surface.geometry(),
                        flowConfig, electricConfig, liquid, gas,units,*initial);
                    nextConducting = std::make_unique<ElectrosprayConductingCandidate>(prepareElectrosprayConductingStep(
                        mesh,*domain,*initial,inputs,workflow,units,*conducting,candidate,time,relaxation));
                }
                else
                {
                    nextElectric = std::make_unique<ElectrosprayElectricState>(updateElectrosprayDielectricState(
                        mesh,*domain,candidate.state.surface.geometry(),flowConfig, electricConfig, liquid, gas, electricLinearTolerance,units,*movingElectric,time));
                    // Direct initialization preserves C++17 copy elision for
                    // Field arrays and retains diagnostics until acceptance.
                    nextCapillary.reset(new HydroFittedCapillaryState(assembleElectrosprayCapillarity(
                        mesh, *domain, candidate.state.surface, liquid.densitySI, surfaceTensionSI,
                        flowConfig.capillaryHalfWidthCells, units, time, HydroForceEvolution::TimeDependent)));
                    nextForces = std::make_unique<ElectrosprayDielectricForces>(assembleElectrosprayDielectricForces(
                        mesh,*domain,*nextElectric,*nextCapillary,liquid.densitySI, electricConfig,units));
                }
                stage = "candidate-checks";
                const auto& scalar = candidate.state.surface.geometry();
                if (candidate.state.surface.stage() != InterfaceState::Stage::Ready || !scalar
                    || !candidate.state.surface.panels() || !candidate.state.surface.panels()->usesSource(*scalar))
                    throw std::logic_error("coupled candidate requires ready interface geometry/panels");
                const auto& nextSource = charged ? nextConducting->state->forces.total.source() : nextForces->total.source();
                if (!nextSource || nextSource->sourceTime() != time
                    || nextSource->evolution() != HydroForceEvolution::TimeDependent
                    || nextSource->interfaceGeometry() != scalar)
                    throw std::logic_error("coupled candidate force source differs from interface/time");
                if (charged)
                {
                    const auto& electrical = nextConducting->state->electrical;
                    if (electrical.completedStep() != step || electrical.time() != time
                        || electrical.electric().potential.source() != scalar
                        || !electrical.charge()->geometry().usesSource(*scalar))
                        throw std::logic_error("coupled charge candidate differs from interface/step/time");
                }
                else if (nextElectric->dielectric.time() != time || nextElectric->dielectric.potential().source() != scalar)
                    throw std::logic_error("coupled dielectric candidate differs from interface/time");
                const auto& nextPhases = charged ? nextMaterials.phaseByCell : nextElectric->phaseByCell;
                const bool electricSuperLU = charged ? nextConducting->state->electrical.electric().report.usedSuperLU : nextElectric->dielectric.report().usedSuperLU;
                const double electricResidual = charged ? nextConducting->state->electrical.electric().report.relativeResidual : nextElectric->dielectric.report().finalResidual;
                const long double volume = (candidate.state.surface.geometry()->usesNativeLeaves()
                    ? std::visit([](const auto& measure) { return measure.summedNegativeVolume; }, candidate.state.surface.phase().volumeMeasure) : InterfacePhaseGeometry::volume(mesh,domain->geometry(),
                        *candidate.state.surface.geometry(),*candidate.state.surface.panels())) * scale.length * scale.length * scale.length;
                const double liquidFlux = static_cast<double>(hydro->fittedPhaseOutwardFlux(mesh, *domain, *candidate.state.surface.geometry()) * scale.velocity * scale.length * scale.length);
                const auto& nextPhaseBudget = candidate.state.surface.phase().phaseBudget;
                const long double fluxFactor = static_cast<long double>(scale.velocity)*scale.length*scale.length;
                const long double ledgerBoundary = nextPhaseBudget.boundaryOutward*fluxFactor;
                const long double pressureDifference = ledgerBoundary-liquidFlux;
                if (std::abs(nextPhaseBudget.assemblyClosure) > flowConfig.phaseFluxBalance*nextPhaseBudget.absoluteFaceFlux
                    || std::abs(pressureDifference) > flowConfig.phaseFluxBalance
                        * (std::abs(ledgerBoundary)+std::abs(liquidFlux)))
                    throw std::runtime_error("shared phase flux ledger failed assembly/pressure-port closure");
                const long double nextExpected = expectedVolumeChange
                    - 0.5L * run.dt * scale.time * (previousLiquidFlux + liquidFlux);
                int changedMaterials = 0;
                for (int k : domain->cellIndexByRow()) changedMaterials += nextPhases[k] != materialPhases()[k];
                stage = "hydro-material-candidate";
                updateElectrosprayHydroMaterials(mesh, *domain, *hydro, nextPhases, liquid, gas, scale,
                    step * run.dt, HydroMaterialEvolution::TimeDependent);
                boundaryFlow = evaluateElectrosprayBoundaryFlow(mesh, *domain, *hydro, flowConfig, scale);
                // Hydro has already advanced locally: failures stop this run,
                // without retry or rollback. Remaining publication exchanges
                // owners and ledgers without allocating/copying field arrays.
                static_assert(std::is_nothrow_swappable_v<ElectrosprayInterfaceState>);
                static_assert(std::is_nothrow_swappable_v<ElectrosprayMaterialState>);
                hydroCandidate.accept(mesh);
                [&]() noexcept {
                    std::swap(*movingInterface, candidate.state);
                    if (charged)
                    {
                        conducting.swap(nextConducting->state);
                        std::swap(conductingMaterials, nextMaterials);
                        forceX = conducting->forces.total.density().faceX; forceY = conducting->forces.total.density().faceY;
                    }
                    else
                    {
                        movingElectric.swap(nextElectric);
                        movingForces.swap(nextForces);
                        forceX = movingForces->total.density().faceX; forceY = movingForces->total.density().faceY;
                    }
                    expectedVolumeChange = nextExpected; previousLiquidFlux = liquidFlux;
                    acceptedVolume = volume;
                    acceptedVolumeDefect = volume-initialVolume-expectedVolumeChange;
                    acceptedSteps = step;
                }();
                stage = "accepted-motion-output";
                const auto& phaseBudget = movingInterface->surface.phase().phaseBudget;
                if (!charged)
                {
                    writeElectrosprayCapillaryDiagnostics(mesh, *nextCapillary, liquid.densitySI, flowConfig.capillaryHalfWidthCells, frameRun, units);
                    writeElectrosprayDielectricForceDiagnostics(mesh, *domain, *movingForces, *nextCapillary, liquid.densitySI, frameRun, units);
                }
                if (writeStep && run.diagnosticFields)
                {
                    writeElectrosprayInterfaceMotion(mesh, *movingInterface, candidate.motion.report, units,
                        fieldPrefix + "_step" + std::to_string(step));
                    if (!charged)
                        writeElectrosprayElectrostatics(mesh, *movingElectric, units, fieldPrefix + "_step" + std::to_string(step),
                            false, run.detailedLog);
                    writeElectrosprayPhaseFluxBudget(*domain, phaseBudget, scale, fieldPrefix + "_step" + std::to_string(step));
                    // The exchanged candidate owns the OLD interface and ledgers.
                    volumeDiagnostics->writeFields(fieldPrefix + "_step" + std::to_string(step), *domain,
                        candidate.state.surface.phase().volumeMeasure, candidate.state.surface.phase().phaseBudget,
                        movingInterface->surface.phase().volumeMeasure, movingInterface->surface.phase().phaseBudget, run.dt);
                }
                if (charged)
                    chargeDiagnostics->writeStep(step,time,*conducting,*sourceCurrent,
                        nextConducting->chargeStep,nextConducting->report);
                ElectrosprayMotionStepValues motionValues;
                motionValues.volumeSI = volume;
                motionValues.initialVolumeSI = initialVolume;
                motionValues.expectedVolumeChangeSI = expectedVolumeChange;
                motionValues.liquidFluxSI = liquidFlux;
                motionValues.changedMaterials = changedMaterials;
                motionValues.electricSuperLU = electricSuperLU;
                motionValues.electricResidual = electricResidual;
                motionValues.ledgerBoundarySI = ledgerBoundary;
                motionValues.fluxFactorSI = fluxFactor;
                motionValues.pressureDifferenceSI = pressureDifference;
                diagnostics.writeAcceptedMotion(step, time, candidate.motion.report, phaseBudget, motionValues, scale);
                const auto& conservation = candidate.motion.conservation;
                volumeDiagnostics->writeStep(step, time*scale.time, conservation);
                if (step == 1 || step % run.logEverySteps == 0 || step == static_cast<int>(stepCount))
                {
                    std::cout << "[电喷雾] 局部体积修正完成：迭代=" << conservation.iterations
                        << "，回退=" << conservation.backtracks << "，相对输运收支="
                        << conservation.initial.relativeBalance << " -> " << conservation.final.relativeBalance << '\n';
                }
                if (run.detailedLog)
                    std::cout << std::setprecision(17) << "[Electrospray] Moving interface accepted: step=" << step
                    << ", changedPhi=" << candidate.motion.report.changedValues << ", maxPhiIncrement_code=" << candidate.motion.report.maxValueIncrement
                    << ", traceRatio=" << candidate.motion.report.maxTraceRatio << ", liquidVolume_m3=" << volume
                    << ", volumeDefect_m3=" << volume-initialVolume-expectedVolumeChange << '\n';
            }
            else
            {
                boundaryFlow = evaluateElectrosprayBoundaryFlow(mesh, *domain, *hydro, flowConfig, scale);
                hydroCandidate.accept(mesh);
                acceptedSteps = step;
            }
            stage = "accepted-hydro-output";
            const double time = step * run.dt;
            diagnostics.writeStep(step, time, report, *boundaryFlow, scale);
            if (step == 1 || step % run.logEverySteps == 0 || step == static_cast<int>(stepCount))
            {
                std::ostringstream message;
                message << std::scientific << std::setprecision(4)
                    << "[电喷雾] 步 " << step << '/' << static_cast<int>(stepCount) << " | t=" << time*scale.time
                    << " s | 最大流速=" << report.maxSpeed*scale.velocity << " m/s | max|div U|="
                    << report.divergenceAfterProjection/scale.time << " 1/s | 动量相对残差="
                    << report.momentumBalanceRelativeResidual << '\n'
                    << "  体积流量[m^3/s]：供液=" << -boundaryFlow->supplyOutward << "，侧面出流=" << boundaryFlow->ambientOutward
                    << "，顶部出流=" << boundaryFlow->outletOutward << "，净向外=" << boundaryFlow->netOutward << '\n';
                if (moving)
                {
                    message << "  液体体积=" << acceptedVolume << " m^3，体积收支缺陷=" << acceptedVolumeDefect << " m^3";
                    if (expectedVolumeChange != 0.0L)
                        message << "，占累计净液体输运量=" << std::scientific << std::setprecision(4)
                            << 100.0L*std::abs(acceptedVolumeDefect/expectedVolumeChange) << '%';
                    else message << "，累计净液体输运量为零（不计算百分比）";
                    message << "；守恒诊断，不代表已通过物理验收。\n";
                }
                std::cout << message.str() << std::flush;
            }
            if (run.writesFieldsAtStep(step,stepCount))
            {
                stage = "accepted-field-output";
                if (run.diagnosticFields)
                    writeElectrosprayHydroFields(mesh, *domain, *hydro, materialPhases(), scale, fieldPrefix + "_step" + std::to_string(step));
                if (charged)
                    tecplot->writeConducting(mesh,*domain,*hydro,materialPhases(),*initial,conducting->electrical,liquid, gas, surfaceTensionSI,units,*movingInterface,conducting->forces,step,time);
                else
                    tecplot->write(mesh, *domain, *hydro, moving ? *movingElectric : *initial, liquid, surfaceTensionSI, units, movingInterface,
                        movingForces.get(), forceX, forceY, step, time, workflow.problem);
                std::cout << "[Electrospray] Tecplot 场输出完成: step=" << step << ", time_s=" << time*scale.time
                    << ", directory=output/fields; 本帧为已完成的计算步。\n";
            }
            outputSteps = step;
            if (acceptedObserver)
            {
                stage = "accepted-state-observer";
                acceptedObserver(step, {mesh, *domain, *initial, *movingInterface, *hydro,
                    *conducting, conductingMaterials});
            }
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error("[Electrospray] step " + std::to_string(step) + " failed at " + stage
                + "; physical steps accepted through " + std::to_string(acceptedSteps)
                + ", output completed through " + std::to_string(outputSteps) + ": " + error.what());
        }
    }
    std::cout << "[电喷雾] 已完成 " << static_cast<int>(stepCount)
        << " 个时间步。场文件：output/fields；逐步诊断：output/series；帧清单：output/meta。\n"
        << (charged ? "[电喷雾] 请求的导电耦合区间已完成；仍需独立物理基准及网格/时间收敛验收。\n"
            : "[电喷雾] 本次为开发模式短程，不能据此判定完整导电电喷雾通过。\n");
}
} // namespace

void runElectrosprayFrozenMaterialHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units)
{
    if (workflow.problem != ElectrosprayProblem::PressureDriven) throw std::invalid_argument("pressure-only Hydro requires explicit frozenMaterialHydro mode");
    runHydro(mesh, domain, initial, inputs, workflow, units, nullptr, nullptr);
}

void runElectrosprayFrozenCapillaryHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units,
    const VEdgeField<double>& forceX, const HEdgeField<double>& forceY,
    std::shared_ptr<const HydroForceSource> source)
{
    if (workflow.problem != ElectrosprayProblem::FixedInterfaceCapillary)
        throw std::invalid_argument("fixed-interface capillary Hydro requires explicit frozenInterfaceCapillary mode");
    runHydro(mesh, domain, initial, inputs, workflow, units, &forceX, &forceY, std::move(source));
}

void runElectrosprayFrozenEhdHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units,
    const VEdgeField<double>& totalForceX, const HEdgeField<double>& totalForceY,
    std::shared_ptr<const HydroForceSource> source)
{
    if (workflow.problem != ElectrosprayProblem::FixedInterfaceDielectric)
        throw std::invalid_argument("fixed-interface EHD Hydro requires explicit frozenInterfaceEhd mode");
    runHydro(mesh, domain, initial, inputs, workflow, units, &totalForceX, &totalForceY, std::move(source));
}

void runElectrosprayMovingDielectricHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& initial, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units, ElectrosprayInterfaceState& interfaceState)
{
    if (workflow.problem != ElectrosprayProblem::MovingDielectric) throw std::invalid_argument("moving dielectric Hydro requires explicit movingDielectric mode");
    runHydro(mesh, domain, initial, inputs, workflow, units, nullptr, nullptr, {}, &interfaceState);
}

void runElectrosprayConductingHydro(QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const ElectrosprayElectricState& binding, const ElectrosprayModuleInputs& inputs,
    const ElectrosprayWorkflow& workflow, const UnitsContext& units, ElectrosprayInterfaceState& interfaceState,
    std::unique_ptr<ElectrosprayConductingState> initialState,
    const Quantum::Electricity::SurfaceChargeRelaxationOptions* relaxation,
    ElectrosprayMeshAdaptation* adaptation,
    const ElectrosprayAcceptedObserver& acceptedObserver)
{
    if (workflow.problem != ElectrosprayProblem::Conducting || !initialState)
        throw std::invalid_argument("conducting Hydro requires explicit coupled mode and an owned initial state");
    runHydro(mesh,domain,binding,inputs,workflow,units,nullptr,nullptr,{},&interfaceState,std::move(initialState),relaxation,adaptation,acceptedObserver);
}
