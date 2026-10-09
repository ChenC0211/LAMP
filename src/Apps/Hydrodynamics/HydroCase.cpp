#include "Hydrodynamics/HydroCase.h"
#include "Hydrodynamics/HydroMomentCase.h"
#include "Hydrodynamics/HydroConservativeInterface.h"
#include "Hydrodynamics/Coupling/HydroCaseCoupling.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroBoundaryVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceBoundary.h"
#include "Modules/Hydrodynamics/HydroInterfaceForce.h"
#include "Modules/Hydrodynamics/HydroFlowDiagnostics.h"
#include "Modules/Interface/Compatibility/InterfaceNodalState.h"
#include "Hydrodynamics/IO/HydroCaseInput.h"
#include "Modules/Hydrodynamics/IO/HydroGeometryLoader.h"
#include "Modules/Hydrodynamics/IO/HydroDiagnosticsWriter.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/OutputLayout.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Field/StaggeredVelocitySampler.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/HydroMidpointStep.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"

#include <cmath>
#include <array>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void updateMaterials(QuadtreeMesh& mesh, Hydrodynamics& hydro,
                     const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
                     const HydrodynamicsConfig& config, bool interfaceMoves, double sourceTime)
{
    const auto domain = HydroStepGeometry::complete(config.coordinateSystem == "axisymmetricRZ");
    if (geometry)
        hydro.updateMaterials(mesh, {domain, sourceTime,
            interfaceMoves ? HydroMaterialEvolution::TimeDependent : HydroMaterialEvolution::Frozen,
            HydroLevelSetMaterial{&geometry->levelSets[interfaceIndex], config.densityNegative, config.densityPositive,
                config.viscosityNegative, config.viscosityPositive, config.interfaceThickness}});
    else
        hydro.updateMaterials(mesh, {domain, sourceTime, HydroMaterialEvolution::Frozen,
            HydroConstantMaterial{config.density, config.viscosity}});
}

bool supportsAxisInterface(const QuadtreeMesh& mesh, const HydrodynamicsConfig& config,
    const InterfaceNodalConfig& interface)
{
    return config.coordinateSystem == "axisymmetricRZ" && mesh.minCoord.x == 0.0 &&
        config.interfaceForceGeometry == "hermite" && config.interfaceForceLayer == "regularized" &&
        interface.reinitializeEvery == 0;
}

}

using HydroFlowDiagnostics::collectDiagnostics;
using HydroFlowDiagnostics::collectForcePressureImbalance;
using HydroFlowDiagnostics::boundaryFlow;

struct HydroCase::Impl
{
    explicit Impl(HydroCaseCoupling* provider);
    int initialize();
    int advance(int step, double forceTime);
    int advanceConservative(int step,double forceTime,double proposedDt);
    int advanceMidpoint(int step,double begin);
    int advanceMoment(int step);
    int advanceNext();
    bool hasNextStep() const;
    void appendReport(int step, HydroSolveReport& report,double actualTime=-1,double actualDt=-1);
    void finish();
    void writeAdaptation(int step,int attempt,double time,int sourceCells,const char* status,
        const HydroConservativeInterface::MeshCandidate*,const HydroSolveReport&);

    enum class Lifecycle { Prepared, Ready, Failed, Complete };
    Lifecycle lifecycle = Lifecycle::Prepared;
    int completedSteps = 0;
    int acceptedStep = 0;
    double acceptedTime=0.0,outputTime=0.0,nextDt=0.0;
    int retryFreeSteps=0;
    HydroCaseCoupling* coupling; // Borrowed from the named App, never from Hydro.
    HydrodynamicsConfig config;
    InterfaceNodalConfig interfaceConfig;
    std::shared_ptr<const HydroCaseSchedule> schedule;
    DictionaryGeometry geometryInput;
    std::optional<HydroGeometryInput::AxisymmetricColumn> initialColumn;
    std::optional<HydroGeometryInput::ParametricLoop> initialParametric;
    std::optional<InterfaceColumnMode> acceptedColumnMode;
    std::optional<InterfaceAdaptationControls> adaptation;
    std::ofstream adaptationHistory;
    int meshChanges=0;
    bool columnStopped=false;
    DictionaryMesh meshInput;
    int interfaceIndex = -1;
    bool rz = false, openBoundary = false, scalarBoundary = false;
    bool hermiteForceGeometry = false, interfaceForceEnabled = false;
    HydroInterfaceLoadGeometry::Options loadGeometryOptions;
    std::unique_ptr<QuadtreeMesh> meshOwner;
    std::unique_ptr<Hydrodynamics> hydroOwner;
    std::unique_ptr<HydroMomentCase> momentCase;
    std::unique_ptr<NativeLeafDomain> midpointDomain;
    std::unique_ptr<InterfaceNodalState> interfaceState;
    std::unique_ptr<HydroConservativeInterface> conservativeState;
    const ImmersedGeometry<QuadtreeMesh>* geometry() const
    { return interfaceState ? &interfaceState->geometry() : nullptr; }
    const InterfaceBoundaryGeometry::Report& boundaryGeometry() const
    {
        static const InterfaceBoundaryGeometry::Report noInterface;
        return interfaceState ? interfaceState->boundaryGeometry() : noInterface;
    }
    bool axisInterface() const { return boundaryGeometry().axisContact; }
    double initialGeometricVolume = 0.0, integratedNegativeInflow = 0.0;
    CellField<double> couplingForceX, couplingForceY;
    VEdgeField<double> couplingFaceForceX;
    HEdgeField<double> couplingFaceForceY;
    std::vector<double> initialLevelSet;
    HydroSolveOptions options;
    HydroFlowDiagnostics::GeometryMaterialDiagnostics geometryDiagnostic;
    HydroInterfaceForce::Report forceDiagnostic;
    HydroFlowDiagnostics::BoundaryFlowDiagnostics flowDiagnostic;
    std::unique_ptr<HydroOutput::DiagnosticsWriter> writer;
};

HydroCase::Impl::Impl(HydroCaseCoupling* provider) : coupling(provider)
{
    auto input = readHydroCaseInput(coupling != nullptr);
    config = std::move(input.flow);
    interfaceConfig = std::move(input.interface);
    schedule = std::move(input.schedule);
    nextDt=schedule->controls.dt;
    interfaceIndex = HydroGeometryInput::load(config, interfaceConfig, geometryInput, meshInput,
        "./system/geometry.dict", "./system/mesh.dict",&initialColumn,&adaptation,&initialParametric);
    if(schedule->controls.columnAmplitudeLimit && !initialColumn)
        throw std::invalid_argument("Time.Stop columnAmplitude requires geometry AxisymmetricColumn");
    if(interfaceConfig.neckBreakup && schedule->controls.columnAmplitudeLimit)
        throw std::invalid_argument("terminal breakup and small-amplitude stopping are separate run modes");
    if(interfaceConfig.neckBreakup && adaptation
        && interfaceConfig.neckBreakup->maximumLevel!=adaptation->policy.maximumLevel)
        throw std::invalid_argument("neckBreakup maximumLevel must match the actual adaptation policy");
    if(schedule->controls.columnAmplitudeLimit && std::abs(initialColumn->relativeAmplitude)
        /std::sqrt(1+initialColumn->relativeAmplitude*initialColumn->relativeAmplitude/2)>=*schedule->controls.columnAmplitudeLimit)
        throw std::invalid_argument("column amplitude stopping limit must exceed its initial amplitude");
    rz = config.coordinateSystem == "axisymmetricRZ";
    if(config.timeScheme=="cchMomentMidpoint"){
        if(coupling || (adaptation&&!initialParametric) || (!initialParametric&&config.materialModel!="constant") || interfaceIndex>=0 || initialColumn
            || (!initialParametric&&interfaceConfig.advect) || !interfaceConfig.boundary.model.empty() || schedule->controls.stepRetry
            || (!initialParametric&&config.surfaceTension!=0) || config.interfaceNormalTraction!=0 || config.interfaceTangentialTraction!=0)
            throw std::invalid_argument("cchMomentMidpoint requires standalone flow and supported parametric boundary ownership");
        validateHydroBoundaryPatches(config.boundaryPatches,geometryInput.domainMin,geometryInput.domainMax,
            {(geometryInput.domainMax.x-geometryInput.domainMin.x)/meshInput.Nx,(geometryInput.domainMax.y-geometryInput.domainMin.y)/meshInput.Ny},rz,true);
        if(rz&&config.initialVelocity=="vortex"&&meshInput.minCoord.x!=0)
            throw std::invalid_argument("RZ moment vortex requires a domain starting at r=0");
        std::cout<<"[Hydrodynamics] cchMomentMidpoint: "<<(initialParametric?"moving two-phase":"constant material")<<(config.boundaryPatches.empty()?", no-slip enclosure ":", explicit boundary patches ")<<(rz?"RZ":"XY")<<" leaves; "<<meshInput.Nx<<'x'<<meshInput.Ny
            <<" root cells, "<<meshInput.refineBoxes.size()<<" static refinement boxes, dynamic AMR="<<bool(adaptation)<<"; "<<schedule->steps<<" steps, dt="<<schedule->controls.dt<<".\n"
            <<"[Hydrodynamics] Velocity belongs to step end; pressure to solved midpoint. Convection="<<config.momentConvection
            <<". Energy history: output/series/hydro_moment_energy.csv.\n";return;
    }
    if(config.timeScheme=="cchMidpoint" && (coupling || adaptation || config.materialModel!="constant"
        || !config.boundaryPatches.empty() || schedule->controls.stepRetry))
        throw std::invalid_argument("cchMidpoint case requires standalone constant-material fixed uniform closed flow");
    openBoundary = !config.boundaryPatches.empty();
    scalarBoundary = !interfaceConfig.boundary.model.empty();
    if(interfaceConfig.transport=="conservative" && coupling && !coupling->supportsConservative())
        throw std::runtime_error("selected coupling has no joint conservative acceptance protocol");
    if(adaptation && coupling && interfaceConfig.transport=="conservative")
        throw std::runtime_error("conservative electric/flow mesh transfer is not yet jointly connected");
    if (scalarBoundary && coupling)
        throw std::runtime_error("prescribedScalar is currently Hydro-only: wall electric transmission and moving contact charge transport are not implemented");
    validateHydroBoundaryPatches(config.boundaryPatches, geometryInput.domainMin, geometryInput.domainMax,
        Vector2D((geometryInput.domainMax.x-geometryInput.domainMin.x)/meshInput.Nx,
            (geometryInput.domainMax.y-geometryInput.domainMin.y)/meshInput.Ny), rz);
    hermiteForceGeometry = config.interfaceForceGeometry == "hermite";
    loadGeometryOptions.axisymmetric = rz;
    loadGeometryOptions.includeTangentialLoads = coupling && coupling->needsTangentialLoadGeometry();
    if (config.interfaceForceLayer == "regularized")
    {
        loadGeometryOptions.layer = HydroInterfaceLoadGeometry::Layer::Regularized;
        loadGeometryOptions.halfWidth = config.interfaceThickness;
    }
    if (scalarBoundary)
        HydroInterfaceBoundarySupport::configureLoads(interfaceConfig.boundary.fixedContacts, config.boundaryPatches, geometryInput.domainMin,
            geometryInput.domainMax, loadGeometryOptions);
    if (hermiteForceGeometry && (interfaceIndex < 0 || meshInput.Nx < 4 || meshInput.Ny < 4))
        throw std::runtime_error("interfaceForceGeometry hermite requires a named twoPhase interface and at least 4x4 uniform cells");
    interfaceForceEnabled = config.surfaceTension != 0.0 ||
        config.interfaceNormalTraction != 0.0 || config.interfaceTangentialTraction != 0.0;
    if (interfaceForceEnabled && interfaceIndex < 0 && !initialColumn)
        throw std::runtime_error("interface forces require a named twoPhase material interface");
    if (coupling)
        coupling->preflight(config, interfaceConfig, geometryInput, meshInput, interfaceIndex, schedule);
    std::cout << "[Hydrodynamics] " << config.coordinateSystem << (adaptation?", adaptive from root ":", fixed uniform ")
              << meshInput.Nx << "x" << meshInput.Ny << ", materialModel=" << config.materialModel
              << ", steps=" << schedule->steps << ", dt=" << schedule->controls.dt << " (code units)\n";
    if(schedule->controls.stepRetry)
        std::cout<<"[Hydrodynamics] Conservative transactional clock: dt is the initial maximum, halve failed physical steps; minDt="
            <<schedule->controls.stepRetry->minDt<<", maxRetries="<<schedule->controls.stepRetry->maxRetries
            <<", recoverAfter="<<schedule->controls.stepRetry->recoverAfter<<" (0 disables 1.25x recovery)"
            <<", exact target time="<<schedule->controls.endTime<<". Step count above is nominal.\n";
    if(adaptation)
        std::cout<<"[Hydrodynamics] Closed conservative AMR: review both interface histories after each accepted step;"
            " joint mesh/flow/interface publication after the next physical solve. Initial root grid must support step 1."
            " Current-geometry/stencil policy, no certified continuous travel bound; no coarsening.\n";
    if (openBoundary)
        std::cout << "[Hydrodynamics] Open boundary patches: prescribed outward normal velocity (negative=inflow),"
                     " or full pressure traction (-pI+2muD)n=-p_out*n; zero tangential traction."
                     " Unpatched sides retain zero normal velocity/free slip. No solid bodies or AMR.\n";
    else if (coupling)
        std::cout << "[Hydrodynamics] Closed domain: zero normal velocity; optional interface traction plus external coupling force, no solid boundaries, evaporation or AMR.\n";
    else if(adaptation)
        std::cout<<"[Hydrodynamics] Closed adaptive domain: zero normal velocity/free slip, conservative capillary loads.\n";
    else
        std::cout << "[Hydrodynamics] Closed domain: zero normal velocity; optional interface traction, no solid boundaries, electric solve, evaporation or AMR.\n";
    if (interfaceIndex >= 0)
        std::cout << "[Hydrodynamics] Interface binding: '" << config.interfaceLevelSet
                  << "', index=" << interfaceIndex << ", negative/positive density="
                  << config.densityNegative << '/' << config.densityPositive
                  << ", advect=" << interfaceConfig.advect << '\n';
    if(initialColumn)
        std::cout<<"[Hydrodynamics] Analytic column: halfWaves="<<initialColumn->halfWaves<<", R="<<initialColumn->radius
            <<", relativeAmplitude="<<initialColumn->relativeAmplitude
            <<", k="<<initialColumn->halfWaves*std::acos(-1.)/(geometryInput.domainMax.y-geometryInput.domainMin.y)
            <<"; south/north reflection sections, moving endpoint radii, no periodic pairing.\n";
    if (scalarBoundary)
        std::cout << "[Hydrodynamics] interfaceBoundary=prescribedScalar: configured affine inlet phi,"
                     " exact fixed boundary zero vertices, shared open-chain geometry and derivative-only refresh."
                     " No contact-angle law, mechanical support-reaction report or wall electric coupling."
                     " geometric_negative_volume is GL8 zero-curve volume, NOT conservative VOF;"
                     " its boundary-flux budget defect is reported without correction.\n";
    if (hermiteForceGeometry)
        std::cout << "[Hydrodynamics] interfaceForceGeometry=hermite, interfaceForceLayer="
                  << config.interfaceForceLayer << ", halfWidth=" << loadGeometryOptions.halfWidth
                  << ": one t^n public snapshot and selected force-layer H;"
                     " capillarity is -sigma*surfaceCurvature*G(H) at the same load frames used by a context-aware provider."
                     " C1 geometry/local projection only; unsupported topology fails without fallback.\n";
    if (hermiteForceGeometry && rz)
        std::cout << "[Hydrodynamics] RZ Hermite: surfaceCurvature=kappa_meridional+n_r/r,"
                     " physical dual measure=2*pi*r_face*dualArea; off-axis loops or regularized axis-to-axis interfaces."
                     " Axis contact does not support prescribed normal/tangential traction."
                     " Meridional PL fraction is not an axisymmetric volume fraction.\n";
    if(interfaceConfig.transport=="conservative")
        std::cout<<"[Hydrodynamics] Conservative native Interface: C-based density/viscosity and sigma*kappa*G(C);"
            " one frozen solved velocity drives trajectories and phase fluxes, first-order consistent mass/momentum."
            " interfaceThickness is the curvature projection search half-width, not material smoothing.\n";
    if (interfaceForceEnabled && interfaceConfig.transport!="conservative")
        std::cout << "[Hydrodynamics] Interface forcing: surfaceTension=" << config.surfaceTension
                  << ", normalTraction=" << config.interfaceNormalTraction
                  << ", tangentialTraction=" << config.interfaceTangentialTraction << '\n'
                  << "[Hydrodynamics] CSV force_time labels the load at the beginning of each momentum step;"
                     " field frames contain the resulting flow/geometry and omit lagged force fields.\n"
                  << (hermiteForceGeometry
                      ? "[Hydrodynamics] Capillarity uses the shared selected force-layer H and pressure-compatible staggered G(H); H is not transported VOF.\n"
                      : "[Hydrodynamics] Capillarity uses balanced staggered sigma*kappa_face*G(Hcap);"
                        " Hcap is a compact geometric snapshot, not material tanh or transported VOF.\n")
                  << "[Hydrodynamics] Existing interface_force/spread_force/max_interface_force_density CSV columns"
                     " are legacy node/Peskin reference quadrature (including sigma), NOT active sigma loads."
                     " Active sigma face integrals/component maximum are appended separately;"
                     " only prescribed non-sigma traction uses active cell spreading."
                     " Axis-to-axis cases mark the unavailable legacy reference columns empty.\n";
    if(config.timeScheme=="cchMidpoint")std::cout<<"[Hydrodynamics] cchMidpoint: shared midpoint mass/momentum, compatible pressure and energy viscosity."
        " Velocity is at step end; pressure at step midpoint. Native mass, kinetic energy and accumulated viscous/external work: output/series/hydro_energy.csv.\n";
    else std::cout << "[Hydrodynamics] force-pressure imbalance is total t^n cell/face load minus G(p^(n+1))"
                 " on interior velocity faces (component max and combined dual-volume RMS)."
                 " It does not subtract inertial/viscous terms; not a full momentum residual."
                 " Dynamic flow need not make it zero.\n"
              << "[Hydrodynamics] Fixed uniform steps additionally use a coupled correction."
                 " Predictor/projection residuals retain their original meanings; coupled-system and"
                 " full discrete momentum-balance residuals are reported separately.\n";
}

void HydroCase::Impl::appendReport(int step, HydroSolveReport& report,double actualTime,double actualDt)
{
    try
    {
        // Failed modules may contain partial work and reject queries. Preserve
        // the original failure instead of attempting another state inspection.
        const double unavailable = std::numeric_limits<double>::quiet_NaN();
        writer->appendReport(step, actualTime>=0?actualTime:step*schedule->controls.dt,
            actualDt>0?actualDt:schedule->controls.dt, report,
            report.success ? HydroFlowDiagnostics::axisVelocity(*meshOwner, *hydroOwner, rz) : unavailable,
            geometryDiagnostic, forceDiagnostic, flowDiagnostic,
            report.success ? boundaryGeometry().negativeVolume : unavailable,
            initialGeometricVolume, integratedNegativeInflow);
    }
    catch (const std::exception& error)
    {
        if (report.success)
        {
            report.success = false;
            report.stage = "diagnostic-output";
            report.reason = error.what();
        }
        else std::cerr << "[Hydrodynamics] Failure diagnostics unavailable: " << error.what() << '\n';
    }
}

int HydroCase::Impl::initialize()
{
    if (lifecycle != Lifecycle::Prepared)
        throw std::logic_error("Hydro case initialization is single-use");
    lifecycle = Lifecycle::Failed;
    meshOwner = std::make_unique<QuadtreeMesh>(meshInput);
    auto& mesh = *meshOwner;
    if(config.timeScheme=="cchMomentMidpoint"){
        try{
            momentCase=std::make_unique<HydroMomentCase>(meshOwner,config,initialParametric,adaptation?&*adaptation:nullptr);
            // Initial AMR can replace the mesh owner; the old borrow above
            // must not be used for the accepted initial output.
            momentCase->write(*meshOwner,schedule->controls.outputPrefix,schedule->writeStepFields(0));
        }
        catch(const std::exception& error){std::cerr<<"[Hydrodynamics] Moment initialization failed: "<<error.what()<<'\n';return 3;}
        lifecycle=Lifecycle::Ready;return 0;
    }
    hydroOwner = std::make_unique<Hydrodynamics>(mesh);
    auto& hydro = *hydroOwner;
    if(config.timeScheme=="cchMidpoint"){
        FaceSegmentGraph graph;graph.build(mesh);
        midpointDomain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(mesh,graph,rz));
    }
    options.closedDomainProjection = true;
    options.linearResidualTolerance = config.linearResidualTolerance;
    options.completeConservation = config.conservation;
    options.accelerationX = config.accelerationX;
    options.accelerationY = config.accelerationY;
    options.boundaryPatches = config.boundaryPatches;
    writer = std::make_unique<HydroOutput::DiagnosticsWriter>(config, scalarBoundary);
    HydroSolveReport initialization;
    try
    {
        initialization.stage = "geometry-initialize";
        if(initialColumn) {
            conservativeState=std::make_unique<HydroConservativeInterface>(mesh,nullptr,config,interfaceConfig,&*initialColumn,adaptation.has_value());
        } else if (interfaceIndex >= 0)
        {
            auto initialGeometry = std::make_unique<ImmersedGeometry<QuadtreeMesh>>(geometryInput, mesh);
            if (scalarBoundary)
                initialGeometry->setLevelSetBoundaryConditions(mesh, interfaceIndex,
                    HydroInterfaceBoundarySupport::scalarConditions(interfaceConfig.boundary, config.boundaryPatches));
            initialGeometry->initializeLevelSets(mesh);
            if(interfaceConfig.transport=="conservative") {
                conservativeState=std::make_unique<HydroConservativeInterface>(mesh,
                    &initialGeometry->levelSets[interfaceIndex],config,interfaceConfig,nullptr,adaptation.has_value());
            } else {
            InterfaceNodalOptions interfaceOptions;
            interfaceOptions.axisymmetric = rz;
            interfaceOptions.allowAxisContact = supportsAxisInterface(mesh, config, interfaceConfig);
            interfaceOptions.prescribedBoundary = scalarBoundary;
            interfaceOptions.reinitializeEvery = interfaceConfig.reinitializeEvery;
            interfaceOptions.fixedContacts = interfaceConfig.boundary.fixedContacts;
            interfaceState = std::make_unique<InterfaceNodalState>(mesh,
                std::move(initialGeometry), interfaceIndex, std::move(interfaceOptions));
            if (!scalarBoundary && axisInterface() &&
                (config.interfaceNormalTraction != 0.0 || config.interfaceTangentialTraction != 0.0))
                throw std::runtime_error("prescribed normal/tangential traction on an axis-contact interface requires a public face-load implementation");
            initialGeometricVolume = boundaryGeometry().negativeVolume;
            forceDiagnostic.legacyReferenceAvailable = !axisInterface() && !scalarBoundary;
            const auto& levelSet = geometry()->levelSets[interfaceIndex];
            initialLevelSet.assign(levelSet.value.get(), levelSet.value.get() + levelSet.size());
            }
        }
        initialization.stage = "material-initialize";
        if(conservativeState) conservativeState->materials(hydro,0);
        else updateMaterials(mesh, hydro, geometry(), interfaceIndex, config, interfaceConfig.advect, 0.0);
        std::unique_ptr<CellField<double>> initialPressure;
        if(config.initialPressureJump!=0) {
            if(!conservativeState || !config.boundaryPatches.empty())
                throw std::invalid_argument("initial pressure jump requires a closed conservative interface");
            const auto& phase=conservativeState->state().conservativePhase();
            initialPressure=std::make_unique<CellField<double>>(mesh,0.0);
            for(int i=0;i<mesh.cellIndexNum();++i)
                initialPressure->value[i]=config.initialPressureJump*phase.liquidVolumes()[i]/phase.capacities()[i];
        }
        hydro.initializeCompleteVelocity(mesh, HydroStepGeometry::complete(rz),
            config.initialVelocity == "zero" ? HydroInitialVelocity::Zero
                :(config.initialVelocity=="vortex"?HydroInitialVelocity::Vortex:HydroInitialVelocity::Smooth),
            config.velocityAmplitude, initialPressure.get());
        if (scalarBoundary || conservativeState || (midpointDomain && config.initialVelocity=="smooth"))
        {
            // Establish the real t=0 field and its velocity histories before
            // the first old/new RK2 average. This is a zero-time constraint
            // projection; its multiplier is not the physical pressure.
            initialization.stage = "velocity-compatibility";
            const auto velocity = hydro.initializeCompatibleVelocity(mesh, options, rz);
            if (!velocity.success)
                throw std::runtime_error(velocity.stage + ": " + velocity.reason);
            const auto precision = std::cout.precision();
            std::cout << std::setprecision(17)
                      << "[Hydrodynamics] Initial velocity compatibility: time=0, SuperLU="
                      << velocity.usedSuperLU << ", linearResidual=" << velocity.linearRelativeResidual
                      << ", divBefore=" << velocity.divergenceBeforeProjection
                      << ", divAfter=" << velocity.divergenceAfterProjection
                      << ", maxSpeed=" << velocity.maxSpeed
                      << "; current/pre/mid=U0, physical pressure unchanged, no time step consumed.\n";
            std::cout.precision(precision);
        }
        if (coupling)
        {
            initialization.stage = "coupling-initialize";
            if(conservativeState)coupling->initializeConservative(mesh,conservativeState->domain(),conservativeState->state());
            else coupling->initialize(mesh, geometry(), interfaceIndex, config);
            initialization.stage = "material-initialize";
        }
        geometryDiagnostic = collectDiagnostics(mesh, hydro, geometry(), interfaceIndex, initialLevelSet, rz, config.interfaceThickness,
            scalarBoundary ? boundaryGeometry().cutCells : -1);
        initialization.stage = "initial-column-diagnostic";
        if(conservativeState)acceptedColumnMode=conservativeState->columnMode(conservativeState->state());
        if(schedule->controls.columnAmplitudeLimit && (!acceptedColumnMode || !acceptedColumnMode->singleRadiusGraph
            || std::abs(acceptedColumnMode->relativeAmplitude)>=*schedule->controls.columnAmplitudeLimit))
            throw std::invalid_argument("initial discrete column must resolve a radial graph below its amplitude limit");
        initialization.stage = "initial-field-output";
        if(midpointDomain)writer->appendMidpoint(0,0,0,0,initialization,HydroFlowDiagnostics::nativeConservation(mesh,hydro,*midpointDomain));
        if(conservativeState)
            conservativeState->write(hydro,0,0,schedule->controls.outputPrefix,schedule->writeStepFields(0),0,nullptr,nullptr,
                acceptedColumnMode?&*acceptedColumnMode:nullptr);
        else if (schedule->writeStepFields(0))
            HydroOutput::writeFields(mesh, hydro, geometry(), interfaceIndex, config, 0, 0 * schedule->controls.dt,
                schedule->controls.outputPrefix, rz && mesh.minCoord.x == 0.0, scalarBoundary,0);
    }
    catch (const std::exception& error)
    {
        initialization.reason = error.what();
        appendReport(0, initialization);
        std::cerr << "[Hydrodynamics] ERROR: initialization failed at " << initialization.stage
                  << ": " << initialization.reason << '\n';
        return 3;
    }

    lifecycle = Lifecycle::Ready;
    return 0;
}

int HydroCase::Impl::advance(int step, double forceTime)
{
    if(schedule->controls.stepRetry) throw std::logic_error("transactional Hydro clock requires advanceNext()");
    if (lifecycle != Lifecycle::Ready || columnStopped || step != completedSteps + 1 ||
        step > schedule->steps || forceTime != (step - 1) * schedule->controls.dt)
        throw std::logic_error("Hydro case requires the next scheduled step on a ready state");
    if(conservativeState) return advanceConservative(step,forceTime,schedule->controls.dt);
    if(momentCase)return advanceMoment(step);
    if(midpointDomain)return advanceMidpoint(step,forceTime);
    // A failed coupling/geometry/output phase must not permit another call just
    // because the NS substep succeeded. No rollback or retry is implied.
    lifecycle = Lifecycle::Failed;
    auto& mesh = *meshOwner;
    auto& hydro = *hydroOwner;
    HydroSolveReport report;
    forceDiagnostic = HydroInterfaceForce::Report{};
    forceDiagnostic.legacyReferenceAvailable = !axisInterface() && !scalarBoundary;
    forceDiagnostic.forceTime = forceTime;
    try
    {
        report.stage = "material-update";
        updateMaterials(mesh, hydro, geometry(), interfaceIndex, config, interfaceConfig.advect, forceDiagnostic.forceTime);
        std::unique_ptr<HydroInterfaceLoadGeometry> loadGeometry;
        if (hermiteForceGeometry)
        {
            report.stage = "interface-load-geometry";
            loadGeometry = std::make_unique<HydroInterfaceLoadGeometry>(mesh, geometry()->levelSets[interfaceIndex], loadGeometryOptions);
        }
        HydroInterfaceForce::Assembly force(mesh, geometry(), interfaceIndex, rz, axisInterface(),
            {config.surfaceTension, config.interfaceNormalTraction, config.interfaceTangentialTraction,
                scalarBoundary, supportsAxisInterface(mesh, config, interfaceConfig)},
            forceDiagnostic, report.stage, loadGeometry.get());
        if (coupling)
        {
            report.stage = "coupling-force";
            couplingForceX.initialize(mesh, 0.0);
            couplingForceY.initialize(mesh, 0.0);
            // The provider sees exactly the geometry used for this
            // momentum load, and writes any electric frames now,
            // before geometry advection changes its time layer.
            if (loadGeometry)
                coupling->computeForceWithGeometry(mesh, geometry(), interfaceIndex,
                    step, forceDiagnostic.forceTime, schedule->controls.dt, *loadGeometry, couplingForceX, couplingForceY);
            else
                coupling->computeForce(mesh, geometry(), interfaceIndex,
                    step, forceDiagnostic.forceTime, schedule->controls.dt, couplingForceX, couplingForceY);
            force.addCouplingCells(couplingForceX, couplingForceY);
            couplingFaceForceX.initialize(mesh, 0.0);
            couplingFaceForceY.initialize(mesh, 0.0);
            // A face provider reads the solve just performed above.
            // There is exactly one computeForce/electric solve per step.
            const bool suppliedFaceForce = loadGeometry
                ? coupling->computeFaceForceWithGeometry(mesh, geometry(), interfaceIndex,
                    step, forceDiagnostic.forceTime, schedule->controls.dt, *loadGeometry, couplingFaceForceX, couplingFaceForceY)
                : coupling->computeFaceForce(mesh, geometry(), interfaceIndex,
                    step, forceDiagnostic.forceTime, schedule->controls.dt, couplingFaceForceX, couplingFaceForceY);
            force.addCouplingFaces(couplingFaceForceX, couplingFaceForceY, suppliedFaceForce);
            // Keep forceDiagnostic capillary/prescribed-traction only;
            // the provider owns separately time-labelled diagnostics.
        }
        if (loadGeometry)
        {
            report.stage = "interface-load-source-check";
            loadGeometry->validateSource(mesh, geometry()->levelSets[interfaceIndex]);
        }
        auto load = force.publish(forceDiagnostic.forceTime);
        if (load.cellX || load.faceX)
            load.source = hydro.captureForceSource(forceDiagnostic.forceTime);
        auto hydroCandidate = hydro.prepareStep(mesh, {HydroStepGeometry::complete(rz),
            {step, forceDiagnostic.forceTime, schedule->controls.dt}, options, load});
        report = hydroCandidate.report();
        if (report.success)
        {
            if (openBoundary)
            {
                report.stage = "boundary-flow-diagnostics";
                flowDiagnostic = boundaryFlow(mesh, hydro, rz, config.conservation.boundaryFluxBalance,
                    config.materialModel == "constant");
            }
            if (scalarBoundary)
                integratedNegativeInflow -= schedule->controls.dt*hydro.completePhaseOutwardFlux(
                    mesh, loadGeometry->snapshot(), rz);
            report.stage = "force-pressure-diagnostics";
            collectForcePressureImbalance(mesh, hydro, force.faceX(), force.faceY(), rz, forceDiagnostic);
            // All t^n consumers and load output have completed. Never
            // retain this snapshot across the following phi update.
            loadGeometry.reset();
            std::unique_ptr<InterfaceNodalMaterialMap> materialMap;
            if (interfaceState && interfaceConfig.advect)
            {
                report.stage = "geometry-advection";
                const bool velocityAxis = axisInterface() || (openBoundary && rz && mesh.minCoord.x == 0.0);
                auto [vx, vy] = hydro.vertexVelocity(mesh, HydroVelocityTimeLevel::Current, openBoundary, velocityAxis);
                auto [previousVx, previousVy] = hydro.vertexVelocity(mesh,
                    HydroVelocityTimeLevel::Previous, openBoundary, velocityAxis);
                // Hydro retained U^n in its Previous layer. The geometry tracer
                // averages these with U^(n+1), matching this used dt.
                LevelSetAdvectionVelocity advectionVelocity;
                if (scalarBoundary)
                {
                    advectionVelocity.current = hydro.boundaryAdvectionVelocity(mesh, HydroVelocityTimeLevel::Current,
                        config.boundaryPatches, velocityAxis);
                    advectionVelocity.previous = hydro.boundaryAdvectionVelocity(mesh, HydroVelocityTimeLevel::Previous,
                        config.boundaryPatches, velocityAxis);
                }
                InterfaceNodalVelocityLayers velocity{std::move(vx), std::move(vy),
                    std::move(previousVx), std::move(previousVy), std::move(advectionVelocity)};
                try
                {
                    materialMap = interfaceState->advance(mesh, step, forceTime, schedule->controls.dt,
                        std::move(velocity), coupling && coupling->needsInterfaceMaterialMap());
                }
                catch (...)
                {
                    report.stage = interfaceState->stage();
                    throw;
                }
                report.stage = "material-output-update";
                updateMaterials(mesh, hydro, geometry(), interfaceIndex, config, interfaceConfig.advect, step * schedule->controls.dt);
            }
            geometryDiagnostic = collectDiagnostics(mesh, hydro, geometry(), interfaceIndex, initialLevelSet, rz, config.interfaceThickness,
                scalarBoundary ? boundaryGeometry().cutCells : -1);
            if (coupling)
            {
                report.stage = "coupling-step-complete";
                coupling->completeStepWithVelocity(mesh, geometry(), interfaceIndex, materialMap.get(),
                    hydro,
                    step, forceDiagnostic.forceTime, schedule->controls.dt);
            }
            // The sole App acceptance point: Hydro, Interface, material and
            // electric substeps/checks have all completed. IO cannot roll back
            // those physical states. No new-time field was published above.
            hydroCandidate.accept(mesh);
            acceptedStep = step;
            acceptedTime=step*schedule->controls.dt;
            if (coupling)
            {
                report.stage = "coupling-output";
                coupling->publishCompletedStep(step, forceDiagnostic.forceTime, schedule->controls.dt);
            }
            if (schedule->writeStepFields(step))
            {
                report.stage = "field-output";
                HydroOutput::writeFields(mesh, hydro, geometry(), interfaceIndex, config, step, step * schedule->controls.dt,
                schedule->controls.outputPrefix, rz && mesh.minCoord.x == 0.0, scalarBoundary);
            }
            report.stage = "complete";
        }
    }
    catch (const std::exception& error)
    {
        report.success = false;
        report.reason = error.what();
    }
    appendReport(step, report);
    if (!report.success)
    {
        std::cerr << "[Hydrodynamics] ERROR: step " << step << " failed at " << report.stage
                  << ": " << report.reason << "; physical steps accepted through " << acceptedStep
                  << ", output completed through " << completedSteps << '\n';
        return 3;
    }
    if (step == 1)
        std::cout << "[Hydrodynamics] 首步实际求解器：动量预测 momentumBackend="
                  << (report.momentumUsedSuperLU ? "SuperLU" : "unverified")
                  << "，压力投影 pressureBackend=" << (report.pressureUsedSuperLU ? "SuperLU" : "unverified")
                  << "，耦合修正 coupledBackend=" << (report.coupledUsedSuperLU ? "SuperLU"
                      : (report.coupledCorrectionUsed ? "unverified" : "not-used"))
                  << "；逐步数值诊断保存在 ./output/series/hydro_diagnostics.csv。\n";
    if (step % schedule->controls.logEverySteps == 0 || step == schedule->steps)
    {
        std::cout << "[Hydrodynamics] 流动推进 step=" << step << '/' << schedule->steps
                  << " | time=" << step * schedule->controls.dt << " (code units)"
                  << " | 最大流速=" << report.maxSpeed
                  << " | 连续性误差(max|divU|)=" << report.divergenceAfterProjection
                  << " | 完整动量相对残差=" << report.momentumBalanceRelativeResidual
                  << " | 界面变化(max|delta_phi|)=" << geometryDiagnostic.interfaceMaxChange
                  << " | 活跃毛细面力密度=" << forceDiagnostic.maxActiveCapillaryFaceDensity << '\n';
        if (scalarBoundary)
        {
            const auto savedPrecision = std::cout.precision();
            std::cout << std::setprecision(17) << "[Hydrodynamics] 显式标量入口体积诊断 step=" << step
                      << '/' << schedule->steps << " (code units)"
                      << " | 负相几何体积=" << boundaryGeometry().negativeVolume
                      << " | 初始体积=" << initialGeometricVolume
                      << " | 累计净流入=" << integratedNegativeInflow
                      << " | 体积收支缺陷="
                      << boundaryGeometry().negativeVolume-initialGeometricVolume-integratedNegativeInflow
                      << " | 这是水平集输运，不是守恒VOF；此行不代表体积验收通过。\n";
            std::cout.precision(savedPrecision);
        }
    }

    completedSteps = step;outputTime=acceptedTime;
    lifecycle = Lifecycle::Ready;
    return 0;
}

int HydroCase::Impl::advanceMoment(int step)
{
    lifecycle=Lifecycle::Failed;const char* stage="moment-advance";
    try{
        stage="moment-advance";
        // Continue the App's accepted clock exactly. Recomputing (step-1)*dt
        // can differ from the preceding begin+dt by one ulp; Interface and
        // Hydro must consume the same interval identity, not two formulas.
        const auto report=momentCase->advance(meshOwner,step,acceptedTime,schedule->controls.dt);
        acceptedStep=momentCase->state().acceptedSteps;acceptedTime=momentCase->state().time;
        stage="moment-output";momentCase->write(*meshOwner,schedule->controls.outputPrefix,schedule->writeStepFields(step),&report);
        completedSteps=acceptedStep;outputTime=acceptedTime;lifecycle=Lifecycle::Ready;
        if(step==1 || step%schedule->controls.logEverySteps==0 || step==schedule->steps)
            std::cout<<std::setprecision(17)<<"[Hydrodynamics] Moment accepted "<<step<<'/'<<schedule->steps<<" time="<<acceptedTime
                <<" pressure_time="<<*momentCase->state().stageTime<<" kinetic="<<report.kineticAfter<<" energy_defect="<<report.energyDefect
                <<" iterations="<<report.nonlinearIterations<<" residual="<<report.equationResidual<<'\n';return 0;
    }catch(const std::exception& error){std::cerr<<"[Hydrodynamics] "<<stage<<" failed: "<<error.what()<<"; accepted="<<acceptedStep<<", output="<<completedSteps<<'\n';return 3;}
}

int HydroCase::Impl::advanceMidpoint(int step,double begin)
{
    lifecycle=Lifecycle::Failed;auto& mesh=*meshOwner;auto& hydro=*hydroOwner;
    const double dt=schedule->controls.dt,end=step*dt;HydroSolveReport report;
    forceDiagnostic=HydroInterfaceForce::Report{};forceDiagnostic.legacyReferenceAvailable=false;forceDiagnostic.forceTime=begin+.5*dt;
    try {
        report.stage="midpoint-advance";
        HydroForceDensity force;force.evaluationTime=forceDiagnostic.forceTime;force.timeLevel=HydroForceTimeLevel::StepMidpoint;
        const auto result=Quantum::Hydrodynamics::advanceSinglePhaseMidpointStep(mesh,*midpointDomain,hydro,
            {step,begin,dt,HydroTimePolicy::ConservativeOneStep},options,force,config.midpoint);
        report=result.hydro;acceptedStep=step;acceptedTime=end;
        report.stage="midpoint-diagnostics";
        geometryDiagnostic=collectDiagnostics(mesh,hydro,nullptr,-1,initialLevelSet,rz,0);
        writer->appendMidpoint(step,end,dt,forceDiagnostic.forceTime,report,HydroFlowDiagnostics::nativeConservation(mesh,hydro,*midpointDomain),
            result.iterations,result.midpointVelocityDefect);
        if(schedule->writeStepFields(step)){
            report.stage="field-output";
            HydroOutput::writeFields(mesh,hydro,nullptr,-1,config,step,end,schedule->controls.outputPrefix,
                rz && mesh.minCoord.x==0,false,forceDiagnostic.forceTime);
        }
        report.stage="complete";
    } catch(const std::exception& error){report.success=false;report.reason=error.what();}
    appendReport(step,report,end,dt);
    if(!report.success){std::cerr<<"[Hydrodynamics] Midpoint step "<<step<<" failed at "<<report.stage<<": "<<report.reason
        <<"; accepted="<<acceptedStep<<", output="<<completedSteps<<'\n';return 3;}
    completedSteps=step;outputTime=end;lifecycle=Lifecycle::Ready;
    if(step==1 || step%schedule->controls.logEverySteps==0 || step==schedule->steps)
        std::cout<<"[Hydrodynamics] Midpoint accepted "<<step<<'/'<<schedule->steps<<" time="<<end<<" pressure_time="<<forceDiagnostic.forceTime
            <<" max_speed="<<report.maxSpeed<<" divU="<<report.divergenceAfterProjection<<" mass_closure="<<report.conservativeMassRelativeClosure<<'\n';
    return 0;
}

bool HydroCase::Impl::hasNextStep() const
{
    if(lifecycle!=Lifecycle::Ready) throw std::logic_error("Hydro clock requires a ready case");
    return !columnStopped && (schedule->controls.stepRetry?acceptedTime<schedule->controls.endTime:completedSteps<schedule->steps);
}

int HydroCase::Impl::advanceNext()
{
    if(!hasNextStep()) throw std::logic_error("Hydro clock has reached its end time");
    if(!schedule->controls.stepRetry) return advance(completedSteps+1,completedSteps*schedule->controls.dt);
    if(!conservativeState || (coupling && !coupling->supportsConservative()))
        throw std::logic_error("transactional clock requires conservative Hydro and a joint coupling protocol");
    return advanceConservative(acceptedStep+1,acceptedTime,nextDt);
}

void HydroCase::Impl::writeAdaptation(int step,int attempt,double time,int sourceCells,const char* status,
    const HydroConservativeInterface::MeshCandidate* candidate,const HydroSolveReport& report)
{
    if(!adaptationHistory.is_open()) {
        adaptationHistory.open(OutputLayout::series()/"hydro_adaptation.csv");
        adaptationHistory<<"step,attempt,source_time,status,source_cells,target_cells,required_leaf_count,planned_leaf_count,interface_level,"
            "current_momentum_relative_closure,previous_momentum_relative_closure,current_remap_divergence,previous_remap_divergence,"
            "current_energy_change,previous_energy_change,pressure_integral_change,current_geometry_inventory_l1,previous_geometry_inventory_l1,"
            "physical_divergence,stage,reason";
        if(interfaceConfig.neckBreakup)adaptationHistory<<",terminal_curvature_requests,maximum_terminal_curvature_cells";
        adaptationHistory<<'\n';
    }
    const bool accepted=std::string(status)=="accepted";
    const bool terminal=std::string(status)=="terminal_model";
    adaptationHistory<<std::setprecision(17)<<step<<','<<attempt<<','<<time<<','<<status<<','<<sourceCells<<','
        <<(accepted || terminal?meshOwner->cellIndexNum():candidate && candidate->mesh?candidate->mesh->cellIndexNum():0);
    if(candidate && (accepted || candidate->hydro)) {
        const auto& f=candidate->flow;
        adaptationHistory<<','<<candidate->plan.requiredLeafCount<<','<<candidate->plan.plannedLeafCount<<','<<candidate->plan.interfaceLevel
            <<','<<f.current.momentumRelativeClosure<<','<<f.previous.momentumRelativeClosure
            <<','<<f.current.divergenceAfter<<','<<f.previous.divergenceAfter
            <<','<<f.current.kineticEnergy-f.current.unprojected.sourceKineticEnergy[0]-f.current.unprojected.sourceKineticEnergy[1]
            <<','<<f.previous.kineticEnergy-f.previous.unprojected.sourceKineticEnergy[0]-f.previous.unprojected.sourceKineticEnergy[1]
            <<','<<f.targetPressureIntegral-f.sourcePressureIntegral
            <<','<<candidate->currentGeometry.geometryInventoryL1<<','<<candidate->previousGeometry.geometryInventoryL1;
    } else for(int i=0;i<12;++i)adaptationHistory<<",nan";
    adaptationHistory<<','<<(accepted || terminal?report.divergenceAfterProjection:std::numeric_limits<double>::quiet_NaN())
        <<','<<report.stage<<",\"";
    for(char ch:report.reason) {if(ch=='"')adaptationHistory<<'"';adaptationHistory<<ch;}
    adaptationHistory<<'"';
    if(interfaceConfig.neckBreakup)adaptationHistory<<','<<(candidate?candidate->plan.terminalCurvatureRequests:0)
        <<','<<(candidate?candidate->plan.maxTerminalCurvatureCells:0);
    adaptationHistory<<'\n';adaptationHistory.flush();
    if(!adaptationHistory)throw std::runtime_error("cannot write conservative Hydro AMR history");
}

int HydroCase::Impl::advanceConservative(int step,double forceTime,double proposedDt)
{
    lifecycle=Lifecycle::Failed;
    const auto& retryPolicy=schedule->controls.stepRetry;
    for(int attempt=0;;++attempt) {
    HydroSolveReport report;forceDiagnostic=HydroInterfaceForce::Report{};
    forceDiagnostic.legacyReferenceAvailable=false;forceDiagnostic.forceTime=forceTime;
    double end=step*schedule->controls.dt,dt=proposedDt;
    if(retryPolicy) {
        const double target=schedule->controls.endTime;
        end=std::min(target,forceTime+proposedDt);
        if(target-end<=8*std::numeric_limits<double>::epsilon()*target) end=target;
        dt=end-forceTime;
        if(!(dt>0)) throw std::runtime_error("Hydro retry step is below clock resolution");
    }
    const InterfaceStepTime time{step,forceTime,end,dt};
    bool writeFields=retryPolicy?schedule->writeVariableStepFields(step,time.begin,time.end):schedule->writeStepFields(step);
    auto stageBegin=std::chrono::steady_clock::now();
    std::array<double,7> stageSeconds{};
    double electricSeconds=0;
    InterfaceConservativeGeometry::AdvanceReport interfaceReport;
    std::optional<HydroConservativeInterface::MeshCandidate> meshCandidate;
    const int sourceCells=meshOwner->cellIndexNum();
    bool adaptationReviewed=false,adaptationPublished=false;
    const auto timeStage=[&](int index) {
        if(!schedule->controls.detailedLog)return;
        const auto end=std::chrono::steady_clock::now();
        stageSeconds[index]=std::chrono::duration<double>(end-stageBegin).count();stageBegin=end;
    };
    try {
        if(adaptation && acceptedStep>0) {
            report.stage="conservative-mesh-preparation";adaptationReviewed=true;
            meshCandidate.emplace(conservativeState->prepareAdaptation(*hydroOwner,*adaptation));
        }
        const bool changing=meshCandidate && meshCandidate->mesh;
        auto& mesh=changing?*meshCandidate->mesh:*meshOwner;
        auto& hydro=changing?*meshCandidate->hydro:*hydroOwner;
        auto* activeInterface=changing?meshCandidate->adapter.get():conservativeState.get();
        auto transaction=hydro.beginStepTransaction(mesh);
        report.stage="conservative-material-source";activeInterface->materials(hydro,time.begin);
        timeStage(0);
        report.stage="conservative-interface-motion";
        HydroConservativeInterface::MaterialMotion materialMotion;
        auto next=activeInterface->prepare(hydro,time,schedule->controls.detailedLog?&interfaceReport:nullptr,
            coupling?&materialMotion:nullptr);
        const auto& previousEvents=activeInterface->state().topology().events();
        const bool newEvent=next->topology().events().size()>previousEvents.size();
        const bool afterEvent=!previousEvents.empty() && previousEvents.back().end==time.begin;
        if(schedule->controls.writeFields && newEvent) {
            report.stage="conservative-event-source-output";
            // During AMR, activeInterface is still an unpublished remap.
            // The before-frame must use the original accepted owners.
            conservativeState->writeEventSourceFrame(*hydroOwner,schedule->controls.outputPrefix);
            if(coupling)coupling->writeConservativeEventSource(acceptedStep,acceptedTime);
            report.stage="conservative-interface-motion";
        }
        writeFields|=schedule->controls.writeFields && (newEvent || afterEvent);
        const auto columnMode=activeInterface->columnMode(*next);
        const auto& amplitudeLimit=schedule->controls.columnAmplitudeLimit;
        if(amplitudeLimit && (!columnMode || !columnMode->singleRadiusGraph))
            throw std::runtime_error("column amplitude stopping requires a resolved single radial graph");
        const bool stopCandidate=amplitudeLimit && std::abs(columnMode->relativeAmplitude)>=*amplitudeLimit;
        writeFields|=stopCandidate && schedule->controls.writeFields && schedule->controls.writeFinal;
        timeStage(1);
        // Drift first with the old solved velocity, then kick momentum with
        // the new interface's capillary force. Using the old interface in
        // both updates makes even an inviscid linear capillary oscillator
        // a forward-Euler pair with amplification modulus greater than one.
        forceDiagnostic.forceTime=time.end;
        if(next->conservativeGeometry().terminalInventoryHandoff)
            std::cout<<"[Interface terminal synchronization] step="<<step
                <<" inventory_handoff=1 gas_strip_cells="<<next->conservativeGeometry().terminalGasStripCells<<'\n';
        if(next->conservativeGeometry().interiorNeckSplit || next->conservativeGeometry().interiorGapProtected)
            std::cout<<"[Interface interior neck] step="<<step
                <<" split="<<next->conservativeGeometry().interiorNeckSplit
                <<" gap_protected="<<next->conservativeGeometry().interiorGapProtected
                <<" liquids="<<next->topology().liquids().size()<<'\n';
        report.stage="conservative-capillary-load";
        Quantum::Hydrodynamics::ConservativeCapillaryReport capillary;
        std::optional<Quantum::Hydrodynamics::HydroConservativePressure> pressure;
        if(mesh.cellIndexNum()!=meshInput.Nx*meshInput.Ny)
            pressure.emplace(Quantum::Hydrodynamics::HydroConservativePressure::build(mesh,activeInterface->domain()));
        const bool recordCapillary=writeFields || !next->eventCurvatureCells().empty();
        std::shared_ptr<const Quantum::Hydrodynamics::ConservativeTractionGeometry> tractionGeometry;
        auto load=pressure?Quantum::Hydrodynamics::buildNativeConservativeCapillaryLoad(mesh,activeInterface->domain(),
            *next,*pressure,config.surfaceTension,config.interfaceThickness,recordCapillary?&capillary:nullptr,coupling?&tractionGeometry:nullptr)
            :Quantum::Hydrodynamics::buildUniformConservativeCapillaryLoad(mesh,activeInterface->domain(),
                *next,config.surfaceTension,config.interfaceThickness,recordCapillary?&capillary:nullptr,coupling?&tractionGeometry:nullptr);
        forceDiagnostic.maxActiveCapillaryFaceDensity=load.maxForce;
        forceDiagnostic.activeCapillaryFaceX=load.totalForceX;forceDiagnostic.activeCapillaryFaceY=load.totalForceY;
        std::unique_ptr<HydroConservativeCouplingCandidate> coupledCandidate;
        if(coupling) {
            report.stage="conservative-coupling-preparation";
            const auto electricBegin=std::chrono::steady_clock::now();
            if(!materialMotion.physical)throw std::logic_error("missing physical surface material transfer");
            coupledCandidate=coupling->prepareConservative(mesh,activeInterface->domain(),activeInterface->state(),*next,
                *materialMotion.physical,materialMotion.event?&*materialMotion.event:nullptr,time,pressure?&*pressure:nullptr,*tractionGeometry);
            if(!coupledCandidate)throw std::logic_error("missing conservative coupling candidate");
            const auto& electric=coupledCandidate->load();
            for(int i=0;i<mesh.vEdgeIndexNum();++i)load.faceForceX.value[i]+=electric.faceForceX.value[i];
            for(int i=0;i<mesh.hEdgeIndexNum();++i)load.faceForceY.value[i]+=electric.faceForceY.value[i];
            electricSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-electricBegin).count();
        }
        timeStage(2);
        stageSeconds[2]-=electricSeconds;
        report.stage="conservative-mass-momentum";
        const auto transport=Quantum::Hydrodynamics::HydroConservativeTransport::build(mesh,activeInterface->domain(),
            activeInterface->state().conservativePhase(),next->conservativePhase(),time,config.densityNegative,config.densityPositive);
        const auto source=hydro.captureForceSource(time.end,HydroForceEvolution::TimeDependent,next->geometry());
        const HydroForceDensity force{nullptr,nullptr,&load.faceForceX,&load.faceForceY,time.end,source,HydroForceTimeLevel::StepEnd};
        timeStage(3);
        auto candidate=hydro.prepareStep(mesh,{HydroStepGeometry::complete(rz),
            {step,time.begin,time.dt,retryPolicy?HydroTimePolicy::ConservativeOneStep:HydroTimePolicy::Fixed},options,force,
            next->geometry().get(),&transport});
        timeStage(4);
        report=candidate.report();
        if(report.success) {
            report.stage="conservative-acceptance-diagnostics";
            flowDiagnostic=boundaryFlow(mesh,hydro,rz,config.conservation.boundaryFluxBalance,false);
            if(pressure)collectForcePressureImbalance(mesh,hydro,load.faceForceX,load.faceForceY,activeInterface->domain(),*pressure,forceDiagnostic);
            else collectForcePressureImbalance(mesh,hydro,load.faceForceX,load.faceForceY,rz,forceDiagnostic);
            double scalarChange=0;
            for(int v=0;v<mesh.vertIndexNum();++v) scalarChange=std::max(scalarChange,
                std::abs(next->displayValue(v)-activeInterface->state().displayValue(v)));
            report.stage="conservative-target-material";activeInterface->materials(hydro,*next,time.end);
            geometryDiagnostic=collectDiagnostics(mesh,hydro,nullptr,-1,{},rz,config.interfaceThickness);
            geometryDiagnostic.negativeVolume=static_cast<double>(next->conservativePhase().report().liquidAfter);
            geometryDiagnostic.interfaceMaxChange=scalarChange;
            geometryDiagnostic.levelSetMin=std::numeric_limits<double>::infinity();
            geometryDiagnostic.levelSetMax=-std::numeric_limits<double>::infinity();
            for(int v=0;v<mesh.vertIndexNum();++v) {
                const double value=next->displayValue(v);
                geometryDiagnostic.levelSetMin=std::min(geometryDiagnostic.levelSetMin,value);
                geometryDiagnostic.levelSetMax=std::max(geometryDiagnostic.levelSetMax,value);
            }
            const auto& phase=next->conservativePhase();
            for(int c:activeInterface->domain().cellIndexByRow())
                if(phase.liquidVolumes()[c]>0 && phase.liquidVolumes()[c]<phase.capacities()[c]) ++geometryDiagnostic.cutCells;
            // All fallible physical preparation precedes joint publication.
            // If acceptance validation throws, the transaction restores Hydro
            // while the old Interface handle is still authoritative. Output
            // follows commit and never rolls back an accepted physical step.
            candidate.accept(mesh);transaction.commit();
            activeInterface->accept(std::move(next));
            if(coupledCandidate)coupledCandidate->commit();
            if(changing) {
                activeInterface->inheritOutput(*conservativeState);
                meshOwner.swap(meshCandidate->mesh);hydroOwner.swap(meshCandidate->hydro);
                conservativeState.swap(meshCandidate->adapter);
                ++meshChanges;adaptationPublished=true;
            }
            acceptedStep=step;acceptedTime=time.end;
            acceptedColumnMode=columnMode;columnStopped=stopCandidate;
            report.stage="conservative-accepted-output";
            timeStage(5);
            if(adaptationPublished)writeAdaptation(step,attempt,time.begin,sourceCells,"accepted",&*meshCandidate,report);
            else if(meshCandidate && meshCandidate->plan.terminalCurvatureRequests)
                writeAdaptation(step,attempt,time.begin,sourceCells,"terminal_model",&*meshCandidate,report);
            conservativeState->write(hydro,step,time.end,schedule->controls.outputPrefix,writeFields,time.end,&report,&capillary,
                acceptedColumnMode?&*acceptedColumnMode:nullptr);
            if(coupling)coupling->publishCompletedStep(step,time.begin,time.dt);
            timeStage(6);
            report.stage="complete";
        }
    } catch(const std::exception& error) {report.success=false;report.reason=error.what();}
    if(!report.success && adaptationReviewed && !adaptationPublished
        && (!meshCandidate || meshCandidate->mesh)) {
        try {writeAdaptation(step,attempt,time.begin,sourceCells,"rejected",meshCandidate?&*meshCandidate:nullptr,report);}
        catch(const std::exception& error) {std::cerr<<"[Hydrodynamics] AMR rejection output unavailable: "<<error.what()<<'\n';}
    }
    if(schedule->controls.detailedLog) {
        const char* names[]={"materials","interface","capillary","mass_flux","NS","acceptance","output"};
        std::cout<<"[Hydrodynamics timing] step="<<step;
        for(int i=0;i<7;++i)std::cout<<' '<<names[i]<<"_seconds="<<stageSeconds[i];
        if(coupling)std::cout<<" electric_seconds="<<electricSeconds;
        std::cout<<std::endl;
        const char* interfaceNames[]={"donors","flux","inventory","predictor","synchronization","surface_transfer"};
        std::cout<<"[Interface timing] step="<<step;
        for(int i=0;i<6;++i)std::cout<<' '<<interfaceNames[i]<<"_seconds="<<interfaceReport.stageSeconds[i];
        std::cout<<" constant_phase_chords="<<interfaceReport.spaceTime.constantPhaseChords;
        std::cout<<" primitive_requests="<<interfaceReport.spaceTime.primitiveRequests
            <<" primitive_evaluations="<<interfaceReport.spaceTime.primitiveEvaluations;
        std::cout<<" arrangement_vertex_requests="<<interfaceReport.spaceTime.arrangementVertexRequests
            <<" arrangement_vertex_builds="<<interfaceReport.spaceTime.arrangementVertexBuilds
            <<" arrangement_edge_root_queries="<<interfaceReport.spaceTime.arrangementEdgeRootQueries;
        std::cout<<" max_spatial_depth="<<interfaceReport.spaceTime.maxAcceptedSpaceDepth
            <<" max_spatial_leaves="<<interfaceReport.spaceTime.maxAcceptedSpaceLeaves;
        std::cout<<std::endl;
    }
    if(!report.success && acceptedStep<step && report.stage!="conservative-event-source-output"
        && retryPolicy && attempt<retryPolicy->maxRetries
        && time.dt*.5>=retryPolicy->minDt) {
        const auto restored=hydroOwner->state().progress;
        if(restored.stopped || restored.pendingAcceptance || restored.acceptedSteps!=acceptedStep
            || restored.completedSteps!=acceptedStep || restored.endTime!=acceptedTime
            || conservativeState->state().progress().endTime!=acceptedTime)
            throw std::runtime_error("Hydro failed step was not jointly restored; retry refused");
        proposedDt=time.dt*.5;
        std::cout<<std::setprecision(17)<<"[Hydrodynamics retry] step="<<step<<" begin="<<time.begin
            <<" rejected_dt="<<time.dt<<" next_dt="<<proposedDt<<" attempt="<<attempt+1
            <<" stage="<<report.stage<<" reason="<<report.reason<<std::endl;
        continue;
    }
    appendReport(step,report,time.end,time.dt);
    if(!report.success) {
        std::cerr<<"[Hydrodynamics] Conservative step "<<step<<" failed at "<<report.stage<<": "<<report.reason
            <<"; accepted="<<acceptedStep<<", output="<<completedSteps<<'\n';return 3;
    }
    completedSteps=step;outputTime=acceptedTime;lifecycle=Lifecycle::Ready;
    if(retryPolicy) {
        nextDt=time.dt;
        if(retryPolicy->recoverAfter>0) {
            retryFreeSteps=attempt==0?retryFreeSteps+1:0;
            if(retryFreeSteps>=retryPolicy->recoverAfter) {
                retryFreeSteps=0;
                const double maximum=schedule->controls.dt;
                // A recovered step is only a proposal. It goes through the
                // same joint preparation, CFL/accuracy gates and rollback;
                // neither a rejected trial nor an output failure earns credit.
                if(!columnStopped && time.end<schedule->controls.endTime
                    && time.dt<maximum*(1-64*std::numeric_limits<double>::epsilon())) {
                    nextDt=time.dt<maximum/1.25?1.25*time.dt:maximum;
                    std::cout<<std::setprecision(17)<<"[Hydrodynamics step recovery] step="<<step
                        <<" time="<<time.end<<" accepted_dt="<<time.dt<<" proposed_dt="<<nextDt
                        <<" retry_free_steps="<<retryPolicy->recoverAfter<<std::endl;
                }
            }
        }
    }
    if(columnStopped || step%schedule->controls.logEverySteps==0 || (retryPolicy?time.end==schedule->controls.endTime:step==schedule->steps)) {
        std::cout<<"[Hydrodynamics] Conservative NS accepted "<<step;
        if(!retryPolicy) std::cout<<'/'<<schedule->steps;
        std::cout<<" time="<<time.end
            <<" divU="<<report.divergenceAfterProjection<<" dual mass closure="<<report.conservativeMassRelativeClosure
            <<" liquid="<<conservativeState->state().conservativePhase().report().liquidAfter<<'\n';
    }
    return 0;
    }
}

void HydroCase::Impl::finish()
{
    if (lifecycle != Lifecycle::Ready || acceptedStep!=completedSteps ||
        (!columnStopped && (schedule->controls.stepRetry?acceptedTime!=schedule->controls.endTime:completedSteps!=schedule->steps)))
        throw std::logic_error("Hydro case cannot finish before all coupled steps complete");
    if(schedule->controls.columnAmplitudeLimit) {
        lifecycle=Lifecycle::Failed; // IO failure retains accepted physics, never a completed run.
        std::ofstream termination(OutputLayout::series()/"termination.json");
        termination<<std::setprecision(17)<<"{\"reason\":\""<<(columnStopped?"columnAmplitude":"endTime")
            <<"\",\"accepted_step\":"<<acceptedStep<<",\"accepted_time\":"<<acceptedTime
            <<",\"requested_end_time\":"<<schedule->controls.endTime
            <<",\"relative_amplitude\":"<<acceptedColumnMode->relativeAmplitude
            <<",\"amplitude_limit\":"<<*schedule->controls.columnAmplitudeLimit<<"}\n";
        termination.close();if(!termination)throw std::runtime_error("cannot write Hydro stopping result");
        if(columnStopped)std::cout<<"[Hydrodynamics] Column amplitude reached its stopping limit at accepted step "
            <<acceptedStep<<", |A|/R="<<std::abs(acceptedColumnMode->relativeAmplitude)<<'\n';
    }
    std::cout << "[Hydrodynamics] 已完成 " << completedSteps << " 个流动时间步，最终时间="
              << acceptedTime << " (code units)。\n"
              << "[Hydrodynamics] 场数据：./output/fields/；逐步诊断：./output/series/"
              <<(momentCase?"hydro_moment_energy.csv":"hydro_diagnostics.csv")<<'\n';

    lifecycle = Lifecycle::Complete;
}

HydroCase::HydroCase(HydroCaseCoupling* coupling) : impl_(std::make_unique<Impl>(coupling)) {}
HydroCase::~HydroCase() = default;
int HydroCase::steps() const { return impl_->schedule->steps; }
double HydroCase::dt() const { return impl_->schedule->controls.dt; }
HydroCaseProgress HydroCase::progress() const noexcept {
    return {impl_->acceptedStep,impl_->completedSteps,impl_->acceptedTime,impl_->outputTime,
        impl_->meshOwner?impl_->meshOwner->cellIndexNum():0,impl_->meshChanges};
}
bool HydroCase::hasNextStep() const {return impl_->hasNextStep();}
int HydroCase::advanceNext() {return impl_->advanceNext();}
int HydroCase::initialize() { return impl_->initialize(); }
int HydroCase::advance(int step, double forceTime) { return impl_->advance(step, forceTime); }
void HydroCase::finish() { impl_->finish(); }
