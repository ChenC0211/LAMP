#include "Electrohydrodynamics/Coupling/DielectricCoupling.h"
#include "Electrohydrodynamics/IO/DielectricDiagnostics.h"
#include "Electrohydrodynamics/IO/SharpDielectricDiagnostics.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Electricity/IO/SurfaceChargeDiagnostics.h"

#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>

#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/DielectricAssemblyDiagnostics.h"
#include "Modules/Electricity/DielectricEvolution.h"
#include "Hydrodynamics/Coupling/HydroCaseCoupling.h"
#include "Hydrodynamics/HydroCase.h"
#include "Modules/Interface/Compatibility/InterfaceNodalState.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Core/Field/StaggeredVelocitySampler.h"
#include "Core/Mesh/QuadtreeMesh.h"

namespace
{
// Application adapter only: owns no mesh, geometry, Hydro state or clock.
// Electricity remains independent of the Hydro class and of spray species.
class DielectricCoupling final : public HydroCaseCoupling
{
    using ChargeMotionContext = Quantum::Electricity::DielectricChargeMotion;
public:
    explicit DielectricCoupling(DielectricConfig config) : electricConfig_(std::move(config)) {}

    bool needsTangentialLoadGeometry() const override { return charging(); }
    bool needsInterfaceMaterialMap() const override { return charging() && interfaceMoves_; }

    void completeStepWithVelocity(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
        const InterfaceNodalMaterialMap* materialMap, const Hydrodynamics& hydro,
        int step, double forceTime, double dt) override
    {
        if (charging() && interfaceMoves_)
        {
            requireInterface(geometry, interfaceIndex);
            // Consume the actual Interface step, not a second velocity
            // reconstruction or another material-trajectory setup.
            if (!electric_ || !electric_->charge()) throw std::runtime_error("moving charge state is not initialized");
            if (!materialMap) throw std::runtime_error("moving charge requires the Interface material map");
            const bool axisContact = !electric_->charge()->geometry->closed();
            const auto mapped = materialMap->motionFor(mesh, geometry, interfaceIndex,
                electric_->evaluated()->electric.fields->levelSet, step, forceTime, dt, axisContact);
            const ChargeMotionContext motion{mesh, mapped.snapshot, mapped.forwardMaterialPoint};
            completeChargeStep(step, forceTime, dt, nullptr, &motion);
            return;
        }
        if (!charging() || !electricConfig_.chargeConvection)
        { completeChargeStep(step, forceTime, dt, nullptr); return; }
        const auto velocity = hydro.averagedNativeVelocity(mesh);
        completeChargeStep(step, forceTime, dt, &velocity);
    }

    void publishCompletedStep(int step, double forceTime, double dt) override
    {
        if (!charging()) return;
        if (!pendingCharge_ || pendingCharge_->sourceTime.step != step ||
            pendingCharge_->sourceTime.forceTime != forceTime || pendingCharge_->sourceTime.dt != dt)
            throw std::logic_error("charge publication requires the matching completed physical substep");
        // Consume once before IO. A failed write cannot be retried or confused
        // with another physical charge update. The owning App also stops.
        auto completed = std::move(*pendingCharge_);
        pendingCharge_.reset();
        const auto& update = completed.transport;
        const auto& q = update.transport;
        const bool writeFields = schedule_->writeStepFields(step, true);
        chargeDiagnostics_.writeStep(*completed.charge, *completed.current, update, step, forceTime, dt, writeFields);
        if (step % schedule_->controls.logEverySteps == 0 || step == schedule_->steps)
            std::cout << "[EHD charge] 表面电荷推进 step=" << step << '/' << schedule_->steps
                      << " | chargeTime=" << forceTime + dt << " (code units)"
                      << " | 最大电荷密度绝对值=" << q.maxAbsDensity
                      << " | 总电荷=" << q.totalChargeAfter << " | 电极收支残差=" << completed.electrodeResidual
                      << " | 对流CFL=" << q.maxAdvectiveCfl
                      << " | 守恒重映射=" << update.remapped << " | 面元数=" << completed.charge->geometry->panels().size()
                      << " | 最大投影距离=" << update.maxProjectionDistance << '\n';
    }

private:
    void completeChargeStep(int step, double forceTime, double dt, const StaggeredVelocitySampler* velocity,
        const ChargeMotionContext* motion = nullptr)
    {
        if (!charging()) return;
        if (!electric_) throw std::runtime_error("electric coupling is not initialized");
        if (pendingCharge_) throw std::logic_error("previous charge result has not been published");
        pendingCharge_ = electric_->completeCharge({step, forceTime, dt}, velocity, motion);
    }

public:
    void preflight(const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface, const DictionaryGeometry&,
        const DictionaryMesh&, int interfaceIndex, std::shared_ptr<const HydroCaseSchedule> schedule) override
    {
        if (!schedule) throw std::invalid_argument("EHD requires the owning Hydro case schedule");
        electricConfig_.axisymmetric = config.coordinateSystem == "axisymmetricRZ";
        // The Hydro RZ App domain starts at r=0; that axis cannot be a
        // Dirichlet electrode. Reject in --check-case, not at the first solve.
        if (electricConfig_.axisymmetric && electricConfig_.electrodeAxis != "y")
            throw std::runtime_error("RZ EHD requires electrodeAxis y; r=0 is a homogeneous-Neumann symmetry axis");
        if (electricConfig_.axisymmetric && (electricConfig_.interfaceDiscretization != "sharpHermiteP3"
            || config.interfaceForceLayer != "regularized"))
            throw std::runtime_error("RZ EHD requires sharpHermiteP3 and regularized loads");
        if (config.materialModel != "twoPhase" || interfaceIndex < 0)
            throw std::runtime_error("Electrohydrodynamics requires twoPhase and a named fluid interface");
        if (electricConfig_.interfaceDiscretization == "sharpHermiteP3" && config.interfaceForceGeometry != "hermite")
            throw std::runtime_error("sharpHermiteP3 requires Hydro interfaceForceGeometry hermite so electric and capillary loads share geometry");
        if (config.interfaceForceLayer == "regularized" && electricConfig_.interfaceDiscretization != "sharpHermiteP3")
            throw std::runtime_error("regularized EHD load layer requires sharpHermiteP3 so electric and capillary loads share G(H)");
        if (charging())
        {
            if (interface.reinitializeEvery != 0)
                throw std::runtime_error("leakyDielectric charge transport currently requires reinitializeEvery 0");
            if (interface.advect && !electricConfig_.chargeConvection)
                throw std::runtime_error("moving leakyDielectric requires chargeConvection true for material charge transport");
            const double relaxationCfl = Quantum::Electricity::DielectricEvolution::explicitRelaxationCfl(electricConfig_, schedule->controls.dt);
            std::cout << "[EHD] leakyDielectric: finite explicit charging, "
                      << (interface.advect ? "moving interface (conservative material remap), " : "fixed interface, ")
                      << "chargeConvection=" << (electricConfig_.chargeConvection ? "true" : "false")
                      << "; no bulk charge, emission or AMR.\n"
                      << "[EHD] Charge relaxation CFL=" << relaxationCfl << " (step restriction, not a convergence certificate).\n";
        }
        interfaceThickness_ = config.interfaceThickness;
        interfaceMoves_ = interface.advect;
        schedule_ = std::move(schedule);
        if (electricConfig_.axisymmetric)
            std::cout << "[EHD] RZ electric operator: revolved FV measures, axisymmetric P3 jets, even scalar at r=0.\n";
        if (!charging()) std::cout << "[EHD] perfectDielectric: no free charge, conduction, evaporation, particles, AMR or restart.\n";
        std::cout << (electricConfig_.interfaceDiscretization == "sharpHermiteP3"
                      ? "[EHD] 物理步骤：共享Hermite界面 -> sharp FV/P3电载荷与毛细力 -> 流动求解 -> 按配置推进界面/电荷。\n"
                      : "[EHD] 物理步骤：共享水平集 -> 弥散介电系数/FV电场 -> 电体力与毛细力 -> 流动求解 -> 按配置推进界面。\n")
                  << "[EHD] 电载荷帧使用受力时刻 force_time=t^n；流场帧使用完成时刻 t^(n+1)，均为代码时间。\n";
        std::cout << "[EHD] 动量方程使用压力相容的原生面电力；单元电力输出仅供诊断。\n";
    }

    void initialize(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
        int interfaceIndex, const HydrodynamicsConfig&) override
    {
        requireInterface(geometry, interfaceIndex);
        electric_ = std::make_unique<Quantum::Electricity::DielectricEvolution>(mesh, electricConfig_, interfaceThickness_);
        diagnostics_.initialize();
        if (electricConfig_.interfaceDiscretization == "sharpHermiteP3") sharpDiagnostics_.initialize();
        if (charging()) chargeDiagnostics_.initialize();
    }

    void computeForce(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
        int interfaceIndex, int step, double forceTime, double dt,
        CellField<double>& fx, CellField<double>& fy) override
    {
        evaluateForce(mesh, geometry, interfaceIndex, step, forceTime, dt, nullptr, fx, fy);
    }

    void computeForceWithGeometry(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
        int interfaceIndex, int step, double forceTime, double dt, const HydroInterfaceLoadGeometry& loadGeometry,
        CellField<double>& fx, CellField<double>& fy) override
    {
        evaluateForce(mesh, geometry, interfaceIndex, step, forceTime, dt, &loadGeometry, fx, fy);
    }

    bool computeFaceForce(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>*,
        int, int step, double forceTime, double dt, VEdgeField<double>& fx, HEdgeField<double>& fy) override
    {
        if (!load_ || completedStep_ < 0 || completedStep_ != step || completedForceTime_ != forceTime || completedDt_ != dt)
            throw std::runtime_error("electric force does not match the completed time layer");
        for (int i = 0; i < mesh.vEdgeIndexNum(); ++i) fx.value[i] = load_->faceForceX.value[i];
        for (int i = 0; i < mesh.hEdgeIndexNum(); ++i) fy.value[i] = load_->faceForceY.value[i];
        completedGeometry_ = nullptr;
        completedStep_ = -1;
        return true;
    }

    bool computeFaceForceWithGeometry(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
        int interfaceIndex, int step, double forceTime, double dt, const HydroInterfaceLoadGeometry& loadGeometry,
        VEdgeField<double>& fx, HEdgeField<double>& fy) override
    {
        requireInterface(geometry, interfaceIndex);
        loadGeometry.validateSource(mesh, geometry->levelSets[interfaceIndex]);
        if (&loadGeometry != completedGeometry_ || step != completedStep_ ||
            forceTime != completedForceTime_ || dt != completedDt_)
            throw std::runtime_error("electric face force does not match the completed shared-geometry step");
        return computeFaceForce(mesh, geometry, interfaceIndex, step, forceTime, dt, fx, fy);
    }

private:
    DielectricConfig electricConfig_;
    double interfaceThickness_ = 0.0;
    bool interfaceMoves_ = false;
    std::shared_ptr<const HydroCaseSchedule> schedule_;
    std::unique_ptr<Quantum::Electricity::DielectricEvolution> electric_;
    std::unique_ptr<Quantum::Hydrodynamics::HydroUniformLoad> load_;
    DielectricDiagnostics diagnostics_;
    SharpDielectricDiagnostics sharpDiagnostics_;
    SurfaceChargeDiagnostics chargeDiagnostics_;
    std::optional<Quantum::Electricity::DielectricChargeUpdate> pendingCharge_;
    int completedStep_ = -1;
    double completedForceTime_ = 0.0, completedDt_ = 0.0;
    // Borrowed only between the two adjacent callbacks, cleared on consumption.
    const HydroInterfaceLoadGeometry* completedGeometry_ = nullptr;

    bool charging() const { return electricConfig_.model == "leakyDielectric"; }

    void evaluateForce(const QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
        int interfaceIndex, int step, double forceTime, double dt, const HydroInterfaceLoadGeometry* loadGeometry,
        CellField<double>& fx, CellField<double>& fy)
    {
        if (pendingCharge_) throw std::logic_error("previous charge result has not been published");
        completedStep_ = -1;
        completedGeometry_ = nullptr;
        requireInterface(geometry, interfaceIndex);
        if (!electric_) throw std::runtime_error("electric coupling is not initialized");
        const bool writeFields = schedule_->writeStepFields(step, true);
        DielectricSharpAssemblyDiagnostics assemblyDiagnostics;
        auto* assemblyOutput = writeFields && electricConfig_.interfaceDiscretization == "sharpHermiteP3"
            ? &assemblyDiagnostics : nullptr;
        const auto evaluation = electric_->evaluate(mesh, geometry->levelSets[interfaceIndex],
            {step, forceTime, dt}, loadGeometry ? &loadGeometry->geometry() : nullptr, assemblyOutput,
            loadGeometry && electricConfig_.interfaceDiscretization == "sharpHermiteP3" ? loadGeometry->traceGeometry() : nullptr);
        const auto& report = evaluation->electric.report;
        const auto& fields = *evaluation->electric.fields;
        if (report.usedSharpInterface)
        {
            if (!loadGeometry) throw std::runtime_error("missing sharp load geometry");
            load_ = std::make_unique<Quantum::Hydrodynamics::HydroUniformLoad>(
                Quantum::Hydrodynamics::buildUniformInterfaceTractionLoad(mesh, geometry->levelSets[interfaceIndex],
                    *loadGeometry, fields.traceGeometry, fields.interfaceTractions(), report.usedSurfaceCharge));
        }
        else
        {
            if (!fields.polarization) throw std::runtime_error("missing solved polarization coefficients");
            const auto& polarization = *fields.polarization;
            load_ = std::make_unique<Quantum::Hydrodynamics::HydroUniformLoad>(
                Quantum::Hydrodynamics::buildUniformScalarGradientLoad(mesh, polarization.meshIdentity,
                    polarization.permittivityByCell, polarization.coefficientByNative));
        }
        if (charging() && step == 1)
            std::cout << "[EHD] Charge integration panels=" << electric_->charge()->geometry->panels().size() << '\n';
        for (int i = 0; i < mesh.cellIndexNum(); ++i)
        {
            // Momentum consumes the face loads below, exactly once.
            fx.value[i] = 0.0;
            fy.value[i] = 0.0;
        }
        diagnostics_.writeStep(report, *load_, step, forceTime, dt);
        if (report.usedSharpInterface)
        {
            if (!loadGeometry) throw std::runtime_error("missing geometry for completed sharp electric solve");
            sharpDiagnostics_.writeStep(mesh, geometry->levelSets[interfaceIndex], *loadGeometry, fields, *load_,
                electricConfig_, report, step, forceTime, dt, interfaceThickness_,
                schedule_->controls.outputPrefix, writeFields, assemblyOutput);
        }
        if (writeFields)
            writeDielectricLoadFields(mesh, step, forceTime, fields, *load_, electricConfig_.axisymmetric, schedule_->controls.outputPrefix);
        if (step == 1)
            std::cout << "[EHD] 首步实际电场求解器 electricBackend=" << (report.usedSuperLU ? "SuperLU" : "unverified")
                      << "；逐步电学诊断：./output/series/electric_diagnostics.csv\n";
        if (step % schedule_->controls.logEverySteps == 0 || step == schedule_->steps)
            std::cout << "[EHD] 电场/界面载荷 step=" << step << '/' << schedule_->steps
                      << " | forceTime=" << forceTime << " (code units)"
                      << " | 电势范围=[" << report.minPotential << ',' << report.maxPotential << ']'
                      << " | maxE=" << report.maxField
                      << " | 电学相对残差=" << report.relativeResidual
                      << " | 边界电通量不平衡=" << report.boundaryFluxImbalance
                      << " | 最大电力密度=" << load_->maxForce << '\n';
        completedStep_ = step;
        completedForceTime_ = forceTime;
        completedDt_ = dt;
        completedGeometry_ = loadGeometry;
    }

    static void requireInterface(const ImmersedGeometry<QuadtreeMesh>* geometry, int index)
    {
        if (!geometry || index < 0 || index >= static_cast<int>(geometry->levelSets.size()))
            throw std::runtime_error("missing shared electric/fluid interface");
    }

};

}

std::unique_ptr<HydroCaseCoupling> makeDielectricCoupling(DielectricConfig config)
{
    return std::make_unique<DielectricCoupling>(std::move(config));
}
