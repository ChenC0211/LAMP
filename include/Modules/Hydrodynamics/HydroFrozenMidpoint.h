#pragma once
#include "Modules/Hydrodynamics/HydroKepConvection.h"
#include "Modules/Hydrodynamics/HydroPhaseViscosity.h"
#include "Modules/Hydrodynamics/HydroSinglePhaseSaddle.h"
#include "Modules/Hydrodynamics/HydroBalancedCsfLoad.h"
#include "Modules/Hydrodynamics/HydroInterfaceNormalLoad.h"
#include "Modules/Hydrodynamics/HydroMaxwellLoad.h"
#include <optional>
#include <cstddef>
#include <variant>

namespace LAMP::Hydrodynamics
{
enum class HydroFrozenMidpointScheme { LaggedConvection, FrozenConvection };
struct HydroFrozenMidpointForces
{
    Vector2D bodyAcceleration{0,0};
    // ACTUAL physical interface-force owners. Sharp balance sign is converted
    // once here; no caller-provided opaque force matrix/array enters this kernel.
    std::variant<std::monostate,HydroBalancedCsfLoad,HydroInterfaceNormalLoad> interfaceLoad;
    // A distinct actual coeval load: capillarity and Maxwell can coexist.
    std::optional<HydroMaxwellLoad> maxwellLoad;
};
struct HydroFrozenMidpointBoundary
{
    HydroAdvectorBoundary walls;
    std::vector<double> normalVelocity;
    std::vector<HydroSaddlePressureOutlet> pressureOutlets;
    static HydroFrozenMidpointBoundary closed(HydroAdvectorBoundary);
    static HydroFrozenMidpointBoundary prescribed(HydroAdvectorBoundary,std::vector<double>,
        std::vector<HydroSaddlePressureOutlet> = {});
};
// One actual Frozen preparation. Hits do not imply values were used.
// Stock is the current bounded pure-index slot, not cumulative work.
// Union hit requires the complete current H route, not only fixed K/P.
// A newly compiled current route may be used with unionPlanHit=false.
struct HydroFrozenMatrixPlanReport
{
    bool requested=false,eligible=false,physicalPlanHit=false,physicalValuesUsed=false;
    bool unionPlanHit=false,unionValuesUsed=false;
    // A newly compiled actual K/diagonal or H-union route may serve this call.
    bool physicalColdCompiled=false,physicalColdUsed=false,unionColdCompiled=false,unionColdUsed=false;
    // Successful actual current-value teams only; fallback is independent.
    bool valuesParallelUsed=false,valuesParallelFallback=false;
    int valuesMaximumWorkers=0;
    // Added optional worker stock peak; excludes original matrices and RSS.
    std::size_t valuesScratchBytes=0;
    std::size_t viscosityNonZeros=0,convectionNonZeros=0,physicalNonZeros=0,velocityNonZeros=0,retainedBytes=0;
};
// Shared native NS block used by actual FirstJet and VOF AC preparation.
// Retains genuine M/F/R/KEP/K/C/load/boundary owners;
// only one original saddle solver and one matrix/RHS arithmetic are used.
class HydroFrozenMidpoint
{
public:
    using Matrix=HydroSaddleVelocityBlock::Matrix;
    // GMRES-only native solve uses initialGuess, or begin when absent. The
    // hint is not retained; the current equations and acceptance gates decide.
    static HydroFrozenMidpoint prepare(const QuadtreeMesh&,HydroDualMaterial,HydroKepConvection,
        HydroKepConvectionResult appliedConvection,HydroPhaseViscosity,HydroPressureLoad,HydroNativeValues begin,
        HydroFrozenMidpointForces,HydroFrozenMidpointBoundary,HydroFrozenMidpointScheme,
        HydroMaterialSources,LAMP::TimeLayer midpointClock,HydroSaddleConfig saddleConfig={},
        const HydroNativeValues* initialGuess=nullptr);
    std::size_t solveCount() const noexcept;
    const HydroDualMaterial& material() const noexcept;
    const HydroKepConvection& convection() const noexcept;
    const HydroKepConvectionResult& appliedConvection() const noexcept;
    const HydroPhaseViscosity& viscosity() const noexcept;
    const HydroPressureLoad& pressureLoad() const noexcept;
    const HydroNativeValues& beginValues() const noexcept;
    const HydroFrozenMidpointForces& forces() const noexcept;
    const HydroFrozenMidpointBoundary& boundary() const noexcept;
    HydroFrozenMidpointScheme scheme() const noexcept;
    const HydroFrozenMatrixPlanReport& matrixAssemblyPlanReport() const noexcept;
    const Matrix& viscosityMatrix() const noexcept;
    const Matrix& physicalVelocityMatrix() const noexcept;
    const Matrix& velocityMatrix() const noexcept;
    const Matrix& frozenConvectionMatrix() const noexcept;
    const Eigen::VectorXd& beginViscousAction() const noexcept;
    const Eigen::VectorXd& beginFrozenAction() const noexcept;
    const std::vector<long double>& bodyLoads() const noexcept;
    const std::vector<double>& interfaceLoads() const noexcept;
    bool hasInterfaceLoad() const noexcept;
    const std::vector<double>& maxwellLoads() const noexcept;
    bool hasMaxwellLoad() const noexcept;
    const std::vector<double>& integralMomentumRhs() const noexcept;
    const HydroSinglePhaseSaddle& saddle() const noexcept;
    const HydroSaddleResult& solveResult() const noexcept;
    void validateFor(const QuadtreeMesh&,HydroMaterialSources,LAMP::TimeLayer) const;
private:
    struct Data;
    static std::shared_ptr<Data> assemble(const QuadtreeMesh&,HydroDualMaterial,HydroKepConvection,
        HydroKepConvectionResult,HydroPhaseViscosity,HydroPressureLoad,HydroNativeValues,
        HydroFrozenMidpointForces,HydroFrozenMidpointBoundary,HydroFrozenMidpointScheme,
        HydroMaterialSources,LAMP::TimeLayer,HydroSaddleConfig);
    explicit HydroFrozenMidpoint(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}


