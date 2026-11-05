#pragma once

#include "Modules/Hydrodynamics/HydroPressureLoad.h"
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
class HydroMomentumMass;
class HydroMomentumGeometry;
class HydroEnergyViscosity;
// None records an unrequested SPD certificate; it is not a proof of positivity.
enum class HydroVelocitySpdTier : unsigned char
{ None=0, Diagonal=1, StrictDiagonalDominance=2, EnergyGram=3, LLT=4 };
struct HydroVelocitySpdReport
{
    HydroVelocitySpdTier tier=HydroVelocitySpdTier::None;
    bool gramAttempted=false,gramBoundsAvailable=false;
    double gramBoundUpper=0,scaledMassLower=0,certificateSeconds=0,lltSeconds=0;
};
// Actual validation-path diagnostics only. The plan stores indices, never values.
// eligible requires strict source/current FP and actual normal/zero S/W.
// planHit means the S pattern already existed before this factory lookup.
// retainedBytes is the participating single plan accounted size, not work.
struct HydroSparseValidationReport
{
    bool requested=false,eligible=false,planHit=false;
    bool exactSymmetryUsed=false,exactSkewUsed=false;
    std::size_t nonZeros=0,retainedBytes=0;
    // Fresh rounded S+=W and final transpose-defect values used existing routing.
    // Scratch is one private scalar buffer for this call, never retained stock.
    bool appendPlanUsed=false;
    std::size_t appendScratchBytes=0;
};
// Per-factory input-A split diagnostics. No physical values are cached.
// planHit describes an existing exact input pattern, not value/FP eligibility.
// Used flags mean fresh S then fresh W values were actually evaluated by it.
struct HydroVelocitySplitPlanReport
{
    bool requested=false,eligible=false,planHit=false;
    bool symmetricValuesUsed=false,skewValuesUsed=false;
    std::size_t inputNonZeros=0,unionNonZeros=0,retainedBytes=0;
    // Actual current Frozen numerical stages; serial execution records one.
    bool parallelUsed=false,parallelFallback=false;
    int maximumWorkers=0;
};
// Per-preparation pure index-plan diagnostics; retained storage is not work.
struct HydroSaddleAssemblyPlanReport
{
    bool requested=false,eligible=false,planHit=false;
    bool systemScatterUsed=false,quotientScatterUsed=false;
    std::size_t systemNonZeros=0,quotientNonZeros=0,retainedBytes=0;
    bool symbolicBuilt=false;
    int symbolicWorkers=0; // Minimum actual team of cold count and fill, zero on reuse.
    bool scatterParallelUsed=false,scatterParallelFallback=false;
    int scatterMaximumWorkers=0;
    bool quotientAliasesSystem=false; // Actual identity quotient; no quotient scatter.
};
// Caller-owned physics prepares A_u; this contract copies and validates it.
class HydroSaddleVelocityBlock
{
public:
    using Matrix = HydroPressureLoad::Matrix;
    static HydroSaddleVelocityBlock build(HydroDofLayout layout, const Matrix& integralVelocity,int tangentialCount=0);
    // Same symmetric/skew validation and S+=W assembly, with one transpose
    // and one mutable data owner for the actual Frozen transport split.
    static HydroSaddleVelocityBlock buildWithTransportSplit(HydroDofLayout layout,
        const Matrix& integralVelocity,int tangentialCount=0);
    // Actual mass/energy sufficient certificate for the same Lower SPD part;
    // failure/unsupported arithmetic retains the original complete LLT gate.
    static HydroSaddleVelocityBlock buildWithEnergyTransportSplit(const QuadtreeMesh&,
        HydroDofLayout,const Matrix& actualVelocity,const HydroMomentumMass&,
        const HydroEnergyViscosity&,double dt);
    // Explicit actual geometry owns one bounded pure transpose-index plan.
    // Foreign owner is rejected; unsupported/nonexact inputs use original gates.
    static HydroSaddleVelocityBlock buildWithEnergyTransportSplit(const QuadtreeMesh&,
        HydroDofLayout,const Matrix& actualVelocity,const HydroMomentumMass&,
        const HydroEnergyViscosity&,double dt,const HydroMomentumGeometry& sparsePlanOwner);
    const HydroVelocitySpdReport& spdReport() const noexcept;
    const HydroSparseValidationReport& sparseValidationReport() const noexcept;
    const HydroVelocitySplitPlanReport& splitPlanReport() const noexcept;
    // Explicit coercive extension, retaining the certified SPD part. Supplied
    // transport must satisfy the original 1e-12 skew tolerance; no projection.
    HydroSaddleVelocityBlock withSkew(const Matrix& integralSkew) const;
    const HydroDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept;
    // Unrequested on the actual Frozen production block: quiet NaN. Public
    // arbitrary-matrix factories and withSkew retain the complete diagnostic.
    double symmetryDefect() const noexcept;
    // Diagnostic wall time only; never used by a physical acceptance gate.
    double validationSeconds() const noexcept;
    int tangentialCount() const noexcept;
private:
    friend class HydroFrozenMidpoint;
    // Actual Frozen momentum A has already been built from validated typed
    // physical inputs in this same callback-free invocation. Move it once;
    // check shape/current finiteness, without an unrequested S/W or SPD audit.
    // true is reserved for its exact local final-A writer, which has already
    // observed every current stored value before this move; no global token.
    static HydroSaddleVelocityBlock buildFrozenActual(HydroDofLayout,Matrix&&,bool currentValuesFinite=false);
    static HydroSaddleVelocityBlock buildEnergyTransportSplitImpl(const QuadtreeMesh&,
        HydroDofLayout,const Matrix&,const HydroMomentumMass&,const HydroEnergyViscosity&,
        double,const HydroMomentumGeometry*);
    struct Data;
    explicit HydroSaddleVelocityBlock(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};

class HydroSaddleRhs
{
public:
    HydroSaddleRhs(HydroDofLayout layout, const std::vector<double>& integralMomentum,
        const std::vector<double>& prescribedNormalVelocity,const std::vector<double>& tangentialMomentum={});
    const std::vector<double>& integralMomentum() const noexcept { return momentum_; }
    const std::vector<double>& prescribedNormalVelocity() const noexcept { return prescribed_; }
    const std::vector<double>& tangentialMomentum() const noexcept {return tangential_;}
    void validateFor(const HydroDofLayout& expected, const QuadtreeMesh& mesh) const;
private:
    HydroDofLayout layout_;
    std::vector<double> momentum_, prescribed_,tangential_;
};

// Passive timings supplied by an actual preparation consumer. Old callers
// leave available=false; no numerical acceptance reads these values.
struct HydroSaddlePreparationTimings
{
    bool available=false;
    double inputValidationSeconds=0,forceValidationSeconds=0;
    double matrixActionSeconds=0,rhsSeconds=0,velocityPrepareSeconds=0;
};
struct HydroSaddleConfig
{
    double algebraicRelativeTolerance=1e-11, momentumRelativeTolerance=1e-11;
    double continuityRelativeTolerance=1e-11, pressureMeanRelativeTolerance=1e-11;
    double boundaryFluxRelativeTolerance=1e-12;
    // Explicit caller selection; unrelated consumers retain their existing path.
    // GMRES failure rejects the candidate instead of invoking a global direct solve.
    bool gmresOnlyParallel=false;
    // Explicit pressure preconditioner only. Requires GMRES-only so no
    // full quotient direct fallback is enabled by this selection.
    bool gmresGlobalPressure=false;
};
struct HydroSaddleReport
{
    double algebraicRelativeResidual=0, freeMomentumDefect=0, maxCellContinuityDefect=0;
    double maxPressureMeanDefect=0, maxBoundaryFluxDefect=0, maxAbsGaugeMultiplier=0;
    int freeVelocityCount=0, pressureGaugeCount=0,freeTangentialCount=0;
    int linearRefinements=0;
    // Internal quotient factor input, while matrix() remains augmented K.
    int factorDimension=0;
    Eigen::Index factorNonZeros=0;
    // Same original equations/gates; backend timings do not control acceptance.
    bool usedSuperLU=false;
    // Krylov only proposes a solution; the original local gates still decide.
    bool usedBlockGmres=false,usedDirectFallback=false;
    int krylovIterations=0,krylovRestarts=0,krylovSolves=0;
    // 0 none; 1 preconditioner qualification; 2 inner; 3 original outer gates.
    int directFallbackReason=0;
    double krylovRelativeResidual=0,preconditionerSeconds=0,directFactorSeconds=0;
    // Private algebraic norm only; physical acceptance still uses original gates.
    double krylovScaledRelativeResidual=0;
    bool krylovStagnated=false;
    int parallelPressureBlockCount=0,parallelPressureMaximumWorkers=0;
    // Explicit GMRES-only global-S pressure action, never full-K direct solve.
    // Factor nnz counts stored strict L (unit diagonal is implicit); D is separate.
    // Storage is actual indexed factor payload, not allocator/RSS peak memory.
    bool pressureGlobalRequested=false,pressureGlobalAvailable=false,pressureGlobalPrepared=false;
    bool pressureGlobalTimingsAvailable=false;
    std::size_t pressureGlobalRows=0,pressureGlobalNonZeros=0,pressureGlobalFactorNonZeros=0;
    std::size_t pressureGlobalStorageBytes=0,pressureGlobalAdmissionBytes=0,pressureGlobalByteLimit=0;
    std::size_t pressureGlobalSolveCalls=0;
    double pressureGlobalSymbolicSeconds=0,pressureGlobalNumericSeconds=0,pressureGlobalSolveSeconds=0;
    double velocityBlockSeconds=0, factorSeconds=0, solveSeconds=0;
    // Passive THIS-solve breakdown, including all attempted inner solves and
    // outer assessments. Global-pressure solve is nested in precondition apply.
    // Orth calls count complete successful two-pass CGS2 calls with actual team>1.
    bool solveBreakdownTimingsAvailable=false;
    double krylovPreconditionApplySeconds=0,krylovSpmvSeconds=0,krylovOrthogonalizationSeconds=0;
    double krylovTrueResidualSeconds=0,physicalAssessmentSeconds=0;
    int orthogonalizationParallelCalls=0,orthogonalizationMaximumWorkers=0;
    HydroSparseValidationReport velocitySparseValidation;
    HydroVelocitySplitPlanReport velocitySplitPlan;
    HydroSaddleAssemblyPlanReport sparseAssemblyPlan;
    int velocitySpdTier=0;
    bool velocitySpdGramAttempted=false,velocitySpdGramBoundsAvailable=false;
    double velocitySpdGramBoundUpper=0,velocitySpdScaledMassLower=0;
    double velocitySpdCertificateSeconds=0,velocitySpdLltSeconds=0;
    bool frozenPreparationTimingsAvailable=false,saddleAssemblyTimingsAvailable=false;
    double frozenInputValidationSeconds=0,frozenForceValidationSeconds=0;
    double frozenMatrixActionSeconds=0,frozenRhsSeconds=0,frozenVelocityPrepareSeconds=0;
    double saddlePartitionSeconds=0,saddleSystemAssemblySeconds=0;
    double saddleQuotientAssemblySeconds=0,saddleRhsEliminationSeconds=0;
    // Current-factor optional ordered-row view; no cross-factor coefficient cache.
    // requested: eligible GMRES factor attempted view preparation; eligible:
    // the complete current view exists. Used/calls are THIS public solve only,
    // including its original outer refinements. Fallback denotes original Eigen
    // main SpMV execution, never a physical/Direct-solver acceptance fallback.
    bool orderedSpmvRequested=false,orderedSpmvEligible=false;
    bool orderedSpmvUsed=false,orderedSpmvFallback=false;
    std::size_t orderedSpmvCalls=0,orderedSpmvParallelCalls=0;
    // Current view stock, not cumulative; workers are successful actual-team max.
    std::size_t orderedSpmvViewBytes=0;
    int orderedSpmvMaximumWorkers=0;
    // Current-factor input observation, not per-solve Direct-use counters.
    // Deferred means the actual complete pre-certificate omitted eager storage;
    // materialized also includes unsupported/canonical-cold eager inputs.
    bool directInputDeferred=false,directInputMaterialized=false;
    // Current pressure-factor preparation only; no output/dictionary schema.
    // Reason codes are Detail/HydroPressureCoarseCorrection::Reason (0..9).
    // Applications count THIS public solve, including original refinements.
    bool pressureCoarsePrepared=false;
    int pressureCoarseDimension=0,pressureCoarseReason=0;
    std::size_t pressureCoarseApplications=0;
    double pressureCoarseConditionUpper=0,pressureCoarseFormationNormUpper=0;
};
class HydroSaddleResult
{
public:
    const HydroNativeValues& values() const noexcept { return values_; }
    const HydroSaddleReport& report() const noexcept { return report_; }
    const std::vector<double>& tangentialVelocity() const noexcept {return tangential_;}
private:
    friend class HydroSinglePhaseSaddle;
    HydroSaddleResult(HydroNativeValues values, HydroSaddleReport report,std::vector<double> tangential={});
    HydroNativeValues values_;
    HydroSaddleReport report_;
    std::vector<double> tangential_;
};

// Total normal traction is -pressure*n_out, not a pressure cell constraint.
// The native normal velocity on this physical outlet remains free.
struct HydroSaddlePressureOutlet
{
    int velocityDof=-1;
    double pressure=0;
};

// Monolithic reference solve. Closed pressure components retain the volume
// gauge; a physical traction outlet anchors its component through original C.
class HydroSinglePhaseSaddle
{
public:
    using Matrix = Eigen::SparseMatrix<double, Eigen::ColMajor, Eigen::Index>;
    static HydroSinglePhaseSaddle buildClosed(HydroPressureLoad pressure,
        HydroSaddleVelocityBlock velocity, HydroSaddleConfig config={},
        const HydroSaddlePreparationTimings& preparation={});
    static HydroSinglePhaseSaddle buildClosed(HydroPressureLoad pressure,
        HydroSaddleVelocityBlock velocity,HydroSaddleConfig,
        const HydroSaddlePreparationTimings&,const HydroMomentumGeometry& sparsePlanOwner);
    static HydroSinglePhaseSaddle buildWithPressureOutlets(HydroPressureLoad pressure,
        HydroSaddleVelocityBlock velocity, const std::vector<HydroSaddlePressureOutlet>&,
        HydroSaddleConfig config={},const HydroSaddlePreparationTimings& preparation={});
    static HydroSinglePhaseSaddle buildWithPressureOutlets(HydroPressureLoad pressure,
        HydroSaddleVelocityBlock velocity,const std::vector<HydroSaddlePressureOutlet>&,
        HydroSaddleConfig,const HydroSaddlePreparationTimings&,const HydroMomentumGeometry& sparsePlanOwner);
    const HydroPressureLoad& pressureLoad() const noexcept;
    const HydroSaddleVelocityBlock& velocityBlock() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    const Matrix& matrix() const noexcept;
    const std::vector<int>& freeVelocityDofs() const noexcept;
    const std::vector<int>& pressureComponents() const noexcept;
    // Global component to volume-gauge index; an outlet component has -1.
    const std::vector<int>& pressureGaugeByComponent() const noexcept;
    const std::vector<HydroSaddlePressureOutlet>& pressureOutlets() const noexcept;
    // Nu integral traction loads T=C_boundary*p_out; exact zero on closed builds.
    const std::vector<double>& pressureOutletLoads() const noexcept;
    int pressureGaugeCount() const noexcept;
    // Adds retained outlet T once; caller momentum excludes this boundary load.
    HydroSaddleResult solve(const QuadtreeMesh& mesh, const HydroSaddleRhs& rhs) const;
    // Native GMRES-only initial iterate; current RHS/boundaries and all gates remain authoritative.
    HydroSaddleResult solve(const QuadtreeMesh& mesh,const HydroSaddleRhs& rhs,
        const HydroNativeValues& initialGuess) const;
private:
    struct Data;
    HydroSaddleResult solveImpl(const QuadtreeMesh&,const HydroSaddleRhs&,const HydroNativeValues*) const;
    explicit HydroSinglePhaseSaddle(std::shared_ptr<const Data> data);
    static HydroSinglePhaseSaddle buildImpl(HydroPressureLoad, HydroSaddleVelocityBlock,
        const std::vector<HydroSaddlePressureOutlet>*, HydroSaddleConfig,
        const HydroSaddlePreparationTimings&,const HydroMomentumGeometry* sparsePlanOwner=nullptr);
    std::shared_ptr<const Data> data_;
};
}
