#pragma once

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <optional>
#include <functional>
#include "Core/Field/StaggeredVelocitySampler.h"
#include "Modules/Hydrodynamics/HydroSolveReport.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "Modules/Hydrodynamics/Detail/HydroMomentumAssembly.h"
#include "Modules/Hydrodynamics/HydroPressureSystem.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroStateInput.h"
#include "Modules/Hydrodynamics/HydroMaterialFields.h"
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/HydroNativePressureStencil.h"
#include "Core/Base/Base.h"
#include "Core/Geometry/LevelSet.h"
#include "Core/Geometry/LevelSetHyperTan.h"
#include "Core/Field/VertField.h"
#include "Core/Field/CellField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Modules/Hydrodynamics/LegacyNodalState.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

class QuadtreeMesh;
class HydroFittedDomain;
class HydroFittedMomentumLayout;
class HydroFittedMeshTransfer;
class NativeMeshIdentity;
class InterfaceState;
namespace LevelSetHermiteGeometry { class Snapshot; }

class Hydrodynamics
{
private:
    // Physical fields, histories and work buffers belong to Hydro. Consumers
    // use state() for synchronous reads or an owning velocity reconstruction.
    VertField<double> vertRho;
    VertField<double> vertMu;
    // Physical accumulated pressure on the standalone closed-domain path.
    // The legacy geometry-aware path retains its nonincremental pressure.
    CellField<double> p;
    CellField<double> rho;
    CellField<double> mu;
    CellField<double> externalFx;
    CellField<double> externalFy;
    CellField<double> directFx;
    CellField<double> directFy;
    CellField<double> evapUx;
    CellField<double> evapUy;
    VEdgeField<double> ux;
    HEdgeField<double> uy;
    VEdgeField<double> uxPre;
    HEdgeField<double> uyPre;
    VEdgeField<double> uxMid;
    HEdgeField<double> uyMid;

    CellField<double> uMag;
    CellField<double> divU;

    CellField<double> dxUx;
    CellField<double> dyUx;
    CellField<double> dxUy;
    CellField<double> dyUy;
    CellField<double> magGradU;

    LevelSetHyperTan hyperTanLS;
    
public:
    template<typename Mesh>
    Hydrodynamics(Mesh& mesh);
    ~Hydrodynamics();

    // Production single-step service. Owns sequence/material/history checks;
    // delegates unchanged complete/fitted assembly to the existing kernels.
    // Complete XY/RZ: closed native 2:1 leaves or uniform boundary patches.
    // Fitted RZ: explicit fluid domain with wall/axis/reservoir boundaries.
    // Ordinary steps retain CN/AB2; conservative transport selects the existing
    // one-step mass/upwind method. Both keep incremental physical pressure
    // and the projection/coupled correction policy of the selected domain.
    // Requires explicit initializeState and current updateMaterials first.
    // No IO, App config, main clock, or automatic retry/rollback.
    HydroSolveReport advance(QuadtreeMesh& mesh, const HydroStepInput& input);
    // Same local solve with explicit caller acceptance. While pending, queries
    // and next-time material preparation are allowed; another step and mesh
    // capture are not. Abandoning it stops this instance unless an enclosing
    // beginStepTransaction guard subsequently restores the old state.
    HydroStepCandidate prepareStep(QuadtreeMesh& mesh, const HydroStepInput& input);
    // Opt-in owning rollback for coupled preparation; legacy calls retain
    // their fail-stop semantics. No mesh mutation while the guard is alive.
    HydroStepTransaction beginStepTransaction(const QuadtreeMesh& mesh);
    HydroStateView state() const;
    // Complete-domain reconstructions of the actual named native time layer.
    // Results own their samples; no writable Hydro storage escapes. They use
    // the existing wall/axis closures and do not change or advance the state.
    std::pair<std::unique_ptr<VertField<double>>, std::unique_ptr<VertField<double>>>
    vertexVelocity(const QuadtreeMesh& mesh, HydroVelocityTimeLevel timeLevel,
        bool explicitBoundaryPatches, bool radialAxisSymmetry = false) const;
    std::function<Vector2D(Vector2D)> boundaryAdvectionVelocity(const QuadtreeMesh& mesh,
        HydroVelocityTimeLevel timeLevel, const std::vector<HydroBoundaryPatch>& patches,
        bool radialAxisSymmetry = false) const;
    // Exactly 0.5*(U^(n-1)+U^n) on native faces, then the original spatial
    // sampler. No wall extension; only valid after a completed Hydro substep.
    StaggeredVelocitySampler averagedNativeVelocity(const QuadtreeMesh& mesh) const;
    // C1, face-moment-preserving snapshot of the named native layer. Explicit
    // closed uniform XY support; no time averaging or state mutation.
    HydroSplineVelocity splineVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        HydroVelocityTimeLevel timeLevel) const;
    double completePhaseOutwardFlux(const QuadtreeMesh& mesh,
        const LevelSetHermiteGeometry::Snapshot& interfaceGeometry, bool axisymmetric) const;
    // Publish provenance when external forces are prepared for this managed
    // mesh. Requires initialization; does not allocate/copy a mesh or a field.
    std::shared_ptr<const HydroForceSource> captureForceSource(double sourceTime,
        HydroForceEvolution evolution = HydroForceEvolution::TimeDependent,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry = {}) const;
    // Captures the actual parameter candidate identity and its OWN clock.
    // Fields assembled from that candidate remain caller-owned and immutable.
    std::shared_ptr<const HydroForceSource> captureParametricForceSource(const QuadtreeMesh&,
        const NativeLeafDomain&,const InterfaceState&) const;
    // Once, before any step/cache history. Validates all operands/diagnostics
    // before publishing U0 to current/pre/mid and the actual initial p.
    // Does not project velocity, set a pressure gauge or advance the clock.
    void initializeState(QuadtreeMesh& mesh, const HydroStepGeometry& geometry,
        const HydroInitialFields& initial = {});
    // Generate the original zero/smooth native field for a complete XY/RZ
    // domain, then initializeState once with the supplied pressure (null: zero). Boundary faces
    // start at zero; boundary compatibility projection remains a later stage.
    // Amplitude must be finite and nonnegative, including for the zero profile.
    void initializeCompleteVelocity(QuadtreeMesh& mesh, const HydroStepGeometry& geometry,
        HydroInitialVelocity profile, double amplitude, const CellField<double>* initialPressure = nullptr);
    // Existing material kernels with explicit geometry/time provenance.
    // Numerical failure stops the instance; no automatic rollback or retry.
    void updateMaterials(QuadtreeMesh& mesh, const HydroMaterialUpdate& input);
    // Owned immutable snapshot using the existing native reconstruction.
    // Current and Previous refer to the module's actual CN/AB2 time layers;
    // no field copy is added beyond the original sampler's owned snapshot.
    HydroFittedVelocitySampler fittedVelocity(const QuadtreeMesh& mesh,
        const HydroFittedDomain& domain, HydroVelocityTimeLevel timeLevel,
        HydroFittedVelocitySampler::Reconstruction reconstruction) const;
    // Actual pressure-boundary negative-phase volume rate in code units.
    // Reuses the original strict phase certificate and native face integral.
    long double fittedPhaseOutwardFlux(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const LevelSetHermiteGeometry::Snapshot& interfaceGeometry) const;
    // Explicit zero-time, mass-weighted minimum-change projection of initial
    // native velocities on a complete uniform planar/RZ mesh, with patches or
    // a closed zero-normal box. The closed multiplier has the cell-zero gauge.
    // Imposes prescribed normal velocities before projecting all free faces;
    // outlets remain free and every FV cell continuity equation is retained.
    // Only before the first accepted step. Success publishes the same U0 to
    // current/pre/mid and refreshes divU/uMag. Failure publishes nothing.
    // Physical pressure, time-step/mode locks and pressure caches are untouched.
    template<typename Mesh>
    HydroVelocityInitializationReport initializeCompatibleVelocity(Mesh& mesh,
        const HydroSolveOptions& options, bool axisymmetric = false);
    // Prepare ONLY a zero-time velocity candidate from this accepted complete
    // conservative state and the actual target C lineage. Reuses checked
    // SuperLU, physical boundary data and the native FV incidence; jointly
    // preserves component momentum and continuity. No target Hydro owner,
    // physical pressure, history, material publication or next step is created.
    // Source fields/reports/caches remain unchanged on success and failure.
    Quantum::Hydrodynamics::HydroConservativeTransport::ProjectedMomentumRemap
        projectConservativeMomentumRemap(const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
            const InterfaceConservativePhase& source,const QuadtreeMesh& targetMesh,
            const NativeLeafDomain& targetDomain,const InterfaceConservativePhase& target,
            double liquidDensity,double gasDensity);
    // Closed complete XY/RZ, immediately after an accepted conservative step.
    // Both Interface time layers must have their actual zero-time C/geometry
    // transfers on the target. Project U^n and U^(n-1) with their respective
    // masses; retain the physical pressure by parent-cell injection, preserve
    // the clock/method/dt, and rebuild materials from target C. The projection
    // multipliers never replace physical pressure. Predictor/force/cache data
    // are scratch and rebuilt. Returns an independent, unaccepted AMR owner;
    // caller publishes it jointly with mesh and BOTH Interface handles only
    // after all preparation succeeds. No mutation of the source on any exit.
    std::unique_ptr<Hydrodynamics> prepareConservativeMeshTransfer(
        const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
        const InterfaceState& current,const InterfaceState& previous,
        QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,
        const InterfaceState& targetCurrent,const InterfaceState& targetPrevious,
        double liquidDensity,double gasDensity,double liquidViscosity,double gasViscosity,
        Quantum::Hydrodynamics::HydroConservativeTransport::StateRemapReport* report=nullptr);
    // Historical nonincremental nodal coupling and its cache remain a separate
    // compatibility capability. They cannot mix with initialized managed state.
    // Complete XY/RZ and fitted fixed-dt callers use advance/prepareStep above,
    // with explicit initialization and time-labelled materials/loads.
    template<typename Mesh>
    void updateProperties(Mesh& mesh, LevelSet& levelSet);
    template<typename Mesh>
    void interpolationVelToNode(Mesh& mesh, const int levelSetI, LegacyNodalState& sharpInterface);
    template<typename Mesh>
    void solve(Mesh& mesh, const double& dt, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void solveRZ(Mesh& mesh, const double& dt, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);

    inline double maxDivURes() const { return maxDivURes_; }
    inline double uMagMax() const { return uMagMax_; }

    template<typename Mesh>
    void writeCache(const Mesh& mesh, const std::string& filePath) const;
    template<typename Mesh>
    void readCache(Mesh& mesh, const std::string& filePath);

private:
    friend class HydroFittedMeshTransfer;
    friend class HydroStepCandidate;
    friend class HydroStepTransaction;
    HydroSolveReport advanceStep(QuadtreeMesh& mesh, const HydroStepInput& input, bool deferAcceptance);
    void discardPendingStep() noexcept;
    std::shared_ptr<HydroStepCandidate::Control> pendingStep_;
    std::shared_ptr<HydroStepTransaction::Control> stepTransaction_;
    void discardStepTransaction() noexcept;
    void validateStepTransactionInput() const;
    void validateCompleteVelocityQuery(const QuadtreeMesh& mesh,
        HydroVelocityTimeLevel timeLevel) const;
    enum class StepEntry { Unset, Legacy, Managed };
    StepEntry stepEntry_ = StepEntry::Unset;
    // Bind on kernel entry, including failures, before changing reports/fields.
    // Rejected cross-family calls leave the last valid report/caches untouched.
    std::optional<HydroSolveReport> claimStepEntry(StepEntry entry);
    HydroStepProgress stepProgress_;
    bool initialStateSet_ = false;
    bool conservativeTransportActive_ = false;
    std::optional<HydroMaterialState> materialState_;
    void validateStateGeometry(const QuadtreeMesh& mesh, const HydroStepGeometry& geometry) const;
    void requireLegacyMaterialEntry() const;
    template<typename Mesh>
    void assignConstantMaterialFields(Mesh& mesh, double density, double viscosity);
    template<typename Mesh>
    void assignLevelSetMaterialFields(Mesh& mesh, const LevelSet& levelSet, double densityNegative,
        double densityPositive, double viscosityNegative, double viscosityPositive, double thickness);
    std::shared_ptr<const NativeMeshIdentity> stepMeshIdentity_;
    // Definition (and Eigen/SuperLU types) stays in the implementation unit.
    // shared_ptr supports this incomplete type in the templated constructor;
    // existing unique_ptr members keep Hydrodynamics noncopyable.
    struct LinearSolveCache;
    std::shared_ptr<LinearSolveCache> momentumLinearCache_;
    std::shared_ptr<LinearSolveCache> pressureLinearCache_;
    std::shared_ptr<LinearSolveCache> coupledLinearCache_;
    std::shared_ptr<const HydroFittedDomain> fittedDomain_;
    std::shared_ptr<const HydroFittedMomentumLayout> fittedMomentumLayout_;
    enum class PressureMode { Unset, Nonincremental, Incremental };
    PressureMode pressureMode_ = PressureMode::Unset;
    // Transient projection potential, never the accumulated physical pressure.
    CellField<double> pressurePotential_;
    using PressureFaceData = HydroNativePressureStencil::FaceData;
    template<typename Mesh>
    PressureFaceData pressureFaceData(Mesh& mesh, const Edge& edge, bool horizontal,
        const CellField<double>& pressure) const;
    template<typename Mesh>
    void addHistoricalPressure(Mesh& mesh, bool axisymmetric, Eigen::VectorXd& rhs);
    template<typename Mesh>
    void addFaceForce(Mesh& mesh, bool axisymmetric,
        const VEdgeField<double>& forceX, const HEdgeField<double>& forceY, Eigen::VectorXd& rhs);
    template<typename Mesh>
    void correctCoupledBalance(Mesh& mesh, bool axisymmetric,
        const Eigen::SparseMatrix<double>& matrix, const Eigen::VectorXd& momentumRhs);
    // Explicit fixed-grid boundary data. Interior rows keep their existing
    // CN/AB2 assembly; only boundary stresses/fluxes are completed here.
    struct BoundaryFaceData
    {
        int owner = -1, inner = -1;
        bool horizontal = false, outlet = false;
        double sign = 0.0, value = 0.0, area = 0.0, innerArea = 0.0;
        double volume = 0.0, normalSpacing = 0.0, tangentSpacing = 0.0;
    };
    std::vector<BoundaryFaceData> boundaryFaces_;
    std::vector<HydroBoundaryPatch> frozenBoundaryPatches_;
    bool boundaryConfigurationFixed_ = false;
    template<typename Mesh> std::vector<BoundaryFaceData> makeBoundaryFaces(Mesh& mesh,
        const std::vector<HydroBoundaryPatch>& patches, bool axisymmetric, bool includeClosed = false) const;
    template<typename Mesh> HydroBoundaryProjection assembleBoundaryProjection(Mesh& mesh,
        double inverseMassScale, bool axisymmetric, const std::vector<BoundaryFaceData>& boundaries,
        const VEdgeField<double>& vx, const HEdgeField<double>& vy) const;
    template<typename Mesh> void prepareBoundaryPatches(Mesh& mesh, bool axisymmetric);
    template<typename Mesh> void addBoundaryMomentum(Mesh& mesh, double dt, bool axisymmetric, HydroMomentumAssembly& assembly,
        bool includeAdvection=true);
    void solveConservativeNativeStep(QuadtreeMesh& mesh,double dt,bool axisymmetric,
        const Quantum::Hydrodynamics::HydroConservativeTransport& transport,
        const VEdgeField<double>* faceX,const HEdgeField<double>* faceY,const HydroMidpointData* midpoint);
    void addConservativeMomentum(const QuadtreeMesh& mesh,double dt,bool axisymmetric,
        const Quantum::Hydrodynamics::HydroConservativeTransport& transport,HydroMomentumAssembly& assembly);
    Quantum::Hydrodynamics::HydroConservativeTransport::ProjectedMomentumRemap
        projectConservativeMomentumLayer(const QuadtreeMesh&,const NativeLeafDomain&,
            const InterfaceConservativePhase&,const std::vector<double>&,
            const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceConservativePhase&,double,double);
    template<typename Mesh> void addNativeRZHoop(Mesh& mesh, double dt, HydroMomentumAssembly& assembly);
    template<typename Mesh> HydroPressureSystem assembleNativeAdaptivePressure(Mesh& mesh, double dt, bool axisymmetric);
    template<typename Mesh> void projectBoundaryVelocity(Mesh& mesh, double dt, bool axisymmetric);
    template<typename Mesh> void checkOutletVelocity(Mesh& mesh, double dt, bool withPredictor = false) const;

    std::unique_ptr<virtualSupCell[]> supCells_;
    std::unique_ptr<virtualSupCell[]> supSonCells_;

    double maxDivURes_ = 0.0;
    double uMagMax_ = 0.0;

    double iterCoeff_ = 2.3;
    double fixedDt_ = 0.0;
    HydroSolveOptions activeOptions_;
    HydroSolveReport activeReport_;

    template<typename Mesh>
    HydroSolveReport solveStep(Mesh& mesh, double dt, bool axisymmetric,
        const HydroSolveOptions& options, LegacyNodalState* sharpInterface,
        ImmersedGeometry<Mesh>* iGeometry, bool fixedTimeStep,
        const CellField<double>* forceX = nullptr, const CellField<double>* forceY = nullptr,
        const VEdgeField<double>* faceForceX = nullptr, const HEdgeField<double>* faceForceY = nullptr,
        StepEntry entry = StepEntry::Legacy,
        const Quantum::Hydrodynamics::HydroConservativeTransport* conservativeTransport = nullptr,
        bool variableConservativeStep = false,const HydroMidpointData* midpoint = nullptr);
    HydroSolveReport solveFittedStep(QuadtreeMesh& mesh, double dt,
        const HydroSolveOptions& options, const HydroFittedDomain& domain,
        const CellField<double>* forceX, const CellField<double>* forceY,
        const VEdgeField<double>* faceForceX, const HEdgeField<double>* faceForceY);
    template<typename Mesh>
    std::pair<double, double> velocityDiagnostics(Mesh& mesh,
        const VEdgeField<double>& vx, const HEdgeField<double>& vy,
        bool axisymmetric, bool storeFields, CellField<double>* divergenceTarget = nullptr,
        CellField<double>* magnitudeTarget = nullptr);
    Eigen::VectorXd checkedLinearSolve(const Eigen::SparseMatrix<double>& matrix,
        const Eigen::VectorXd& rhs, const std::string& stage, double& residual,
        std::optional<double> relativeTolerance = std::nullopt,int maximumRefinements = 0);

    template<typename Mesh>
    void searchSupCoeffs(Mesh& mesh);

    // template<typename Mesh, typename EdgeField>
    // double checkLenthInterpolation(Mesh& mesh, const EdgeField& edgeField, const Vert& vert, const bool& isHorizontal);
    template<bool Axisymmetric, typename Mesh>
    void buildNativeMomentumMatrix(Mesh& mesh, const double& dt, HydroMomentumAssembly& assembly,
        const std::vector<unsigned char>* selectedRows = nullptr);
    template<bool Axisymmetric, typename Mesh>
    void buildNativeMomentumRHS(Mesh& mesh, const double& dt, HydroMomentumAssembly& assembly,
        const std::vector<unsigned char>* selectedRows = nullptr,
        const CellField<double>* bodyX = nullptr, const CellField<double>* bodyY = nullptr,
        bool includeAdvection=true);
    template<bool Axisymmetric, typename Mesh>
    HydroPressureSystem assembleCellPressure(Mesh& mesh, const double& dt);
    template<bool Axisymmetric, typename Mesh>
    void correctNativeVelocity(Mesh& mesh, const double& dt);

    
    template<typename Mesh>
    void calculateExternalForce(Mesh& mesh, const int levelSetI, LegacyNodalState& sharpInterface);
    template<typename Mesh>
    void applyEvapVel(Mesh& mesh, LegacyNodalState& sharpInterface);
    template<typename Mesh>
    void applyDirectForce(Mesh& mesh, const double& dt, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
};

#include "Modules/Hydrodynamics/Hydrodynamics_impl.h"
#include "Modules/Hydrodynamics/HydrodynamicsMomentum_impl.h"
#include "Modules/Hydrodynamics/HydrodynamicsPressure_impl.h"
#include "Modules/Hydrodynamics/HydrodynamicsCoupled_impl.h"
#include "Modules/Hydrodynamics/HydroBoundary_impl.h"
#include "Modules/Hydrodynamics/HydroInitialization_impl.h"
