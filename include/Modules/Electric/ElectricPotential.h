#pragma once

#include "Core/Base/Vector2D.h"
#include "Core/Base/TimeLayer.h"
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>
#include <Eigen/SparseCore>

class QuadtreeMesh;
namespace LAMP::Interface { class HermiteGalsPublication; }

namespace LAMP::Electric
{
class ElectricPhaseGeometry;
class SurfaceChargeField;
class ElectricPanelChargeLoad;
class ElectricPanelCurrentOperator;
class ElectricFixedGeometryChargeSolver;
class ElectricMovingGeometryChargeStep;
class SurfaceChargeMaterialTransfer;
struct SurfaceChargeFieldStep;
struct SurfaceChargeStepReport;
class ElectricSurfaceChargeTimeState;
class ElectricFixedGeometryChargeStep;
// P0 sharp-interface electrostatics numerical building block on the actual complete, nonperiodic
// UniformLeaves Cartesian mesh. It does not select an App, rebuild interface
// geometry, supply a charge evolution law, or convert old FV potential data.
struct ElectricCartesianBoundary
{
    enum class Type { Dirichlet, Neumann };
    Type type = Type::Dirichlet;
    // Dirichlet V, or outward epsilon*grad(V).n for Neumann.
    std::function<double(Vector2D)> value;
};
struct ElectricCartesianProblem
{
    bool axisymmetric = false; // m=1 (XY), m=2*pi*r (RZ), r=x.
    std::vector<double> epsilonByCell; // finite positive, actual native cell indices.
    // Strong source: -div(m*epsilon*grad(V))=m*f. Empty means f=0.
    std::function<double(Vector2D)> source;
    // Actual outer sides in W,S,E,N order; no inferred periodic or wall law.
    std::array<ElectricCartesianBoundary,4> boundaries;
};
struct ElectricCartesianControls
{
    double relativeLinearTolerance = 1e-11;
    double relativeSymmetryTolerance = 1e-12;
    // Independent physical conservation target; cannot exceed the E-R5 bound.
    double relativeConservationTolerance = 1e-10;
};
// Passive elapsed wall times in seconds. Fields do not enter numerical gates.
struct ElectricCartesianTimingReport
{
    double boundarySeconds = 0, localOperatorsSeconds = 0, localLoadsSeconds = 0;
    double assemblySeconds = 0, factorSeconds = 0, solveSeconds = 0;
    double recoverySeconds = 0, gaussAffineSeconds = 0;
};
struct ElectricCartesianReport
{
    ElectricCartesianTimingReport timing;
    std::size_t cells = 0, physicalFaces = 0, axisFaces = 0;
    int freeSkeletonDofs = 0, prescribedSkeletonDofs = 0;
    double relativeLinearResidual = 0, relativeSymmetryDefect = 0;
    // Zero only when no free skeleton DOFs require a global factorization.
    double minimumSkeletonPivot = 0;
    double maxRelativeCellResidual = 0;
    double maxRelativeFreeFaceBalance = 0;
};
class ElectricCartesianResult
{
public:
    using Matrix = Eigen::SparseMatrix<double>;
    const ElectricCartesianReport& report() const;
    bool axisymmetric() const;
    void validateMesh(const QuadtreeMesh&) const;
    std::size_t nativeCellCount() const;
    // Basis is {1,X,Y,X^2,XY,Y^2,X^3,X^2Y,XY^2,Y^3}, with
    // X=(x-center.x)/scale, Y=(y-center.y)/scale; P2 uses first six.
    Vector2D cellCenter(int nativeCellIndex) const;
    double cellScale(int nativeCellIndex) const;
    std::array<double,6> cellPolynomialP2(int nativeCellIndex) const;
    std::array<double,10> reconstructedPotentialP3(int nativeCellIndex) const;
    double potentialInCell(int nativeCellIndex,Vector2D point) const;
    Vector2D gradientInCell(int nativeCellIndex,Vector2D point) const;
    // Native edge index in its actual horizontal/vertical family, never label.
    // One shared trace per positive-measure physical face, including prescribed
    // Dirichlet faces. P2 Legendre basis uses zeta increasing x/y on horizontal/
    // vertical faces. Axis queries reject: r=0 has no trace unknown.
    std::array<Vector2D,2> nativeFaceBounds(bool horizontal,int nativeEdgeIndex) const;
    std::array<double,3> nativeFaceTrace(bool horizontal,int nativeEdgeIndex) const;
    // Full local K_fc*Vc+K_ff*Vf minus the local Neumann face load.
    // W,S,E,N side; axis has no face row and rejects. Interior sums vanish.
    std::array<double,3> outwardFaceResidual(int nativeCellIndex,int side) const;
    // Only FREE shared face modes survive; Dirichlet data were eliminated.
    // These read-only matrices belong to this result and contain no cell DOFs.
    const Matrix& skeletonMatrix() const;
    const std::vector<double>& skeletonRhs() const;
    const std::vector<double>& skeletonSolution() const;
private:
    struct Data;
    explicit ElectricCartesianResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend ElectricCartesianResult solveElectricPotential(
        const QuadtreeMesh&,const ElectricCartesianProblem&,const ElectricCartesianControls&);
};
// Eight-point physical Gauss rules, P2 weighted face projection, local P3
// weighted Green reconstruction and symmetric local stabilization. Only cell
// DOFs are condensed. Anchored free skeleton is solved by SPD LDLT, without a
// fallback or any nonlocal P3 correction. Axis requires homogeneous Neumann.
ElectricCartesianResult solveElectricPotential(const QuadtreeMesh&,
    const ElectricCartesianProblem&,const ElectricCartesianControls& = {});
struct ElectricSurfaceChargeMomentControls
{
    double relativeQuadratureTolerance = 1e-9;
    int maximumAdditionalRefinements = 4;
};
struct ElectricSurfaceChargeMomentReport
{
    std::size_t panels = 0;
    int additionalRefinements = 0;
    double maxRelativeMomentDifference = 0;
    double maxRelativePanelLengthDifference = 0;
};
// Geometry-bound linear P_Q, not another charge state. Fields are supplied in
// publication component order; flat Q uses that order and each native chain's
// actual material-panel/graph order. The existing Field owns no physical clock:
// coevalChargeLayer is an explicit caller declaration, checked field-for-field
// against this publication, rather than an independent history certification.
class ElectricSurfaceChargeMomentOperator
{
public:
    using Matrix = Eigen::SparseMatrix<double>;
    static ElectricSurfaceChargeMomentOperator prepare(const ElectricPhaseGeometry&,
        const std::vector<const SurfaceChargeField*>& fieldsByComponent,
        const ::LAMP::TimeLayer& coevalChargeLayer,
        const ElectricSurfaceChargeMomentControls& = {});
    const ElectricPhaseGeometry& geometry() const;
    ::LAMP::TimeLayer timeLayer() const;
    const ElectricSurfaceChargeMomentReport& report() const;
    std::size_t panelCount() const;
    int interfaceFace(std::size_t panel) const;
    int component(std::size_t panel) const;
    int graph(std::size_t panel) const;
    // P_Q: (3*N_panel) x N_panel; constant rows are exact identities.
    const Matrix& projectionMatrix() const;
    // B_Gamma: (3*N_shared_face) x (3*N_panel), before BC elimination.
    const Matrix& interfaceInjectionMatrix() const;
    ElectricPanelChargeLoad apply(const std::vector<double>& integratedCharge) const;
private:
    struct Data;
    static ElectricSurfaceChargeMomentOperator prepareRetained(const ElectricPhaseGeometry&,
        const std::vector<const SurfaceChargeField*>&,const ::LAMP::TimeLayer&,
        const ElectricSurfaceChargeMomentControls&,std::shared_ptr<const void>);
    explicit ElectricSurfaceChargeMomentOperator(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend class ElectricPanelChargeLoad;
    friend class ElectricSurfaceChargeTimeState;
};
// Immutable candidate integrals and their derived moments on the bound
// geometry. Applying trial Q neither publishes it to Field nor advances time.
class ElectricPanelChargeLoad
{
public:
    const ElectricSurfaceChargeMomentOperator& momentOperator() const;
    const std::vector<double>& panelCharge() const;
    const std::vector<std::array<double,3>>& moments() const;
    const std::vector<std::array<double,3>>& faceMoments() const;
    void validateGeometry(const ElectricPhaseGeometry&) const;
private:
    struct Data;
    explicit ElectricPanelChargeLoad(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend class ElectricSurfaceChargeMomentOperator;
};
// P1 sharp dielectric problem on the one prepared Hermite phase geometry.
// Indices [0,1] mean negative-phi/liquid and positive-phi/gas respectively.
// The three-argument P1 solve has no surface charge. The explicit P2 overload
// consumes a bound panel-charge load on the same shared P2 Gamma trace.
struct ElectricPhaseProblem
{
    std::array<double,2> epsilonByPhase{{1.0,1.0}};
    std::array<std::function<double(Vector2D)>,2> sourceByPhase;
    std::array<ElectricCartesianBoundary,4> boundaries;
};
struct ElectricPhaseReport
{
    ElectricCartesianReport skeleton;
    std::size_t interfaceFaces = 0;
    double maxRelativeInterfaceBalance = 0;
    double maxRelativeLocalGaussResidual = 0;
    double globalGaussResidual = 0;
    double relativeGlobalGaussResidual = 0;
};
class ElectricPhaseResult
{
public:
    using Matrix = Eigen::SparseMatrix<double>;
    const ElectricPhaseGeometry& geometry() const;
    const ElectricPhaseReport& report() const;
    // Actual material coefficients retained by this immutable solve.
    const std::array<double,2>& epsilonByPhase() const;
    // Actual electric solution clock; fixed geometry retains its separate source clock.
    ::LAMP::TimeLayer timeLayer() const;
    void validateSource(const QuadtreeMesh&,
        const LAMP::Interface::HermiteGalsPublication&) const;
    bool axisymmetric() const;
    std::size_t phaseElementCount() const;
    int nativeCellIndex(std::size_t element) const;
    bool liquid(std::size_t element) const;
    std::vector<std::size_t> phaseElements(int nativeCellIndex,bool liquid) const;
    Vector2D elementCenter(std::size_t element) const;
    double elementScale(std::size_t element) const;
    std::array<double,6> cellPolynomialP2(std::size_t element) const;
    std::array<double,10> reconstructedPotentialP3(std::size_t element) const;
    // Explicit one-sided element polynomial query within its closed native
    // cell. The caller selects the phase element; no material relocation is
    // performed. Same-source Gamma points support both one-sided traces.
    double potentialInElement(std::size_t element,Vector2D point) const;
    Vector2D gradientInElement(std::size_t element,Vector2D point) const;
    std::array<double,3> skeletalFaceTrace(std::size_t face) const;
    // Local face order is geometry.elements()[element].faces; values are
    // K_fc*Vc+K_ff*Vf minus the local physical Neumann load. Raw reconstructed
    // gradients are not substituted for these equilibrated residual moments.
    std::array<double,3> outwardFaceResidual(std::size_t element,std::size_t localFace) const;
    // Physical D.out restores the prescribed Neumann load before negating
    // the primal flux residual. Both Gamma D traces use n=negative->positive.
    std::array<double,3> physicalOutwardDisplacementMoments(std::size_t element,std::size_t localFace) const;
    std::array<double,3> normalDisplacementMoments(std::size_t interfaceFace,bool liquid) const;
    std::array<double,3> normalDisplacementPolynomialP2(std::size_t interfaceFace,bool liquid) const;
    std::array<double,3> interfaceChargeMoments(std::size_t interfaceFace) const;
    const ElectricPanelChargeLoad* surfaceChargeLoad() const;
    // Actual solve permittivities are retained; only nonnegative sigma is supplied.
    // C columns use precisely skeletonSolution()'s free trace ordering.
    ElectricPanelCurrentOperator preparePanelCurrent(const ElectricSurfaceChargeMomentOperator&,
        const std::array<double,2>& conductivityByPhase) const;

    const Matrix& skeletonMatrix() const;
    const std::vector<double>& skeletonRhs() const;
    const std::vector<double>& skeletonSolution() const;
private:
    struct Data;
    explicit ElectricPhaseResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    ElectricPanelCurrentOperator preparePanelCurrentBound(const ElectricSurfaceChargeMomentOperator&,
        const std::array<double,2>&,const ElectricPanelCurrentOperator*) const;
    static ElectricPhaseResult recoverPrepared(const ElectricPhaseResult&,const ElectricPanelChargeLoad&,
        const std::vector<double>&,const ::LAMP::TimeLayer&);
    static ElectricPhaseResult solvePrepared(const ElectricPhaseGeometry&,const ElectricPhaseProblem&,
        const ElectricPanelChargeLoad*,const ElectricCartesianControls&);
    friend class ElectricFixedGeometryChargeSolver;
    friend ElectricPhaseResult solveElectricPhaseElectrostatics(
        const ElectricPhaseGeometry&,const ElectricPhaseProblem&,const ElectricCartesianControls&);
    friend ElectricPhaseResult solveElectricPhaseElectrostatics(
        const ElectricPhaseGeometry&,const ElectricPhaseProblem&,const ElectricPanelChargeLoad&,const ElectricCartesianControls&);
};
struct ElectricPanelCurrentReport
{
    std::size_t panels = 0;
    double totalPanelCurrent = 0;
    double exteriorCurrent = 0;
    double bulkCurrent = 0;
    // Sum(I_i) + exterior J.out - sum(rho_phase*cellLoad_0).
    double globalCurrentBalance = 0;
    double relativeGlobalCurrentBalance = 0;
    double maxRelativeAffineResidual = 0;
};
// Immutable Ohmic functional of the existing equilibrated displacement flux.
// It owns no charge state and performs no potential solve or time advance.
class ElectricPanelCurrentOperator
{
public:
    using Matrix = Eigen::SparseMatrix<double>;
    const ElectricSurfaceChargeMomentOperator& momentOperator() const;
    const ElectricPhaseGeometry& geometry() const;
    ::LAMP::TimeLayer timeLayer() const;
    const ElectricPanelCurrentReport& report() const;
    const std::array<double,2>& epsilonByPhase() const;
    const std::array<double,2>& conductivityByPhase() const;
    // Native material panel order is exactly momentOperator()'s flat Q order.
    const Matrix& matrix() const; // C: N_panel x N_free_skeleton.
    const std::vector<double>& offset() const; // c_D: prescribed trace and volume source.
    const std::vector<double>& panelCurrent() const; // Direct original residual functional.
    std::vector<double> apply(const std::vector<double>& freeSkeleton) const;
    // Both one-sided currents use n=negative->positive; no second E is fitted.
    std::array<double,3> normalCurrentMoments(std::size_t panel,bool liquid) const;
    std::array<double,3> normalCurrentPolynomialP2(std::size_t panel,bool liquid) const;
private:
    struct Data;
    explicit ElectricPanelCurrentOperator(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend class ElectricPhaseResult;
};

enum class ElectricChargeTimeMethod { BackwardEuler, VariableBdf2 };
struct ElectricChargeTimeControls
{
    ElectricChargeTimeMethod method = ElectricChargeTimeMethod::VariableBdf2;
    double relativeSchurTolerance = 1e-11;
    double relativeChargeBalanceTolerance = 1e-10;
    int maximumIterations = 200, restartLength = 30;
    double maximumStepRatio = 2.0;
};
// schurElectricSolveSeconds is included in schurSeconds, not an extra stage.
struct ElectricChargeTimingReport
{
    double totalSeconds = 0, schurSeconds = 0, schurElectricSolveSeconds = 0;
    double momentLoadSeconds = 0, finalElectricSolveSeconds = 0;
    double recoverySeconds = 0, currentSeconds = 0, balanceSeconds = 0;
};
struct ElectricChargeTimeReport
{
    ElectricChargeTimingReport timing;
    int temporalOrder = 0, iterations = 0, electricSolves = 0;
    double a0 = 0, a1 = 0, a2 = 0, dt = 0, previousDt = 0;
    double relativeSchurResidual = 0;
    double maxRelativePanelEquationResidual = 0;
    double integratedEquationResidual = 0, relativeTotalEquationResidual = 0;
};
// A candidate owns read-only Q/V/current at one caller-supplied physical clock.
// On the fixed path geometry/P_Q retain their original publication clock.
// A moving transaction supplies the actual target publication at this physical
// clock, then uses the same frozen-geometry solver for its current update.
class ElectricFixedGeometryChargeStep
{
public:
    const std::vector<double>& integratedCharge() const;
    const ElectricPhaseResult& potential() const;
    const ElectricPanelCurrentOperator& current() const;
    ::LAMP::TimeLayer timeLayer() const;
    ::LAMP::TimeLayer geometrySourceTimeLayer() const;
    const ElectricChargeTimeReport& report() const;
private:
    struct Data;
    explicit ElectricFixedGeometryChargeStep(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend class ElectricFixedGeometryChargeSolver;
    friend class ElectricSurfaceChargeTimeState;
};
// Fixed geometry/material/BC/source prepared once. Each matrix-free Schur
// action uses the same SPD electric factor and the actual sparse P_Q/B_Gamma/C.
// It owns no mutable charge and never accepts a candidate on its own.
class ElectricFixedGeometryChargeSolver
{
public:
    static ElectricFixedGeometryChargeSolver prepare(const ElectricPhaseResult&,
        const ElectricSurfaceChargeMomentOperator&,const std::array<double,2>& conductivityByPhase);
    ::LAMP::TimeLayer geometrySourceTimeLayer() const;
    const ElectricSurfaceChargeMomentOperator& momentOperator() const;
    // Read-only initial current already prepared with this fixed geometry/material.
    const ElectricPanelCurrentOperator& initialCurrent() const;
    ElectricFixedGeometryChargeStep advance(const ElectricSurfaceChargeTimeState&,
        const ::LAMP::TimeLayer&,const ElectricChargeTimeControls& = {}) const;
private:
    struct Data;
    explicit ElectricFixedGeometryChargeSolver(std::shared_ptr<const Data>);
    const Data& checked() const;
    ElectricFixedGeometryChargeStep advanceCandidate(const ElectricSurfaceChargeTimeState&,
        const ::LAMP::TimeLayer&,const ElectricChargeTimeControls&,bool midpoint) const;
    std::shared_ptr<const Data> data_;
    friend class ElectricSurfaceChargeTimeState;
};
// Immutable moving candidate. Historical physical clocks stay unchanged while
// mapped integrals use target geometry. All target Fields are retained by its
// P_Q, so copies of the resulting solver/potential/current stay valid. A Midpoint
// stage is read-only and cannot be committed, including through its inner step.
class ElectricMovingGeometryChargeStep
{
public:
    const ElectricFixedGeometryChargeStep& step() const;
    const std::vector<double>& mappedIntegratedCharge() const;
    const std::vector<double>* mappedPreviousIntegratedCharge() const;
    ::LAMP::TimeLayer sourceTimeLayer() const;
    ::LAMP::TimeLayer sourceGeometryTimeLayer() const;
    const ::LAMP::TimeLayer* previousSourceTimeLayer() const;
    ::LAMP::TimeLayer geometrySourceTimeLayer() const;
    std::size_t componentCount() const;
    // Retains the actual target geometry owner. After acceptance its Q may be
    // updated by that state; the candidate's own numerical receipts stay fixed.
    const SurfaceChargeField& field(std::size_t component) const;
    const std::vector<SurfaceChargeFieldStep>& remapReports() const;
    const std::vector<std::vector<SurfaceChargeStepReport>>& historyRemapReports() const;
private:
    struct Data;
    explicit ElectricMovingGeometryChargeStep(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
    friend class ElectricSurfaceChargeTimeState;
};
// Initially borrows the unique mutable component Fields; successful moving
// commits retain the new component owners. Physical/history and representation
// clocks are separate immutable receipts; outside Field edits invalidate this
// binding. All component charges/history publish atomically.
class ElectricSurfaceChargeTimeState
{
public:
    static ElectricSurfaceChargeTimeState initialize(const ElectricFixedGeometryChargeSolver&,
        const std::vector<SurfaceChargeField*>&,const ::LAMP::TimeLayer&);
    ~ElectricSurfaceChargeTimeState();
    ElectricSurfaceChargeTimeState(ElectricSurfaceChargeTimeState&&) noexcept;
    ElectricSurfaceChargeTimeState& operator=(ElectricSurfaceChargeTimeState&&) noexcept;
    ElectricSurfaceChargeTimeState(const ElectricSurfaceChargeTimeState&) = delete;
    ElectricSurfaceChargeTimeState& operator=(const ElectricSurfaceChargeTimeState&) = delete;
    ::LAMP::TimeLayer timeLayer() const;
    const ::LAMP::TimeLayer* previousTimeLayer() const;
    const std::vector<double>& integratedCharge() const;
    const std::vector<double>* previousIntegratedCharge() const;
    void commit(const ElectricFixedGeometryChargeStep&);
    const ElectricFixedGeometryChargeSolver& solver() const;
    ::LAMP::TimeLayer geometrySourceTimeLayer() const;
    std::size_t componentCount() const;
    const SurfaceChargeField& field(std::size_t component) const;
    // Geometry owner validates completed-native transfer provenance before this
    // synchronous call. Electric binds its real source/target panels and
    // interval to this accepted state and the target publication's full clock.
    // No external history arrays are accepted. Same map for Qn/Qn-1; only the
    // successful commit replaces Fields, preparation and physical receipts.
    ElectricMovingGeometryChargeStep prepareMoving(const QuadtreeMesh&,
        const ElectricPhaseGeometry&,const std::vector<SurfaceChargeMaterialTransfer>&,
        const ::LAMP::TimeLayer&,const ElectricPhaseProblem&,
        const std::array<double,2>& conductivityByPhase,
        const ElectricSurfaceChargeMomentControls& = {},const ElectricCartesianControls& = {},
        const ElectricChargeTimeControls& = {}) const;
    // Read-only next Midpoint candidate on its actual publication. Its charge
    // increment is target.time()-accepted.time(), not the full clock interval.
    // Both histories retain their physical clocks after the same material map;
    // neither this preparation nor public commit can accept a stage as history.
    ElectricMovingGeometryChargeStep prepareMovingStage(const QuadtreeMesh&,
        const ElectricPhaseGeometry&,const std::vector<SurfaceChargeMaterialTransfer>&,
        const ::LAMP::TimeLayer&,const ElectricPhaseProblem&,
        const std::array<double,2>& conductivityByPhase,
        const ElectricSurfaceChargeMomentControls& = {},const ElectricCartesianControls& = {},
        const ElectricChargeTimeControls& = {}) const;
    void commit(const ElectricMovingGeometryChargeStep&);
private:
    struct Data;
    explicit ElectricSurfaceChargeTimeState(std::unique_ptr<Data>);
    const Data& checked() const;
    ElectricMovingGeometryChargeStep prepareMovingCandidate(const QuadtreeMesh&,
        const ElectricPhaseGeometry&,const std::vector<SurfaceChargeMaterialTransfer>&,
        const ::LAMP::TimeLayer&,const ElectricPhaseProblem&,
        const std::array<double,2>&,const ElectricSurfaceChargeMomentControls&,
        const ElectricCartesianControls&,const ElectricChargeTimeControls&,bool midpoint) const;
    std::unique_ptr<Data> data_;
    friend class ElectricFixedGeometryChargeSolver;
};

ElectricPhaseResult solveElectricPhaseElectrostatics(const ElectricPhaseGeometry&,
    const ElectricPhaseProblem&,const ElectricCartesianControls& = {});
ElectricPhaseResult solveElectricPhaseElectrostatics(const ElectricPhaseGeometry&,
    const ElectricPhaseProblem&,const ElectricPanelChargeLoad&,const ElectricCartesianControls& = {});


}
