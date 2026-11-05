#pragma once
#include "Modules/Interface/GeometricVofPublication.h"
#include "Core/Base/Vector2D.h"
#include <cstdint>
#include <array>
#include <memory>
#include <vector>

namespace LAMP::Electric
{
class ElectricVofAdvectedCharge;
class ElectricVofChargeMeshTransfer;
class ElectricVofField;
class ElectricVofCharge;
class ElectricVofChargeBoundary;
class ElectricVofChargeTransport;
namespace Detail
{
// Private proof for one callback-free Electric public entry. No public entry
// accepts/returns it, and no retained Data stores it. Only a completed original
// Publication/Transport full live validation can establish its actual binding.
class ElectricVofValidationBatch
{
    friend class LAMP::Electric::ElectricVofField;
    friend class LAMP::Electric::ElectricVofCharge;
    friend class LAMP::Electric::ElectricVofChargeBoundary;
    friend class LAMP::Electric::ElectricVofAdvectedCharge;
    friend class LAMP::Electric::ElectricVofChargeTransport;
    ElectricVofValidationBatch()=default;
    ElectricVofValidationBatch(const ElectricVofValidationBatch&)=delete;
    ElectricVofValidationBatch& operator=(const ElectricVofValidationBatch&)=delete;
    void rememberValidated(const QuadtreeMesh&,const NativeLeafDomain&);
    bool sameBinding(const QuadtreeMesh&,const NativeLeafDomain&) const;
    void validateTransport(const QuadtreeMesh&,const ElectricVofCharge&,
        const Interface::GeometricVofTransport&);
    const QuadtreeMesh* mesh_=nullptr;
    std::shared_ptr<const NativeMeshIdentity> identity_;
    const Interface::GeometricVofTransport* transport_=nullptr; // this stack only
};
}
struct ElectricVofMaterials
{
    double epsilonLiquid=1,epsilonGas=1;
    double conductivityLiquid=0,conductivityGas=0;
};
struct ElectricVofBoundary
{
    enum class Type { Interior,Dirichlet,HomogeneousNeumann };
    Type type=Type::Interior;
    double value=0; // Static actual face potential; zero for Interior/Neumann.
};
enum class ElectricVofChargeAdvectionScheme { DonorCell, LimitedLinearFct };
struct ElectricVofControls
{
    double relativeLinearTolerance=1e-11; // Cannot exceed the original Electric bound.
    double relativeConservationTolerance=1e-10;
    // Existing direct module callers retain their qualified low-order default.
    // The A App reader explicitly selects LimitedLinearFct by default.
    ElectricVofChargeAdvectionScheme chargeAdvection=ElectricVofChargeAdvectionScheme::DonorCell;
    // Strict integrated backends share the same operator and physical gates;
    // each selected backend rejects failure without an alternative fallback.
    enum class Backend { Gmres, SuperLU };
    Backend backend=Backend::Gmres;
};
// Integrated signed Eulerian cell charge, not an HG surface charge or density.
// Initial public input must be bound to a genuine accepted VOF mesh/clock.
// Native-cell arrays remain dense; cells outside the actual selected complete-
// leaf domain must carry zero sentinels, never an implicit gas/charge reservoir.
class ElectricVofCharge
{
public:
    static ElectricVofCharge fromIntegrated(const QuadtreeMesh&,
        const Interface::GeometricVofPublication& acceptedVof,const std::vector<double>&);
    std::uint64_t identity() const;
    LAMP::TimeLayer timeLayer() const;
    const std::vector<double>& integratedCharge() const;
    void validateFor(const QuadtreeMesh&,LAMP::TimeLayer) const;
    // Same-stage consumers additionally require the actual retained geometry
    // owner and complete geometry/charge clock, never merely equal array sizes.
    void validateFor(const QuadtreeMesh&,const Interface::GeometricVofPublication&,LAMP::TimeLayer) const;
    bool sameInstance(const ElectricVofCharge&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofCharge(std::shared_ptr<const Data>);
    friend class ElectricVofField;
    friend class ElectricVofChargeMeshTransfer;
    friend class Detail::ElectricVofValidationBatch;
    void validateForInBatch(const QuadtreeMesh&,LAMP::TimeLayer,
        const Detail::ElectricVofValidationBatch&) const;
    void validateSelectedDomain(const NativeLeafDomain&) const;
    static ElectricVofCharge fromPrepared(const QuadtreeMesh&,
        const Interface::GeometricVofPublication& meshWitness,LAMP::TimeLayer,std::vector<double>);
};
struct ElectricVofGradientCoefficient { int nativeCell=-1;double value=0; };
struct ElectricVofFace
{
    int nativeRow=-1,negativeCell=-1,positiveCell=-1;
    double physicalMeasure=0,epsilon=0,conductivity=0,implicitCoefficient=0;
    // H is integrated -grad(phi).n in canonical +r/x or +z/y direction.
    // Geometry coefficients include physical Af exactly once, including RZ.
    // The affine action consumes the exact sum from Field.potentialTerms(),
    // with one RN to H; potential() is its nearest scalar output projection.
    std::vector<ElectricVofGradientCoefficient> gradientCoefficients;
    double gradientConstant=0;
    double integratedField=0,displacement=0,current=0;
};
// Exact factorized receipt, bound to this Field's geometry/source/full clock.
// For each canonical faces() row: r=exact(duration*current),
// t=eta*nearestEvenInteger(r/eta), eta=2^quantumExponent. t may have
// more than 53 significant bits; no public double is its exact projection.
struct ElectricVofOhmicTransferReceipt
{
    enum class Rule { MinimumBinary64QuantumNearestEven };
    Rule rule=Rule::MinimumBinary64QuantumNearestEven;
    double duration=0;
    int quantumExponent=-1074;
};
struct ElectricVofReport
{
    std::size_t cells=0,nativeFaces=0,nonorthogonalFaces=0,axisFaces=0;
    double relativeLinearResidual=0,maximumRelativeGaussResidual=0;
    bool usedAdvection=false;
    bool electrostaticSnapshot=false; // No charge advancement; J is instantaneous only.
    double maximumAffineGradientDefect=0,maximumPhaseCapacityDefect=0;
    double geometrySeconds=0,assembleSeconds=0,factorSeconds=0,solveSeconds=0,totalSeconds=0;
    bool usedSuperLU=false;
    bool usedGMRES=false;
    int linearIterations=0,linearRestarts=0,linearMaximumWorkers=0;
    // Optional qualified post-solve execution only; the PDE is always solved once.
    bool postprocessParallelUsed=false,postprocessParallelFallback=false;
    std::size_t postprocessWorkers=0,postprocessScratchBytes=0;
    // Actual single-pass stage time; parent includes setup and necessary replay.
    bool postprocessTimingAvailable=false;
    double postprocessFaceSeconds=0,postprocessCellSeconds=0,postprocessFinalizeSeconds=0,postprocessSeconds=0;
    bool integratedExpandedResidualAvailable=false,integratedExpandedResidualParallelUsed=false,integratedExpandedResidualParallelFallback=false;
    std::size_t integratedExpandedResidualWorkers=0,integratedExpandedResidualRows=0,integratedExpandedResidualScratchBytes=0;
    bool integratedExpandedResidualTimingAvailable=false;
    double integratedExpandedResidualRouteSeconds=0,integratedExpandedResidualSeconds=0;
};
// Frozen-material Ohmic field substep. The accepted-source overload takes
// Qadv=source; the explicit typed transport overload takes actual stage Qadv.
// Charge advection/remap have their own owners; Maxwell/App remain separate.
class ElectricVofField
{
public:
    enum class CellGradientRecovery { Recover,Omit };
    // An optional source offers only its immutable mesh/selected-domain/BC-kind
    // preparation. Current materials, boundary values, charge and clocks are
    // always revalidated and recomputed; no old Field/state owner is retained.
    // Actual accepted Initial/End electrostatics, including same-End AMR.
    // No artificial interval: a=epsilon, dt=0, charge remains the SAME owner.
    static ElectricVofField prepareElectrostaticSnapshot(const QuadtreeMesh&,
        const Interface::GeometricVofPublication& acceptedVof,const ElectricVofCharge& acceptedCharge,
        const ElectricVofMaterials&,const std::vector<ElectricVofBoundary>& boundaryByNativeFace,
        const ElectricVofControls& = {},CellGradientRecovery = CellGradientRecovery::Recover,
        const ElectricVofField* preparedGeometrySource=nullptr);
    bool isElectrostaticSnapshot() const;
    static ElectricVofField prepareOhmic(const QuadtreeMesh&,
        const Interface::GeometricVofPublication& frozenVof,const ElectricVofCharge& source,
        LAMP::TimeLayer candidateLayer,const ElectricVofMaterials&,
        const std::vector<ElectricVofBoundary>& boundaryByNativeFace,const ElectricVofControls& = {},CellGradientRecovery = CellGradientRecovery::Recover,
        const ElectricVofField* preparedGeometrySource=nullptr);
    // Actual transport receipt supplies Qadv at its Midpoint/End, while
    // sourceCharge remains the original accepted Begin and sets Ohmic dt.
    static ElectricVofField prepareOhmic(const QuadtreeMesh&,
        const Interface::GeometricVofPublication& frozenVof,const ElectricVofAdvectedCharge&,
        const ElectricVofMaterials&,const std::vector<ElectricVofBoundary>& boundaryByNativeFace,
        const ElectricVofControls& = {},CellGradientRecovery = CellGradientRecovery::Recover,
        const ElectricVofField* preparedGeometrySource=nullptr);
    const Interface::GeometricVofPublication& geometry() const;
    const ElectricVofCharge& sourceCharge() const;
    const ElectricVofAdvectedCharge* advection() const;
    const std::vector<double>& implicitInputCharge() const;
    const ElectricVofCharge& chargeCandidate() const;
    LAMP::TimeLayer geometryTime() const;
    LAMP::TimeLayer chargeSourceTime() const;
    LAMP::TimeLayer candidateTime() const;
    double effectiveDt() const;
    // Dense native outputs: inactive phi/material/gradient are zero sentinels.
    // faces() retains dense native rows; both-excluded rows have no measure,
    // material or flux and accept only Interior/zero boundary input.
    // Nearest binary64 projection of the same owned physical potential sum.
    const std::vector<double>& potential() const;
    // Exact sum of the actual configured solve value and its original two fixed-operator
    // defect deltas. Native indexed, inactive entries {0,0,0}; same owner/full clock.
    // Face H and every LS gradient consume these same three terms before RN.
    const std::vector<std::array<double,3>>& potentialTerms() const;
    // Native-cell indexing. The original coefficient-independent LS applied
    // to THIS solved three-term potential; E=-gradient. Not face H/A or a new solve.
    // Omit skips this output; this getter rejects an unrequested recovery.
    const std::vector<Vector2D>& cellPotentialGradients() const;
    const std::vector<double>& epsilonByCell() const;
    const std::vector<double>& conductivityByCell() const;
    const std::vector<ElectricVofFace>& faces() const;
    const ElectricVofOhmicTransferReceipt& ohmicTransferReceipt() const;
    const ElectricVofReport& report() const;
    void validateFor(const QuadtreeMesh&,const Interface::GeometricVofPublication&,
        const ElectricVofCharge& source,LAMP::TimeLayer candidateLayer) const;
private:
    // The immutable stress owner may recognize its exact retained Field Data
    // only inside one callback-free validation; no public validation bypass.
    friend class ElectricVofMaxwellStress;
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofField(std::shared_ptr<const Data>);
    static ElectricVofField prepareImpl(const QuadtreeMesh&,
        const Interface::GeometricVofPublication&,const ElectricVofCharge&,LAMP::TimeLayer,
        const std::vector<double>& implicitInput,const ElectricVofAdvectedCharge*,
        const ElectricVofMaterials&,const std::vector<ElectricVofBoundary>&,const ElectricVofControls&,
        bool electrostaticSnapshot=false,Detail::ElectricVofValidationBatch* batch=nullptr,
        CellGradientRecovery = CellGradientRecovery::Recover,const ElectricVofField* preparedGeometrySource=nullptr);
};
}
