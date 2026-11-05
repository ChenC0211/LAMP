#pragma once

#include "Modules/Hydrodynamics/HydroPhaseGeometry.h"
#include "Core/Geometry/LevelSetMappedInterface.h"
#include <array>
#include <cstdint>
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace LAMP::Mesh { class NativeLeafOverlap; }
namespace LAMP::Hydrodynamics
{
namespace Detail { struct HydroRegularOperatorTemplateBank; class HydroRegularAdvectorGeometry; struct HydroSparseValidationPlanCache; struct HydroRegularRemapPreparation; struct HydroFrozenMatrixPlanCache; struct HydroStationaryNoSlipStrainCache; }
// Immutable native incidence only; owned by one genuine geometry, never a state.
namespace Detail
{
struct HydroOpenConvectionSupport
{
    struct Pair { int first=-1,second=-1; };
    struct Incidence { std::uint32_t id=0; int sign=1; };
    std::vector<Pair> pairs;
    std::vector<int> segmentPairSlot,boundaryRowSlot,exteriorVelocityDofs,boundaryRows;
    std::vector<std::uint32_t> pairOffsets,boundaryOffsets,boundaryItems;
    std::vector<Incidence> pairItems;
    std::size_t bytes=0;
};
}
struct HydroMomentumOverlap
{
    int pressure = -1; // Same native pressure layout, never a virtual cell.
    Vector2D lower{0,0}, upper{0,0};
    long double volume = 0, cellFraction = 0;
};
struct HydroMomentumVolume
{
    Vector2D lower{0,0}, upper{0,0};
    long double volume = 0;
    std::vector<HydroMomentumOverlap> overlaps;
    std::array<std::vector<int>,4> sides; // -x/r,+x/r,-y/z,+y/z
};
struct HydroMomentumSegment
{
    int negativeVelocity = -1, positiveVelocity = -1; // -1 outside selection.
    int normalAxis = 0;
    Vector2D begin{0,0}, end{0,0};
    long double area = 0;
    // Exactly one origin: an actual primal face or an active cell interior.
    int primalVelocity = -1, ownerPressure = -1;
};
struct HydroMomentumGeometryReport
{
    long double maxRelativeCellPartition = 0, maxRelativeBoxClosure = 0;
    std::size_t overlaps = 0, sharedSegments = 0, boundarySegments = 0;
};
// Cumulative actual preparation work for this owned, bounded mathematical bank.
// Counts include attempted preparations, not just accepted physical steps.
struct HydroRegularTemplateStatistics
{
    std::size_t entries=0,rowEntries=0,routingEntries=0,retainedBytes=0,entryLimit=0,byteLimit=0;
    std::size_t preparations=0,lookups=0,rowQueries=0,rowHits=0,rowBuilds=0,routingHits=0;
    std::size_t portQueries=0,portHits=0,portBuilds=0,fallbacks=0,actualSvdCalls=0,publishedEntries=0;
};
// One actual target geometry construction; never per-step work or a bank sum.
struct HydroRegularPreparationReuseReport
{
    bool requested=false,enabled=false;
    std::size_t targetDatums=0,reusedRows=0,reusedPorts=0,recomputedRows=0,recomputedPorts=0;
    std::size_t stencilRejectedRows=0,portSpanRejectedRows=0,columnRejectedRows=0,fpRejectedRows=0;
    std::size_t retainedMetadataBytes=0;
};
// Derived rectangles and common sides of the EXISTING native face velocities.
// Actual half-cell intersections; no global additional velocity unknowns,
// second primal divergence, pressure gauge, material state or boundary law.
class HydroMomentumGeometry
{
public:
    static HydroMomentumGeometry build(const QuadtreeMesh&, HydroPhysicalDivergence);
    // Fresh target geometry/layout, sharing only a bounded mathematical library.
    // No source mesh, state, native index or field is retained by that library.
    static HydroMomentumGeometry buildWithRegularTemplates(const QuadtreeMesh&,
        HydroPhysicalDivergence, const HydroMomentumGeometry& sourceTemplateOwner);
    // Authenticated source preparation, consumed once by the real target
    // Regular build. Source/target live native identities and overlap validate.
    static HydroMomentumGeometry buildWithRemappedRegularPreparation(
        const QuadtreeMesh& sourceMesh,const Detail::HydroRegularAdvectorGeometry& source,
        const QuadtreeMesh& targetMesh,HydroPhysicalDivergence targetMetrics,
        const Mesh::NativeLeafOverlap& actualOverlaps);
    HydroRegularPreparationReuseReport regularPreparationReuseReport() const;
    HydroRegularTemplateStatistics regularTemplateStatistics() const;
    const HydroPhysicalDivergence& metrics() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    const std::vector<HydroMomentumVolume>& volumes() const noexcept; // velocity DOFs
    const std::vector<HydroMomentumSegment>& segments() const noexcept; // shared once
    const HydroMomentumGeometryReport& report() const noexcept;
    bool sameIdentity(const HydroMomentumGeometry&) const noexcept;
    void validateMesh(const QuadtreeMesh&) const;
private:
    friend class Detail::HydroRegularAdvectorGeometry;
    friend class HydroSaddleVelocityBlock;
    friend class HydroSinglePhaseSaddle;
    friend class HydroFrozenMidpoint;
    friend class HydroEnergyViscosity;
    friend class HydroKepConvection;
    std::shared_ptr<const Detail::HydroOpenConvectionSupport> openConvectionSupport(const QuadtreeMesh&) const noexcept;
    static std::shared_ptr<Detail::HydroRegularOperatorTemplateBank> makeRegularTemplates();
    std::shared_ptr<Detail::HydroRegularOperatorTemplateBank> regularTemplates() const noexcept;
    static std::shared_ptr<Detail::HydroSparseValidationPlanCache> makeSparseValidationPlans() noexcept;
    std::shared_ptr<Detail::HydroSparseValidationPlanCache> sparseValidationPlans() const noexcept;
    static std::shared_ptr<Detail::HydroFrozenMatrixPlanCache> makeFrozenMatrixAssemblyPlans() noexcept;
    std::shared_ptr<Detail::HydroFrozenMatrixPlanCache> frozenMatrixAssemblyPlans() const noexcept;
    std::shared_ptr<const Detail::HydroRegularRemapPreparation> takeRegularRemapPreparation(bool& requested) const;
    void publishRegularPreparationReuseReport(const HydroRegularPreparationReuseReport&) const;
    static std::shared_ptr<Detail::HydroStationaryNoSlipStrainCache> makeStationaryNoSlipStrainCache() noexcept;
    std::shared_ptr<Detail::HydroStationaryNoSlipStrainCache> stationaryNoSlipStrainCache() const noexcept;
    struct Data;
    explicit HydroMomentumGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
struct HydroMomentumMassReport
{
    std::array<HydroPhaseMeasure,2> totalPhaseVolumes{}; // one partition per component
    std::array<long double,2> totalMass{};
    long double maxRelativeCellPartition = 0, maxRelativeLedgerDefect = 0;
    // LocalPlic only: retain raw per-phase relative error above, plus absolute
    // error/cell capacity and fraction of the declared mixed-scale allowance.
    long double maxLocalPlicCapacityScaledCellPartition = 0, maxLocalPlicPartitionBudgetFraction = 0;
    long double maxRelativeMomentEstimate = 0; // estimate, not a rigorous bound
    std::size_t cutOverlaps = 0;
    // Explicit finite P0 phase partitions only; source stocks are unmodified.
    std::size_t finiteP0PhasePartitions=0,finiteP0PhasePieces=0,finiteP0ZeroPhasePieces=0;
    std::array<std::array<long double,2>,2> maxRelativeP0OwnPhaseCellPartition{}; // component, liquid/gas
    // Passive caller wall segments and nested actual integration worker sums.
    // Each report belongs to this mass's actual source layer; no reused begin work.
    double sourcePreparationSeconds = 0, regionBatchSeconds = 0, massAssemblySeconds = 0;
    double cutPairWorkerSeconds = 0, cutFallbackWorkerSeconds = 0;
    std::size_t cutPairCalls = 0, cutFallbackCalls = 0;
    double cutPairRegionGaugeWorkerSeconds = 0, cutPairCutPreparationWorkerSeconds = 0,
        cutPairCoreWorkerSeconds = 0, cutPairRetryWorkerSeconds = 0, cutPairAssemblyWorkerSeconds = 0;
    std::size_t cutPairCutConstructions = 0, cutPairCoreAttempts = 0, cutPairFaceCutRefinements = 0;
    LevelSetMappedInterface::NativeMomentWorkProfile cutPairCoreWork;
};
// Read-only mass from the SAME retained phase publication. Named factories
// select sharp clipping or cell-average projection; neither updates inventory.
// Geometric VOF consumes its genuine retained geometry owner.
// Typed identities stay distinct from the FirstJet consumer ports.
class HydroMomentumMass
{
public:
    enum class DensityRepresentation { SharpGeometry, CellAverageProjection };
    static HydroMomentumMass build(const QuadtreeMesh&, HydroMomentumGeometry,
        HydroPhaseGeometry, double liquidDensity, double gasDensity);
    // FirstJet: its existing cell-average geometry projection. VOF: P0 material
    // from both genuine retained phase inventories and static overlaps.
    // These represented moments are not sharp phase clipping.
    static HydroMomentumMass buildCellAverageProjection(const QuadtreeMesh&, HydroMomentumGeometry,
        HydroPhaseGeometry, double liquidDensity, double gasDensity);
    DensityRepresentation representation() const noexcept;
    HydroPhaseGeometry::Model model() const noexcept;
    Interface::GeometricVofIdentity vofIdentity() const;
    const HydroMomentumGeometry& geometry() const noexcept;
    const HydroPhaseGeometry& phaseGeometry() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    LAMP::TimeLayer timeLayer() const;
    double liquidDensity() const noexcept;
    double gasDensity() const noexcept;
    // True sharp clipping only; projected mass never publishes fake geometry.
    const std::vector<HydroPhaseMeasure>& phaseVolumes() const;
    // Inner ordering matches geometry().volumes()[u].overlaps; sharp only.
    const std::vector<std::vector<HydroPhaseMeasure>>& overlapPhaseVolumes() const;
    // Material moments in the selected density representation, not inventory.
    const std::vector<HydroPhaseMeasure>& representedPhaseVolumes() const noexcept;
    // Same ordering as geometry().volumes()[u].overlaps, in the SELECTED
    // density representation. Sharp returns its original measured overlaps;
    // P0 returns alpha-weighted static volumes, never fake sharp geometry.
    const std::vector<std::vector<HydroPhaseMeasure>>& representedOverlapPhaseVolumes() const noexcept;
    // Genuine VOF/P0 only. pressure uses this native layout; component and
    // phase are0/1 (x/r,y/z and liquid/gas). Same u/k overlap table as above.
    // A flagged phase uses the shared geometric-CDF finite representation;
    // its canonical endpoint delta must be read from that SAME overlap table.
    bool usesFiniteP0Partition(int pressure,int component,int phase) const;
    const std::vector<double>& masses() const noexcept; // positive finite physical M_u
    const HydroMomentumMassReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, Interface::HermiteGalsIdentity,
        LAMP::TimeLayer) const;
    void validateVofFor(const QuadtreeMesh&, Interface::GeometricVofIdentity,
        Interface::GeometricVofTimeLayer) const;
private:
    static HydroMomentumMass buildVolumeFractionGeometry(const QuadtreeMesh&, HydroMomentumGeometry,
        HydroPhaseGeometry, double liquidDensity, double gasDensity);
    struct Data;
    explicit HydroMomentumMass(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

class HydroAdvectorReconstruction;
class HydroSpatialAdvector;
struct HydroKineticHodgePoint
{
    int velocityDual=-1, pressure=-1;
    Vector2D point{0,0};
    double physicalWeight=0;
    bool liquid=false;
};
struct HydroKineticHodgeReport
{
    std::size_t markedDuals=0, points=0, sourceWeights=0, nonZeros=0,tangentialModes=0;
    long double maxRelativePhaseVolumeDifference=0, maxRelativeMomentEstimate=0;
    double symmetryDefect=0;
};
struct HydroKineticHodgeEnergyChange
{
    long double beginEnergy=0, endEnergy=0, energyChange=0;
    long double midpointInertialWork=0, metricVariation=0;
    long double endpointQuadratureRemainder=0, identityResidual=0;
};
// Local kinetic metric of the supplied SAME prepared reconstruction, integrated
// on the authoritative phase-cut native duals. This assembly does not upgrade
// a C1 reconstruction into a phase-jump field or certify its physical accuracy.
// No velocity/geometry state, mass repair, clock owner or solver fallback.
class HydroKineticHodge
{
public:
    using Matrix=Eigen::SparseMatrix<double,Eigen::RowMajor,Eigen::Index>;
    static HydroKineticHodge build(const QuadtreeMesh&,HydroMomentumMass,
        HydroAdvectorReconstruction);
    const HydroMomentumMass& materialMass() const noexcept;
    const HydroAdvectorReconstruction& reconstruction() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    LAMP::TimeLayer timeLayer() const;
    const Matrix& matrix() const noexcept;
    // Ordering [native normal U, independent tangential q], with ALL cross blocks.
    const Matrix& augmentedMatrix() const noexcept;
    const std::vector<bool>& markedDuals() const noexcept;
    const std::vector<HydroKineticHodgePoint>& points() const noexcept;
    // One complete physical-volume partition, for same-source phase gradients.
    const std::vector<HydroKineticHodgePoint>& strainQuadrature() const noexcept;
    const HydroKineticHodgeReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&,Interface::HermiteGalsIdentity,
        LAMP::TimeLayer) const;
    std::vector<double> apply(const QuadtreeMesh&,const HydroNativeValues&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer) const;
    long double kineticEnergy(const QuadtreeMesh&,const HydroNativeValues&,
        Interface::HermiteGalsIdentity,LAMP::TimeLayer) const;
    std::vector<double> apply(const QuadtreeMesh&,const HydroSpatialAdvector&) const;
    long double kineticEnergy(const QuadtreeMesh&,const HydroSpatialAdvector&) const;
    // Force per physical volume at strainQuadrature(), tested by the SAME
    // full [U,q] velocity basis. No density multiplication or provider owner.
    std::vector<double> bodyForceLoad(const QuadtreeMesh&,const std::vector<Vector2D>&) const;
    // Exact endpoint algebra; the nonzero temporal remainder is retained.
    // This diagnostic does not replace the material F/R/uR energy terms.
    HydroKineticHodgeEnergyChange energyChange(const QuadtreeMesh&,
        const HydroKineticHodge& begin,const HydroKineticHodge& end,
        const HydroNativeValues& beginVelocity,const HydroNativeValues& endVelocity) const;
    HydroKineticHodgeEnergyChange energyChange(const QuadtreeMesh&,
        const HydroKineticHodge& begin,const HydroKineticHodge& end,
        const HydroSpatialAdvector& beginVelocity,const HydroSpatialAdvector& endVelocity) const;
private:
    struct Data;
    explicit HydroKineticHodge(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
