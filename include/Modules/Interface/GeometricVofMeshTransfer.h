#pragma once
#include "Modules/Interface/GeometricVofPublication.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include <array>
#include <memory>
#include <vector>

namespace LAMP::Interface
{
struct GeometricVofTransferOverlap
{
    enum class Kind { Unchanged, Refined, Coarsened };
    Kind kind=Kind::Unchanged;
    int sourceCell=-1,targetCell=-1;
    long double capacity=0,liquidMeasured=0,gasMeasured=0;
    long double inventoryContribution=0,closureCorrection=0; // Historical fan; LocalPlic aliases liquid only.
    bool localPlic=false,usedP0=false,geometricMeasureAvailable=false;
    std::array<long double,2> phaseContribution{},phaseCorrection{},phaseRoundoffClosure{};
};
struct GeometricVofTransferReport
{
    std::size_t refinedDonors=0,refinedTargets=0,coarsenedTargets=0,unchangedTargets=0;
    std::size_t sourceComponents=0,targetComponents=0,sourceJunctions=0,targetJunctions=0;
    long double maximumRelativeCoverageResidual=0,maximumRelativeRawInventoryResidual=0;
    long double maximumRelativeClosureCorrection=0,maximumRelativeDonorClosureCorrection=0,maximumRelativeMinorPhaseCorrection=0,maximumRelativeInventoryResidual=0;
    long double maximumRelativeGeometryResidual=0;
    long double sourceInventory=0,targetInventory=0,sourceGasInventory=0,targetGasInventory=0,sourceMeasured=0,targetMeasured=0;
    long double closureCorrection=0,conversionResidual=0;
    // Independent LocalPlic ledger; topology counters above are unavailable.
    bool localPlic=false;
    std::size_t resolvedRefinedDonors=0,unresolvedRefinedDonors=0,p0RefinedPieces=0;
    std::size_t correctedRefinedDonors=0,localPlaneCuts=0,roundoffClosedPhases=0;
    std::array<long double,2> sourcePhaseInventory{},targetPhaseInventory{};
    std::array<long double,2> maximumRelativeDonorPhaseResidual{},relativeGlobalPhaseResidual{};
    std::array<long double,2> maximumRelativePhaseCorrection{},maximumRelativePhaseRoundoffClosure{};
    long double maximumRelativeLocalGeometryResidual=0,maximumRelativeTargetPhaseClosure=0;
};
// Prepares ONE immutable regrid candidate at the SAME accepted End clock.
// Both meshes remain caller-owned and unchanged. Failure cannot replace an
// accepted handle. Cell overlaps are physical-volume coverage, NOT Gamma/Q
// history coverage; changed component graphs do not require a bijection.
class GeometricVofMeshTransfer
{
public:
    static GeometricVofMeshTransfer prepare(const QuadtreeMesh& sourceMesh,
        const GeometricVofPublication& source,const QuadtreeMesh& targetMesh,
        NativeLeafDomain targetDomain,GeometricVofTimeLayer acceptedTime,
        const GeometricVofConfig& = {});
    GeometricVofIdentity sourceIdentity() const noexcept { return sourceIdentity_; }
    const GeometricVofPublication& candidate() const noexcept { return candidate_; }
    const std::vector<GeometricVofTransferOverlap>& overlaps() const noexcept { return overlaps_; }
    const GeometricVofTransferReport& report() const noexcept { return report_; }
    // LocalPlic only; the same actual coverage used by both phase ledgers.
    // Historical fan callers retain overlaps(); requesting this port rejects.
    const Mesh::NativeLeafOverlap& nativeOverlaps() const;
private:
    GeometricVofMeshTransfer(GeometricVofIdentity,GeometricVofPublication,
        std::vector<GeometricVofTransferOverlap>,GeometricVofTransferReport,
        std::shared_ptr<const Mesh::NativeLeafOverlap> = {});
    GeometricVofIdentity sourceIdentity_;
    GeometricVofPublication candidate_;
    std::vector<GeometricVofTransferOverlap> overlaps_;
    GeometricVofTransferReport report_;
    std::shared_ptr<const Mesh::NativeLeafOverlap> nativeOverlaps_;
};
}
