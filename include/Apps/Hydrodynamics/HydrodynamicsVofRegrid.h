#pragma once
#include "Core/Base/TimeLayer.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>
namespace LAMP::Hydrodynamics { class HydroTwoPhaseState; }
// App policy, not a physical remap or accepted clock owner.
struct HydrodynamicsVofRegridPolicy
{
    int minLevel=0,maxLevel=1,interfaceBandLayers=2;
    double vorticityTolerance=.5,coarsenRatio=.25;
    // Absent means the global minimum. Floors classify actual PLIC geometry,
    // not the mere presence of arbitrarily small independent phase stocks.
    std::optional<int> liquidMinLevel,gasMinLevel;
    // Positive accepted physical-step interval; one retains the original cadence.
    int interval=1;
};
// Uses the complete global End step, including continuation after restart.
bool hydrodynamicsVofRegridDue(const HydrodynamicsVofRegridPolicy&,const LAMP::TimeLayer&);
// Frozen App planning requirements, not another physical state or remap.
// Owns the minimum-level payload and retains the genuine source handle for
// this synchronous planning transaction. No foreign physics dependency or
// caller-supplied mesh/clock certificate is accepted.
// Requires exact selected native maps at a complete End, with one integer
// per native leaf in [0,allocated maxLevel]. Consumption enforces the active
// policy range; excluded leaves must carry exactly policy.minLevel.
class HydrodynamicsVofRegridRequirements
{
public:
    static HydrodynamicsVofRegridRequirements fromCellMinimumLevels(
        const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        std::vector<int>);
    const std::vector<int>& minimumLevels() const;
    LAMP::TimeLayer sourceClock() const;
    void validateFor(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const HydrodynamicsVofRegridPolicy&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HydrodynamicsVofRegridRequirements(std::shared_ptr<const Data>);
};
struct HydrodynamicsVofRegridReport
{
    LAMP::TimeLayer sourceClock;
    // mixedCells is still the independently positive STOCK count.
    // Geometry seeds use the existing immutable resolved/measured certificate;
    // unresolvedMixedCells retains and counts stocks, without inventing a shape.
    std::size_t sourceLeaves=0,targetLeaves=0,mixedCells=0,pureJumpFaces=0;
    std::size_t geometricMixedCells=0,unresolvedMixedCells=0;
    std::size_t interfaceSeedCells=0,interfaceBandCells=0,flowRefineCells=0,flowCoarsenCells=0;
    // Only requirements that raise the original interface/flow target count.
    // Refine exceeds the actual current level; Hold prevents its coarsening.
    std::size_t additionalRefineCells=0,additionalHoldCells=0;
    double maximumHAbsVorticity=0,indicatorSeconds=0,planningSeconds=0,adaptationSeconds=0,totalSeconds=0;
    int clones=0,sweeps=0;
    bool requested=false,changed=false;
    // This one successful requested decision, including its final settle sweep.
    // Nested elapsed diagnostics only, never a physical acceptance condition.
    struct TreePreparation {
        bool available=false;
        double cloneSeconds=0,graphSeconds=0,domainSeconds=0,overlapSeconds=0,nativeAdaptSeconds=0,unattributedSeconds=0;
        std::size_t cloneCalls=0,graphCalls=0,domainCalls=0,overlapCalls=0,nativeAdaptCalls=0;
    } treePreparation;
    // Stock immediately before and after the one actual clone, before adaptation.
    // Order is vertices, horizontal edges, vertical edges, cells; never cumulative.
    struct CloneStorage {
        bool available=false;
        std::array<std::size_t,4> logicalCapacity{},sourceConstructed{},cloneConstructed{};
    } cloneStorage;
};
// Domain belongs to this owned final mesh. Separate Interface/Hydro transfer
// and coupled acceptance are mandatory before this can replace a physical grid.
struct HydrodynamicsVofGridCandidate
{
    std::unique_ptr<QuadtreeMesh> mesh;
    NativeLeafDomain domain;
};
struct HydrodynamicsVofRegridDecision
{
    std::optional<HydrodynamicsVofGridCandidate> candidate;
    HydrodynamicsVofRegridReport report;
};
// Synchronous App geometry factory; never initializes moving phase stocks.
// If supplied, its source maps must exactly match the genuine publication.
// Each actual candidate is requalified, then NativeLeafOverlap certifies the
// same selected physical coverage. Required for a noncomplete source domain.
using HydrodynamicsVofTargetDomainProvider=std::function<NativeLeafDomain(const QuadtreeMesh&)>;
HydrodynamicsVofRegridDecision prepareHydrodynamicsVofRegrid(
    const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
    const HydrodynamicsVofRegridPolicy&,
    const HydrodynamicsVofRegridRequirements* = nullptr,
    const HydrodynamicsVofTargetDomainProvider& = {});
