#pragma once
#include "Core/Base/Vector2D.h"
#include <cstddef>
#include <vector>
#include <memory>
#include <optional>
#include <cstdint>

class QuadtreeMesh;
class NativeLeafDomain;
class InterfaceState;
namespace InterfaceParametricProjection {struct NeckSite;}

// Read-only refinement policy for conservative Interface states. Targets are
// native leaf indexed on THIS source mesh; remap them through the existing
// native adaptation machinery, then transfer the complete conservative state.
// This policy neither changes C/geometry nor accepts an AMR or physical step.
namespace InterfaceResolutionPlan
{
struct Options
{
    int maximumLevel=0;
    // Maximum absolute principal curvature times local h; never use their
    // signed sum, which can cancel in a neck. RZ includes the azimuthal part.
    double maxCurvatureCells=.125;
    double minimumOpposedDistanceCells=8;
    double maximumGeometryFractionDifference=.02;
    double distanceBandCells=4;
    double maximumTravel=0; // Caller-certified physical distance until next review.
    double reserveBandCells=0; // Extra target-level cells only when the required band changes mesh.
    // Explicit physical-consumer opt-in, not dictionary-wide AMR permission.
    // Only curvature requests on the contiguous terminal RZ radial neck may
    // exceed the finest level; gap/quality/inventory defects still reject.
    std::optional<double> terminalAxialSection;
    // The coupled consumer explicitly enables the canonical parametric neck
    // model. Its terminal radius is finest max(hx,hy)/maxCurvatureCells.
    // Only bounded axial sublevel branches and event-born poles may use it;
    // all their support is refined to maximumLevel. No scalar-path exception.
    double parametricNeckRadiusCells=0;
};
struct TerminalRange {
    std::uint64_t component=0,birthParent=0;
    double begin=0,end=0,lowerZ=0,upperZ=0,time=0;
    bool pole=false;
};
struct Result
{
    struct UnresolvedFeature {
        enum class Kind {Curvature,OpposedGap,Inventory};
        Kind kind=Kind::Curvature;
        int cell=-1;
        std::uint64_t component=0,otherComponent=0;
        double begin=0,end=0,observed=0,allowed=0,time=0;
        Vector2D point,otherPoint;
    };
    // Actual witnesses for requests beyond the configured cap. Parametric
    // curvature includes polynomial stationary points and leaf intersections;
    // opposed-gap observations remain sampled normal-ray intersections.
    // Read-only diagnostics, never a new exception to the resolution policy.
    std::vector<UnresolvedFeature> unresolved;
    std::vector<int> targetLevels;
    // Parametric geometry requests local levels; interfaceLevel is their
    // maximum (including retained leaves), not a uniform interface target.
    bool localInterfaceLevels=false;
    int interfaceLevel=0;
    // At saturation this stops at maximumLevel+1: a lower bound on the
    // required level, not a claim that exactly one further level suffices.
    int requiredInterfaceLevel=0;
    std::size_t interfaceCells=0,refinedCells=0,curvatureRequests=0,gapRequests=0;
    std::size_t qualityRequests=0,inventoryGeometryRequests=0;
    std::size_t azimuthalRequests=0;
    std::size_t terminalCurvatureRequests=0;
    std::size_t terminalGapRequests=0;
    double maxTerminalCurvatureCells=0;
    double terminalRadius=0,finestNeckScale=0,terminalGeometricVolume=0;
    std::vector<TerminalRange> terminalRanges;
    double maxObservedCurvatureCells=0;
    double maxObservedAzimuthalCells=0;
    // Raw scalar geometry remains diagnostic where the physical
    // height/local-fit provider is valid, not a second curvature authority.
    double maxRawCurvatureCells=0;
    std::size_t rawOnlyCurvatureCells=0;
    std::size_t requiredLeafCount=0,plannedLeafCount=0;
    bool reserveApplied=false,reserveBudgetLimited=false;
    // Numerical extrema are not a curvature-error certificate; sampled
    // normal rays do not bound the continuous minimum gap.
    // Saturation is an unmet request, not success at a silently clamped level.
    bool saturated=false;
    bool changesMesh() const noexcept {return refinedCells!=0;}
};
// Native XY/RZ fluid domains, including declared fixed contacts. Parametric
// states read curvature/gap/inventory directly from the unique Hermite curve;
// the common band/target/budget machinery and Core allocator are unchanged.
// Parametric curvature observes numerical polynomial extrema on each leaf;
// gaps remain sampled normal rays. Neither is a certified interval bound.
// Parametric geometry uses local physical scales and mixed-leaf reconstruction;
// the older scalar reconstruction retains its common interface level.
// Never coarsens: information-loss budgets and hysteresis require a separate
// coarsening policy. The source must already have a usable donor support band.
// A positive cellBudget limits optional reserve only; if the mandatory plan
// does not fit, retain it so refine() reports the original resource failure.
// Zero means no reserve budget here; refine() always enforces its hard budget.
Result evaluate(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&,const Options&,int cellBudget=0);
// Adaptive topology events must belong to the reviewed terminal branch and
// use two finest native axis cells, never an early coarse-grid cutoff.
bool admitsNeck(const QuadtreeMesh&,const InterfaceState&,const Result&,int maximumLevel,
    const InterfaceParametricProjection::NeckSite&);
// Union the support of both physical histories. Parametric target levels are
// combined pointwise; only the older scalar representation needs a common level.
// Leaf counts/refinedCells describe the union; request counters sum visits in
// the two layers. Optional reserve must fit the combined budget.
Result evaluate(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState& current,
    const InterfaceState& previous,const Options&,int cellBudget=0);
// Construct an unaccepted native refinement candidate from ALL source leaves,
// including solid support at reentrant corners. The source is never mutated.
// Reuses Core's 2:1 allocation. Cell budget/capability failures throw length_error.
// The caller must rebind the domain/boundaries and transfer physical state.
std::unique_ptr<QuadtreeMesh> refine(const QuadtreeMesh&,const Result&,int cellBudget);

struct CoarseningPlan {
    std::vector<int> targetLevels;
    std::size_t families=0;
    std::size_t plannedLeafCount=0;
};
// One generation of complete sibling families, with parent-size geometric
// support and stricter curvature/gap thresholds (0 < hysteresis < 1).
// Physical curves and inventory are not modified. Protected points retain
// boundary partitions. The consumer must still assess ALL physical fields;
// this is geometric eligibility, not permission to accept a transfer.
CoarseningPlan planCoarsening(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&,
    const Options&,double hysteresis,const std::vector<Vector2D>& protectedPoints={});
std::unique_ptr<QuadtreeMesh> coarsen(const QuadtreeMesh&,const CoarseningPlan&);
}
