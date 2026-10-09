#pragma once
#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include <optional>

namespace InterfaceParametricProjection {
// Explicit unresolved RZ neck closure, not a smooth characteristic map or a
// resolved singularity. Only two adjacent axis cells may change geometry.
// The original independent cell inventory is partitioned, never reseeded.
struct NeckSplitOptions {
    std::uint64_t parent=0;
    double section=0;
    double maximumRadiusCells=0.1;
    double inventoryFractionTolerance=1e-10;
};
struct NeckSplitReport {
    std::uint64_t parent=0;
    std::array<std::uint64_t,2> children{};
    std::array<int,2> cells{{-1,-1}};
    double section=0,lowerSection=0,upperSection=0;
    double neckRadius=0,neckPosition=0,sectionRadius=0,localScale=0;
    double lowerPole=0,upperPole=0;
    double parentVolume=0,modelRegionVolume=0;
    std::array<double,2> childVolumes{};
    double maximumInventoryFractionError=0,maximumCellBalance=0;
    double seamPosition=0,seamTangent=0,seamCurvature=0;
    std::array<bool,2> capAdmissible{};
    std::array<double,2> seamRadialTangent{},seamSignedCurvature{},capLength{},capDiscriminant{};
    std::array<std::array<double,2>,2> retainedParentIntervals{},capIntervals{};
};
struct NeckEvent {NeckSplitReport report;int completedSteps=0;double time=0;};
struct NeckSplit {Candidate candidate;NeckSplitReport report;};
enum class NeckReadiness {ScalePending,CapPending,Ready};
struct NeckAssessment {
    NeckReadiness readiness=NeckReadiness::ScalePending;
    NeckSplitReport report;
    std::optional<NeckSplit> split;
};
// A finite regular neck can remain outside the event feasible set. Only a
// Ready assessment owns a split. Invalid geometry/numerical failures still
// throw; they are never converted into a successful or deferred event.
NeckAssessment assessNeckSplit(const Candidate&,const NeckSplitOptions&);
NeckSplit prepareNeckSplit(const Candidate&,const NeckSplitOptions&);
struct NeckSite {NeckSplitOptions options;double radius=0,position=0,localScale=0;};
std::vector<NeckSite> inspectNecks(const Candidate&,double maximumRadiusCells);
// Inspect actual radial minima. Select the nearest shared axis-cell plane;
// Radius screening alone is not event acceptance; assessNeckSplit also tests
// the selected section and the conservative cap feasible set.
std::vector<NeckSplitOptions> detectNecks(const Candidate&,double maximumRadiusCells);
}
