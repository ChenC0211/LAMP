#pragma once
#include "Core/Geometry/AxisymmetricNeckCap.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include <array>
#include <memory>

class NativeLeafDomain;
class QuadtreeMesh;

// Zero-time closure of the declared four-cell interior neck. The ordinary
// material trajectory ends on sourceState; only the two local half-necks are
// remapped in normalized physical area to their respective child caps.
// Every remote panel is an exact native identity. This is a numerical event
// model, not an inferred resolved trajectory through the singular neck.
struct InterfaceInteriorMaterialTransfer
{
    struct Target {std::uint64_t surfaceId=0,liquidId=0;std::shared_ptr<const LevelSetInterfacePanels> panels;};
    struct Weight {std::size_t targetSurface=0,targetPanel=0;long double share=0;};
    std::shared_ptr<const InterfaceState> sourceState;
    std::vector<Target> targets;
    std::vector<std::vector<Weight>> weights;
    InterfaceStepTime time;
    std::array<int,4> cells{};
    double section=0;
};

// A geometry proposal bound to an actual conservative source state.
// Currently supports an upper axial reflection section, two equal native
// cells at r=0, one radial graph, and certified geometric gas outside those
// radial cells in the event band. Its sub-resolution inventory tails retain
// their original C and remain in geometry/inventory diagnostics; resolved
// liquid there is rejected. No minority inventory is erased.
// The complete source scalar trace at the left seam is retained; a monotone
// one-parameter fit preserves its first-cell inventory. Changed non-event
// cells must remain certified gas. The read-only reconstructed state retains
// C and clock, but records a numerical section-detachment proposal at that
// clock, not an estimated physical pinch time. A move port returns a module
// candidate for an uncharged consumer's transaction, never an accepted step.
// Its optional material closure is local to the two event cells. It follows
// the ordinary physical material step and introduces no current or time step.
class InterfaceNeckCap
{
public:
    static InterfaceNeckCap prepare(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,double upperSection,double maximumNeckRadiusCells);
    // Consumes an unpublished uncharged candidate. At the declared finest
    // level, explicitly tag the contiguous unresolved neck/cap up to the first
    // resolved radial section, assessed without any inherited model tags.
    // Its ordinary endpoint/midpoint knots must satisfy the same availability
    // and normal-orientation requirements as their force interpolation.
    // Repeating selection on the same geometry/inventory preserves its scope.
    // After witnessed section detachment the implicit
    // cap may turn; this is not a new radial cap reconstruction. C, geometry,
    // IDs and time are not modified. Force-projection ambiguity remains strict.
    static std::unique_ptr<InterfaceState> prepareTerminalCurvature(const QuadtreeMesh&,
        const NativeLeafDomain&,std::unique_ptr<InterfaceState>,double upperSection,int maximumLevel,
        bool includeInteriorNecks=false);
    // Called inside conservative synchronization BEFORE ordinary identity
    // matching. A shrinking interior section must be witnessed by the old
    // physical scalar and its transported predictor at the declared cutoff.
    // The current synchronized patch supplies the remote scalar seams; the
    // independently advanced inventory supplies all four separate volumes.
    // On success, initializes destination with an explicit 1:2 split and
    // parent/child IDs. The caller attaches that SAME inventory/diagnostics
    // and the NS transaction must still succeed. When material is requested,
    // retain the connected synchronized endpoint and local zero-time closure;
    // the ordinary physical material leg must still be constructed.
    // False leaves patch/destination untouched; inadmissible candidates throw.
    static bool prepareInterior(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState& previous,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor,
        const InterfaceConservativePhase&,const InterfaceStepTime&,
        Quantum::Geometry::ActiveHermitePatch& synchronizedPatch,
        double maximumNeckRadiusCells,int maximumLevel,InterfaceState& destination,
        std::shared_ptr<const InterfaceInteriorMaterialTransfer>* material=nullptr);
    // Prospective finest-grid neck corridor for AMR's explicit terminal model
    // accounting. Does not change curvature or approve gap/quality failures.
    static std::vector<int> terminalResolutionCells(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,double upperSection,int maximumLevel,double maxCurvatureCells);
    const auto& sourceGeometry() const noexcept {return source_;}
    const auto& sourceInventory() const noexcept {return inventory_;}
    const auto& localGeometry() const noexcept {return local_;}
    const auto& localPatch() const noexcept {return patch_;}
    const auto& cells() const noexcept {return cells_;}
    const auto& model() const noexcept {return model_;}
    const auto& nativeVolumes() const noexcept {return nativeVolumes_;}
    const InterfaceState& reconstructedState() const noexcept {return *complete_;}
    struct MaterialTransfer {
        struct Weight {std::size_t targetPanel=0;long double share=0;};
        std::shared_ptr<const LevelSetInterfacePanels> source,target;
        std::vector<std::vector<Weight>> weights;
        std::array<int,2> cells{{-1,-1}};
        std::uint64_t sourceSurfaceId=0,targetSurfaceId=0,liquidId=0;
        InterfaceStepProgress clock;
        double beginZ=0,upperSection=0;
        long double sourceLocalArea=0,targetLocalArea=0;
    };
    // Zero-time numerical closure: remote panels retain their exact material
    // identity. Within the local cap, preserve order in normalized physical
    // surface area and distribute each old panel integral by positive overlap.
    // This is an explicit unresolved-event model, not a resolved trajectory.
    MaterialTransfer materialTransfer() const;
    std::unique_ptr<InterfaceState> takeUnchargedCandidate(const InterfaceStepTime&) &&;
    // The caller must consume this matching material closure before jointly
    // accepting its material fields and this geometry. Neither port commits.
    std::unique_ptr<InterfaceState> takeMaterialCandidate(const InterfaceStepTime&,const MaterialTransfer&) &&;
    double scalarScale() const noexcept {return scalarScale_;}
    double beginZ() const noexcept {return beginZ_;}
    double sourceNeckRadius() const noexcept {return neckRadius_;}
    double cutoffRadiusCells() const noexcept {return cutoffRadiusCells_;}
    // Actual retained source/candidate scalar mismatch at the left seam.
    // Even zero is not permission to publish a physical event.
    double seamScalarDefect() const noexcept {return seamScalarDefect_;}
private:
    std::unique_ptr<InterfaceState> takeCandidate(const InterfaceStepTime&);
    InterfaceNeckCap(const InterfaceState&,Quantum::Geometry::AxisymmetricNeckCap);
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_,local_;
    std::shared_ptr<const LevelSetInterfacePanels> sourcePanels_;
    std::uint64_t sourceSurfaceId_=0,liquidId_=0;
    InterfaceConservativePhase inventory_;
    Quantum::Geometry::AxisymmetricNeckCap model_;
    Quantum::Geometry::ActiveHermitePatch patch_;
    std::unique_ptr<InterfaceState> complete_;
    std::array<int,2> cells_{{-1,-1}};
    std::array<long double,2> nativeVolumes_{};
    double beginZ_=0,neckRadius_=0,cutoffRadiusCells_=0,seamScalarDefect_=0,scalarScale_=1;
};
