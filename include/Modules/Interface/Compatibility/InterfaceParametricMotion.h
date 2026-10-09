#pragma once
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include "Modules/Interface/Compatibility/InterfaceParametricNeck.h"
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"

class NativeLeafDomain;
class QuadtreeMesh;
// The parametric representation uses the same InterfaceState handle and clock,
// with an explicit geometry branch. It owns no legacy scalar snapshot or nodal
// level set. Native mesh/domain ownership is bound at initialization.
class InterfaceParametricMotion {
public:
    // The same native face charts and positive/signed/product integrals as
    // interface motion, with the whole domain occupied by one background
    // material. Returns native-row order; no InterfaceState or curve is made.
    static std::vector<InterfaceParametricTransport::FaceVolumes::Face> backgroundFaceVolumes(
        const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceNormalFluxView&,double duration);
    static InterfaceState initialize(const QuadtreeMesh&,const NativeLeafDomain&,
        InterfaceParametricProjection::Curve initial);
    static InterfaceState initialize(const QuadtreeMesh&,const NativeLeafDomain&,
        std::vector<InterfaceParametricProjection::Curve> initial);
    // Explicit zero-time model operation. Keeps the accepted physical clock,
    // binds parent/children provenance, and publishes no invented face flux.
    // Coupled Hydro must also transfer its material momentum before accepting.
    static InterfaceState prepareNeckSplit(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,const InterfaceParametricProjection::NeckSplitOptions&);
    struct NeckAssessment {
        InterfaceParametricProjection::NeckReadiness readiness;
        InterfaceParametricProjection::NeckSplitReport report;
        std::optional<InterfaceState> state;
    };
    static NeckAssessment assessNeckSplit(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,const InterfaceParametricProjection::NeckSplitOptions&);
    struct MeshTransfer {
        InterfaceState state;
        InterfaceParametricProjection::MeshTransferReport report;
    };
    // Fresh target mesh/domain, same physical domain and clock. Old face
    // fluxes are not reinterpreted on new faces; the next physical step
    // computes new ones. Preparation leaves both accepted owners untouched.
    static MeshTransfer prepareMeshTransfer(const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
        const InterfaceState&,const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,
        const InterfaceParametricProjection::Options&);
    struct Result {
        InterfaceState state;InterfaceParametricTransport::SourceReport transport;
        InterfaceParametricTransport::FluxReport flux;double inventoryBalance=0;
        // The actual transient trajectory before constrained reconstruction;
        // available for source/accepted-geometry error attribution only.
        const std::vector<InterfaceParametricProjection::Source>& transportSources() const noexcept {return sources_;}
    private:
        // Transient data for this prepared motion only, never a second
        // persistent interface or a retained velocity callback.
        std::vector<InterfaceParametricProjection::Source> sources_;
        std::shared_ptr<const InterfaceParametricProjection::Boundary> projected_;
        InterfaceParametricProjection::Options projection_;
        friend class InterfaceParametricMotion;
    };
    // Enrich the material approximation of the SAME prepared trajectory and
    // independent inventory. Existing candidate knots and every original
    // error budget are retained. No Flow/flux callback is evaluated again.
    static void refineMaterialSpace(Result&,const std::vector<double>& minimumKnots);
    // Finalize an explicitly absorbing outlet after coupled acceptance checks,
    // before material/interface owners are committed together by the caller.
    // Retain physical clock, face flux and previous-geometry lineage; only
    // fully absent components move into the immutable numerical remainder log.
    // This is not enabled by ordinary motion or by a small-volume threshold.
    static void retireExteriorComponents(Result&);
    // Physical time and characteristic parameter duration are separate. For a
    // separable velocity a(t)U(x), duration is integral a(t)dt and both supplied
    // callbacks describe U. The caller must split at a change of sign of a(t).
    static Result prepare(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState& old,
        const InterfaceStepTime&,double parameterDuration,const InterfaceParametricTransport::Flow&,
        const InterfaceNormalFluxView& spatialFlux,const InterfaceParametricProjection::Options&);
    // Inventory slots mapped to actual mesh cell indices, never labels.
    // Uniform layouts retain x-major ordering; mixed leaves bind explicit slots.
    static const std::vector<int>& cellIndices(const InterfaceState&);
    static const InterfaceParametricTransport::FaceVolumes::Face& faceVolumes(const InterfaceState&,int nativeRow);
    static void validateDomain(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&);
    // Both endpoint fluxes and Flow use actual elapsed physical time over
    // this step; velocity may change direction within it. No separable-time
    // substitution or frozen average is performed by this entry point.
    static Result prepareTimeLinear(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState& old,
        const InterfaceStepTime&,const InterfaceParametricTransport::Flow&,
        const InterfaceNormalFluxView& beginFlux,const InterfaceNormalFluxView& endFlux,
        const InterfaceParametricProjection::Options&);
private:
    static InterfaceState bind(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::function<InterfaceParametricProjection::Candidate(const InterfaceParametricProjection::Grid&)>&,
        InterfaceStepProgress);
    static Result prepareInterval(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&,
        const InterfaceStepTime&,double,const InterfaceParametricTransport::Flow&,
        const InterfaceNormalFluxView&,const InterfaceNormalFluxView&,const InterfaceParametricProjection::Options&);
};
