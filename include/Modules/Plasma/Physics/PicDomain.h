#pragma once
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Boundary/LegacyGridMotionEvolution.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <filesystem>
#include <optional>

class PicCheckpointReader;

// Owns one PIC mesh, geometry, derived views and refinement policy. Object
// identities stay stable through restore/AMR; borrowers must not outlive this
// nonmovable owner. Leaf arrays and derived views still need the existing
// explicit transfer/rebuild sequence after topology changes. Geometry input
// borrowed by an attached refiner must outlive its use.
class PicDomain
{
public:
    using Mode = PicEmbeddedGeometry::Mode;
    using Geometry = PicEmbeddedGeometry::Geometry;
    struct Restore
    {
        std::optional<std::filesystem::path> directory;
        // Set only for a checkpoint, not for an initial geometry snapshot.
        std::optional<int> nextStep;
    };
    PicDomain(const DictionaryMesh& meshInput, const DictionaryMesh& refinementInput, Mode mode);
    PicDomain(const PicDomain&) = delete;
    PicDomain& operator=(const PicDomain&) = delete;

    // Run startup has two stages so field-boundary mapping retains its original
    // position after mesh restoration. Each stage is attempted at most once;
    // false/exception stops the session, without rollback or retry.
    bool prepareMesh(const PicCheckpointReader* restart,
        std::optional<std::filesystem::path> initialSnapshot = std::nullopt);
    struct InitialState { bool embeddedEnabled, refinementActive; };
    // Initialize/restore geometry, attach refinement, perform startup AMR,
    // save an optional implicit snapshot, then build the graph and cut metrics.
    // MovingEB is still initialized later by its owner; no live wall exists yet.
    // Input borrowed by the refiner must outlive the domain.
    std::optional<InitialState> initializeForRun(DictionaryGeometry* input,
        bool embeddedEnabled, bool refinementEnabled, bool geometryChangesDuringRun,
        const std::vector<int>& lineToBoundaryId,
        const LegacyGridMotionConfig* scheduledMotion = nullptr);

    // False/exception stops the caller; partially initialized state is not
    // rolled back or reusable. A successful empty geometry is a valid state.
    bool initializeGeometry(DictionaryGeometry* input, bool& embeddedEnabled,
        bool& refinementEnabled, const Restore& restore,
        const LegacyGridMotionConfig* scheduledMotion = nullptr);
    void refineInitial(const DictionaryGeometry* input, MovingEBController* moving,
        bool embeddedEnabled, bool phaseActive);
    void writeGeometryCache(const std::string& prefix) const;

    QuadtreeMesh& mesh() { return mesh_; }
    const QuadtreeMesh& mesh() const { return mesh_; }
    QuadtreeMeshRefiner& refiner() { return refiner_; }
    const QuadtreeMeshRefiner& refiner() const { return refiner_; }
    FaceSegmentGraph& faceGraph() { return faceGraph_; }
    const FaceSegmentGraph& faceGraph() const { return faceGraph_; }
    // Call at the existing graph rebuild stages, after geometry/composite
    // refresh. This does not transfer fields or rebuild cut metrics/solvers.
    void rebuildFaceGraph();

    Geometry* geometry() { return geometry_.get(); }
    const Geometry* geometry() const { return geometry_.get(); }
    PicEmbeddedGeometry& embedded() { return embedded_; }
    const PicEmbeddedGeometry& embedded() const { return embedded_; }
    bool hasScheduledMotion() const { return scheduledMotion_ != nullptr; }
    // The prescribed-motion runtime service takes this independent owner after
    // startup; raw geometry and all derived views remain owned here.
    std::unique_ptr<LegacyGridMotionEvolution> takeScheduledMotion();
private:
    friend class PicSpeciesState;
    friend class ImplicitPicRefinement;
    friend class ExplicitPicRefinement;
    QuadtreeMesh mesh_;
    const DictionaryMesh refinementInput_;
    const Mode mode_;
    bool meshPreparationAttempted_ = false, meshPrepared_ = false;
    Restore startupRestore_;
    std::optional<std::filesystem::path> initialSnapshot_;
    bool initialSnapshotLoaded_ = false;
    bool initializationAttempted_ = false;
    std::unique_ptr<Geometry> geometry_;
    // Destroy borrowers before their original level-set storage.
    PicEmbeddedGeometry embedded_;
    std::unique_ptr<LegacyGridMotionEvolution> scheduledMotion_;
    // The refiner borrows our policy and may borrow our geometry.
    QuadtreeMeshRefiner refiner_;
    FaceSegmentGraph faceGraph_;
};
