#pragma once
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include <memory>
#include <optional>

class QuadtreeMeshRefiner;
class PicPoissonFieldSolve;
class ImplicitPicFieldSolve;
namespace ExplicitPIC { namespace Pusher { struct EmbeddedImpact; } }

// Owns the wall controller, lazy per-step impact lookup and last update result.
// initialize is called once at startup; its dictionary/units must outlive us.
// Other geometry/solver objects are borrowed only during the corresponding call.
class PicMovingGeometry
{
public:
    using Mode = PicEmbeddedGeometry::Mode;
    using Geometry = PicEmbeddedGeometry::Geometry;
    struct Refresh
    {
        QuadtreeMesh& mesh;
        Geometry& geometry;
        PicEmbeddedGeometry& embedded;
        FaceSegmentGraph& graph;
        const std::vector<int>& lineToBoundaryId;
        bool embeddedEnabled;
        QuadtreeMeshRefiner* refiner; // null when geometry refinement is inactive
        PicPoissonFieldSolve& poisson;
        ImplicitPicFieldSolve* polarization; // implicit only
    };
    explicit PicMovingGeometry(Mode mode) : mode_(mode) {}
    PicMovingGeometry(const PicMovingGeometry&) = delete;
    PicMovingGeometry& operator=(const PicMovingGeometry&) = delete;
    bool initialize(QuadtreeMesh& mesh, Geometry* geometry,
        const DictionaryGeometry* dictionary, const MovingEBConfig& config,
        const UnitsContext* units, const NormalizationContext* normalization,
        const std::optional<std::filesystem::path>& restartDirectory);
    const MovingEBController* controller() const { return controller_.get(); }
    // Does not reset updated(): existing diagnostics observe the previous update
    // until advance is called at its original end-of-step position.
    void beginStep();
    void onMeshChanged(const QuadtreeMesh& mesh);
    void ingestImpact(const QuadtreeMesh& mesh, const Geometry& geometry,
        const PicEmbeddedGeometry::Metrics* metrics, const std::vector<int>& lineToBoundaryId,
        const ExplicitPIC::Pusher::EmbeddedImpact& impact, int speciesId,
        double mass, double charge, double defaultWeight, double velocityScaleMps, double timeSI);
    bool advance(int step, double timeSI, double dtSI, QuadtreeMesh& mesh, Geometry* geometry);
    bool updated() const { return updated_; }
    // The implicit driver culls particles between advance and this refresh.
    // Failures stop the run, as before; controller/geometry changes are not rolled back.
    void refreshAfterUpdate(int step, Refresh bindings);
    void writeCache(const std::filesystem::path& directory) const;
private:
    friend class ImplicitPicRefinement;
    friend class ExplicitPicRefinement;
    void verifyCutMetrics(int step, const Refresh& bindings) const;
    const Mode mode_;
    std::unique_ptr<MovingEBController> controller_;
    std::unique_ptr<MovingEBController::ImpactResolveContext> impactContext_;
    bool restartEnabled_ = false;
    bool updated_ = false;
};
