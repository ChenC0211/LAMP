#pragma once
#include "Modules/Plasma/Boundary/LegacyGridMotionEvolution.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include <functional>

class QuadtreeMeshRefiner;
class PicPoissonFieldSolve;
class ImplicitPicFieldSolve;

// Owns prescribed motion and its current observed state. Bound module objects
// must outlive this service; their identities remain stable across AMR. A
// changed/forced update rebuilds geometry before accepting the motion. Failure
// stops the caller: partially rebuilt geometry is not rolled back.
class PicScheduledGeometry
{
public:
    struct Bindings
    {
        QuadtreeMesh& mesh;
        PicEmbeddedGeometry::Geometry* geometry;
        PicEmbeddedGeometry& embedded;
        FaceSegmentGraph& graph;
        const std::vector<int>& lineToBoundaryId;
        const MovingEBController* moving;
        QuadtreeMeshRefiner& refiner;
        PicPoissonFieldSolve& poisson;
        ImplicitPicFieldSolve& polarization;
    };
    struct Options
    {
        bool enabled = false;
        bool embeddedEnabled = false;
        bool refinementActive = false;
        // Synchronous diagnostics only; the observer must outlive this service
        // and must not advance or mutate geometry. Exceptions stop the caller
        // at this event without rolling back prior geometry work.
        std::function<void(int, const std::string&)> trace;
    };
    PicScheduledGeometry(std::unique_ptr<LegacyGridMotionEvolution> motion,
        Bindings bindings, Options options);
    PicScheduledGeometry(const PicScheduledGeometry&) = delete;
    PicScheduledGeometry& operator=(const PicScheduledGeometry&) = delete;
    bool hasSchedule() const { return motion_ != nullptr; }
    // True means geometry was rebuilt (including force with unchanged state).
    // The caller uses that result to issue its original AMR request.
    bool advance(int step, double timeSI, bool force);
    const LegacyGridMotionState& current() const { return current_; }
private:
    void trace(int step, const std::string& message);
    std::unique_ptr<LegacyGridMotionEvolution> motion_;
    const Bindings bindings_;
    const Options options_;
    LegacyGridMotionState current_;
};
