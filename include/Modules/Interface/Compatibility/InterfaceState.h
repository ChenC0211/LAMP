#pragma once
#include "Modules/Interface/Compatibility/InterfaceCellMeasureState.h"
#include <cstdint>
#include <memory>
#include <optional>
#include <variant>
#include <vector>
#include "Core/Geometry/ActiveHermitePatch.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"

class QuadtreeMesh;
class InterfaceStateAccess;
class InterfaceTopology;
class InterfaceTransportDomain;
class InterfaceNormalFluxView;
class InterfaceConservativePhase;
class NativeLeafDomain;
struct InterfaceConservativeGeometryReport;
struct InterfacePhaseState;
class InterfaceParametricMotion;
namespace InterfaceParametricProjection {class Candidate;struct NeckEvent;}
namespace InterfaceParametricTransport {struct FaceVolumes;}
// The App supplies its sole clock and the original numerical dt separately:
// end may be step*dt, avoiding a second accumulated clock or a changed dt.
struct InterfaceStepTime
{
    int index = 1;
    double begin = 0.0, end = 0.0, dt = 0.0;
};
// Completed local Interface state, not coupled-step acceptance. Initialization
// starts at zero; zero-time mesh transfer preserves this record exactly.
struct InterfaceStepProgress
{
    int completedSteps = 0;
    double endTime = 0.0;
};
using InterfaceScalarPatch = std::variant<Quantum::Geometry::ActiveHermitePatch,
    Quantum::Geometry::ActiveQuarticHermitePatch>;
struct InterfaceInitialContactJet
{
    Vector2D inwardTangent, normal;
    double dx = 0.0, dy = 0.0, dxx = 0.0, dxy = 0.0, dyy = 0.0;
    double meridionalCurvature = 0.0, axisymmetricCurvature = 0.0;
    int admittedSamples = 0;
};
// Material provenance after the declared interior split. These inventories
// partition the independently transported liquid, including geometric gas
// tails; they are not reconstructed from the current connected regions.
struct InterfaceComponentInventory
{
    std::uint64_t liquidId=0;
    std::vector<long double> liquidVolumes;
};

// Module-owned scalar, constraints, display data and reconstruction workspace.
// App owns this move-only handle and its case IDs, never writable physical
// arrays. Prepared distance data are not yet a usable Hermite interface.
// Producers fill a fresh destination; failed preparation is unusable.
// Success is an Interface candidate, not coupled-step/AMR acceptance.
class InterfaceState
{
public:
    enum class Stage { Empty, Prepared, Ready, Failed };
    explicit InterfaceState(const QuadtreeMesh& mesh);
    ~InterfaceState();
    InterfaceState(const InterfaceState&) = delete;
    InterfaceState& operator=(const InterfaceState&) = delete;
    InterfaceState(InterfaceState&&) noexcept;
    InterfaceState& operator=(InterfaceState&&) noexcept;
    Stage stage() const noexcept;
    InterfaceStepProgress progress() const;
    bool hasParametricGeometry() const noexcept;
    const InterfaceParametricProjection::Candidate& parametricGeometry() const;
    const InterfaceParametricTransport::FaceVolumes& parametricFluxes() const;
    const std::vector<InterfaceParametricProjection::NeckEvent>& parametricNeckEvents() const;
    bool parametricModelCurvature(std::uint64_t component,double materialParameter) const;
    // Birth event of this modeled material, followed through later splits.
    // A missing value means original material; it is not a resolution waiver.
    std::optional<std::uint64_t> parametricModelEvent(std::uint64_t component,double materialParameter) const;
    bool parametricFollows(const InterfaceState& previous,const InterfaceStepTime&) const noexcept;
    void validateMesh(const QuadtreeMesh& mesh) const;
    // Bind initial/AMR phase measures and fluxes to this ready geometry once.
    // Failure leaves geometry and the absent phase state unchanged. The flux
    // provider supplies the actual velocity layer; it is not retained.
    void initializePhase(const QuadtreeMesh& mesh, const InterfaceTransportDomain& domain,
        const InterfaceNormalFluxView& flux);
    bool hasPhase() const noexcept;
    const InterfacePhaseState& phase() const;
    // Opt-in initial conservative inventory, seeded from this state's measured
    // geometry only at time zero: certified pure cells use exact 0/capacity,
    // mixed cells use the raw integral (strict bounds, no clipping). The raw
    // geometric measure remains in phase() for comparison. Synchronization must never
    // overwrite it from phase(). Legacy motion/AMR ports reject this state
    // until their conservative consumers are implemented.
    void initializeConservativePhase(const QuadtreeMesh& mesh, const NativeLeafDomain& domain);
    bool hasConservativePhase() const noexcept;
    const InterfaceConservativePhase& conservativePhase() const;
    const std::vector<InterfaceComponentInventory>& componentInventories() const;
    bool hasConservativeGeometry() const noexcept;
    const InterfaceConservativeGeometryReport& conservativeGeometry() const;
    // Native cells using the declared terminal neck/cap curvature provider.
    // Their curvature is model data, not a resolved height-function estimate.
    // This tag alone does not bypass conservative geometry synchronization;
    // transported-cap geometry protection additionally needs a recorded event.
    // Only module producers can install this scope; no global raw fallback.
    const std::vector<int>& eventCurvatureCells() const;
    // Synchronous read-only borrows expire on handle replacement/destruction.
    // Retaining geometry() by value owns the immutable scalar.
    int vertexCount() const;
    double displayValue(int vertex) const;
    const std::vector<std::uint8_t>& activeVertices() const;
    const std::vector<int>& pinnedVertices() const;
    const std::vector<Vector2D>& pinnedCoordinates() const;
    const std::vector<InterfaceInitialContactJet>& initialContactJets() const;
    const InterfaceScalarPatch& scalarPatch() const;
    const std::vector<Quantum::Geometry::ActiveHermitePatch::NodalDof>& nodalData() const;
    const std::vector<Quantum::Geometry::ActiveHermitePatch::CornerSector>& cornerSectors() const;
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& geometry() const;
    const LevelSetInterfacePanels* panels() const;
    bool hasTopology() const noexcept;
    const InterfaceTopology& topology() const;
private:
    InterfaceState(); // Only the parametric producer creates this representation.
    struct Data;
    std::unique_ptr<Data> data_;
    struct ParametricData;
    std::unique_ptr<ParametricData> parametric_;
    friend class InterfaceStateAccess;
    friend class InterfaceParametricMotion;
};
