#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"

// Read-only topology owned by InterfaceState. IDs are independent of vector
// order and material phase labels. This inventory is not a material remap.
struct InterfaceSurfaceComponent
{
    std::uint64_t id = 0, liquidId = 0;
    LevelSetInterfacePanels geometry;
    double surfaceMeasure = 0;
    std::vector<std::uint64_t> parents; // Previous-frame surface IDs, not material remapping weights.
};
struct InterfaceLiquidComponent
{
    std::uint64_t id = 0;
    std::vector<int> cells, pureCells;
    bool touchesBoundary = false;
    // Volume INSIDE the active mesh, including actual negative wall/port
    // intervals. This geometric measure is not a transport flux budget.
    double volume = 0;
};
struct InterfaceTopologyEvent
{
    enum class Kind { Split, Merge, Reconnection, SectionDetachment };
    Kind kind;
    double begin = 0, end = 0; // Bracket, not an estimated physical pinch time.
    std::vector<std::uint64_t> parents, children;
    // A section event may retain one liquid ID on both sides. Surface IDs
    // still change, and the lost reflection section is explicit. This is
    // neither a charge/material transfer nor proof of physical pinch timing.
    std::vector<LevelSetInterfacePanels::SymmetryPlane> detachedSections;
    // Source-bound, four-cell RZ interior-neck model. Both opposed caps keep
    // their own two cell inventories and share one C1 positive gas trace.
    // Radii refer to the same native section in the old physical state and
    // its actual transported predictor, not to a synthetic pinch-time fit.
    struct InteriorNeck {
        double section=0,radialWidth=0,axialWidth=0,cutoffCells=0;
        double sourceRadius=0,predictedRadius=0,gasValue=0;
        std::vector<int> cells;
    };
    std::vector<InteriorNeck> interiorNecks;
};
class InterfaceTopology
{
public:
    const auto& surfaces() const { return surfaces_; }
    const auto& liquids() const { return liquids_; }
    const auto& events() const { return events_; }
    const auto& cellPhase() const { return cellPhase_; }
    const auto& symmetryPlanes() const { return symmetryPlanes_; }
private:
    std::vector<InterfaceSurfaceComponent> surfaces_;
    std::vector<InterfaceLiquidComponent> liquids_;
    std::vector<InterfaceTopologyEvent> events_;
    std::vector<int> cellPhase_;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetryPlanes_;
    std::uint64_t nextLiquidId_ = 1, nextSurfaceId_ = 1;
    friend class InterfaceGeometryReconstruction;
    friend class InterfaceMeshTransferTopology;
    friend class InterfaceNeckCap;
};

class InterfaceGeometryReconstruction
{
public:
    // Publish a new immutable scalar and its COMPLETE topology atomically.
    // Generic geometry initialization; source DOFs are copied into the state.
    // Optional outer symmetry sections retain movable endpoints and exact
    // zero-normal-derivative traces. Geometry-only continuation inherits them.
    // Conservative transport still requires its own reflected stencil support.
    static void initialize(const QuadtreeMesh& mesh, InterfaceScalarPatch patch,
        bool axisymmetric, const std::vector<Vector2D>& boundaryContacts,
        InterfaceState& destination,
        const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes={});
    // Manufactured geometry sequence only. Match liquid regions through
    // shared certified pure negative cells on an unchanged mesh. Ambiguous or
    // absent overlap rejects, never invents a birth/death or nearest match.
    // This deliberately does not assert volume/charge conservation, transfer
    // material histories, or authorize a physical topology-changing time step.
    static void prepareManufactured(const QuadtreeMesh& mesh, InterfaceScalarPatch patch,
        bool axisymmetric, const std::vector<Vector2D>& boundaryContacts,
        const InterfaceState& previous, const InterfaceStepTime& time,
        InterfaceState& destination);
private:
    friend class InterfaceComponentMotion;
    friend class InterfaceConservativeGeometry;
    friend class InterfaceConservativeGals;
    static void prepare(const QuadtreeMesh& mesh, InterfaceScalarPatch patch,
        bool axisymmetric, const std::vector<Vector2D>& boundaryContacts,
        const InterfaceState* previous, const InterfaceStepTime* time,
        InterfaceState& destination,
        const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes={},
        // Internal producer may reuse its already constructed immutable
        // scalar for THIS patch. Discovery/lineage/clock checks still run.
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retainedScalar={});
};
