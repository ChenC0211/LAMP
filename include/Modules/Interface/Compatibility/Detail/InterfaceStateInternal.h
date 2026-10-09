#pragma once
// Implementation-only storage access, not exported to Apps.
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Modules/Interface/Compatibility/InterfaceParametricProjection.h"
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include "Modules/Interface/Compatibility/InterfaceParametricNeck.h"
#include <stdexcept>
#include <optional>

namespace ActiveQuarticHermiteReconstruction { class StencilCache; }
struct InterfaceState::ParametricData {
    struct Binding {
        const QuadtreeMesh* mesh=nullptr;
        const NativeLeafDomain* domain=nullptr;
        std::shared_ptr<const NativeMeshIdentity> identity;
        std::vector<int> cellIndices,verticalRows,horizontalRows;
        std::vector<std::pair<int,int>> faceSlots;
    };
    std::shared_ptr<const Binding> binding;
    InterfaceParametricProjection::Candidate candidate;
    InterfaceStepProgress progress;
    InterfaceStepTime stepTime;
    std::shared_ptr<const InterfaceParametricProjection::Boundary> previousGeometry;
    InterfaceParametricTransport::FaceVolumes fluxes;
    std::vector<InterfaceParametricProjection::NeckEvent> neckEvents;
    ParametricData(std::shared_ptr<const Binding> b,InterfaceParametricProjection::Candidate c,InterfaceStepProgress p={})
        :binding(std::move(b)),candidate(std::move(c)),progress(p){}
};
struct InterfaceState::Data
{
    explicit Data(const QuadtreeMesh& mesh);
    Stage stage = Stage::Empty;
    InterfaceStepProgress progress;
    VertField<double> freeLevelSet;
    std::vector<std::uint8_t> activeVertex;
    std::vector<int> pinnedVertexIndices;
    std::vector<Vector2D> pinnedCoordinates;
    std::vector<InterfaceInitialContactJet> initialContactJets;
    InterfaceScalarPatch scalarPatch;
    std::shared_ptr<ActiveQuarticHermiteReconstruction::StencilCache> transportStencils;
    std::shared_ptr<const NativeMeshIdentity> meshIdentity;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    std::unique_ptr<LevelSetInterfacePanels> panelData;
    std::unique_ptr<InterfaceTopology> topology;
    std::optional<InterfacePhaseState> phase;
    std::optional<InterfaceConservativePhase> conservativePhase;
    std::vector<InterfaceComponentInventory> componentInventories;
    std::optional<InterfaceConservativeGeometryReport> conservativeGeometry;
    std::vector<int> eventCurvatureCells;
    const auto& geometry() const { return source; }
    const LevelSetInterfacePanels* panels() const;
    void retainGeometry(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> geometry,
        std::unique_ptr<LevelSetInterfacePanels>&& panels);
};
class InterfaceStateAccess
{
public:
    using Data = InterfaceState::Data;
    static Data& get(InterfaceState& state) {
        if (state.parametric_) throw std::logic_error("scalar Interface consumer cannot use parametric geometry");
        if (!state.data_) throw std::logic_error("moved-from Interface state");
        return *state.data_;
    }
    static const Data& get(const InterfaceState& state) {
        if (state.parametric_) throw std::logic_error("scalar Interface consumer cannot use parametric geometry");
        if (!state.data_) throw std::logic_error("moved-from Interface state");
        return *state.data_;
    }
};
class InterfacePreparation
{
public:
    explicit InterfacePreparation(InterfaceStateAccess::Data& state) : state_(state) {}
    ~InterfacePreparation() { if (!completed_) { state_.stage = InterfaceState::Stage::Failed;
        state_.source.reset(); state_.panelData.reset(); state_.topology.reset(); state_.phase.reset();
        state_.conservativePhase.reset();state_.conservativeGeometry.reset();state_.eventCurvatureCells.clear();
        state_.componentInventories.clear(); } }
    void complete(InterfaceState::Stage stage) { state_.stage = stage; completed_ = true; }
private:
    InterfaceStateAccess::Data& state_;
    bool completed_ = false;
};

// Shared zero-time identity binding for remeshing and inventory-driven scalar
// synchronization. No event or physical ancestry edge may be manufactured.
class InterfaceMeshTransferTopology
{
public:
    static void retain(const InterfaceTopology& source, InterfaceTopology& target,
        const std::vector<std::vector<int>>& sourceCellsByTarget);
};
