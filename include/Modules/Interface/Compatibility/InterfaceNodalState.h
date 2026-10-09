#pragma once

#include <functional>
#include <memory>
#include <vector>
#include "Core/Field/VertField.h"
#include "Core/Geometry/LevelSetBoundaryConditions.h"
#include "Modules/Interface/Compatibility/InterfaceBoundaryGeometry.h"

class QuadtreeMesh;
template<typename Mesh> class ImmersedGeometry;
namespace LevelSetHermiteGeometry { class Snapshot; }

// Existing complete-uniform-domain nodal level set, not the masked native
// Hermite representation used by preparePinnedToAxis. Both belong to Interface;
// neither is converted into a second authoritative scalar to share a name.
struct InterfaceNodalOptions
{
    bool axisymmetric = false;
    bool allowAxisContact = false;
    bool prescribedBoundary = false;
    int reinitializeEvery = 0;
    std::vector<Vector2D> fixedContacts;
};

struct InterfaceNodalVelocityLayers
{
    // Ownership of the SAME old/new vertex fields used for transport. A
    // requested material map retains them without copying/reconstructing them.
    std::unique_ptr<VertField<double>> currentX, currentY, previousX, previousY;
    // Required only for prescribed scalar boundaries; owning callbacks, as in
    // the existing Core trace contract. No Hydro types enter this module.
    LevelSetAdvectionVelocity boundaryVelocity;
};

class InterfaceNodalMaterialMap
{
public:
    ~InterfaceNodalMaterialMap();
    InterfaceNodalMaterialMap(const InterfaceNodalMaterialMap&) = delete;
    InterfaceNodalMaterialMap& operator=(const InterfaceNodalMaterialMap&) = delete;

    struct Motion
    {
        const LevelSetHermiteGeometry::Snapshot& snapshot;
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint;
    };
    // Validate the consumer's retained old nodal scalar against the actual
    // transport history, then borrow the final snapshot and SAME trajectory.
    // No second history copy or consumer-built destination geometry is needed.
    // Validate once before a synchronous batch of material queries. This map,
    // the mesh and Interface state must stay alive and unchanged for the batch.
    // A copied callback checks its lease on every query: later advance/failure
    // or state destruction invalidates it. No per-point full mesh scan.
    Motion motionFor(const QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
        const std::vector<double>& sourceScalar, int step, double sourceTime,
        double dt, bool axisContact) const;
private:
    struct Impl;
    explicit InterfaceNodalMaterialMap(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;
    friend class InterfaceNodalState;
};

class InterfaceNodalState
{
public:
    // Adopt an initialized single liquid level set and its existing histories;
    // refresh/check derived geometry here. Geometry initialization/scatter IO
    // remains at the input boundary. No clone of phi, history or mesh is made.
    // The supplied fixed uniform mesh must outlive this state and its queries.
    InterfaceNodalState(QuadtreeMesh& mesh,
        std::unique_ptr<ImmersedGeometry<QuadtreeMesh>> geometry, int interfaceIndex,
        InterfaceNodalOptions options);
    ~InterfaceNodalState();
    InterfaceNodalState(const InterfaceNodalState&) = delete;
    InterfaceNodalState& operator=(const InterfaceNodalState&) = delete;

    // Read-only borrowed geometry, including scalar/history/derivative arrays.
    // The borrow expires on advance or destruction; it is not a frozen past
    // state. New geometry/boundary queries reject a failed numerical substep.
    const ImmersedGeometry<QuadtreeMesh>& geometry() const;
    const InterfaceBoundaryGeometry::Report& boundaryGeometry() const;
    const char* stage() const noexcept;
    bool failed() const noexcept;

    // Same legacy two-layer transport -> scheduled exact redistance -> derived
    // geometry -> topology checks. App supplies its sole clock; this object
    // only checks continuity. Preflight errors do not consume the step; once
    // numerical work begins an exception is fail-stop, with no rollback.
    // Success is an Interface substep, NOT coupled-step acceptance.
    // Material maps are unsupported for prescribedScalar (same prior scope).
    std::unique_ptr<InterfaceNodalMaterialMap> advance(QuadtreeMesh& mesh,
        int step, double sourceTime, double dt, InterfaceNodalVelocityLayers velocity,
        bool retainMaterialMap = false);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
