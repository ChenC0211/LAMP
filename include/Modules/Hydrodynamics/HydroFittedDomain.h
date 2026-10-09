#pragma once

#include <Eigen/Sparse>
#include <vector>
#include "Core/Mesh/NativeLeafDomain.h"

class QuadtreeMesh;
class FaceSegmentGraph;

// Fixed face-fitted fluid domain on complete Cartesian leaves, including 2:1
// faces whose native pressure interpolation remains entirely in fluid.
// Original mesh labels/parents are retained. No cut-cell/material is inferred.
// Integrated continuity is D*u. On same-level faces C=-D^T; coarse/fine C
// integrates the shared native virtual-child gradient over its fine dual.
class HydroFittedDomain
{
public:
    enum class BoundaryType { Interior, NoSlipWall, Symmetry, PressureReservoir, PressureOutlet };
    struct Boundary
    {
        BoundaryType type = BoundaryType::Interior;
        int id = -1; // Caller metadata; -1 is permitted, never used to infer a BC.
        double pressure = 0.0;
        double incomingDensity = 0.0;
        double incomingViscosity = 0.0;
    };
    // Hydro extends the shared native geometry with its own unknown map,
    // physical BC and virtual-child momentum/pressure dual. No second face
    // geometry array is retained by the Core domain.
    struct NativeFace : NativeLeafFace
    {
        int velocityRow = -1; // Compact free row; -1 for solid/constrained faces
        Boundary boundary;
        double dualVolume = 0.0;
    };

    HydroFittedDomain() = default; // Empty until assigned the result of build().
    static HydroFittedDomain build(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
        const std::vector<Boundary>& bySegment, bool axisymmetric);

    const std::vector<NativeFace>& faces() const noexcept { return faces_; }
    const NativeLeafDomain& geometry() const noexcept { return geometry_; }
    const std::vector<int>& cellRowByIndex() const noexcept { return geometry_.cellRowByIndex(); }
    const std::vector<int>& cellIndexByRow() const noexcept { return geometry_.cellIndexByRow(); }
    const std::vector<int>& velocityRowByNative() const noexcept { return velocityRowByNative_; }
    const std::vector<int>& nativeByVelocityRow() const noexcept { return nativeByVelocityRow_; }
    const std::vector<double>& cellVolumes() const noexcept { return geometry_.cellVolumes(); }
    const Eigen::SparseMatrix<double>& divergence() const noexcept { return divergence_; }
    const Eigen::SparseMatrix<double>& pressureGradientLoad() const noexcept { return pressureGradientLoad_; }
    bool axisymmetric() const noexcept { return geometry_.axisymmetric(); }
    bool uniformLeaves() const noexcept { return geometry_.uniformLeaves(); }
    Vector2D spacing() const; // Same-level callers only; mixed meshes reject.

    // Exact stored geometry/topology comparison, not a hash. Phase/cut flags,
    // flow fields and temporary solver interpolation caches are not identity.
    // Call before reusing this immutable domain after any possible mesh change.
    void validateMeshIdentity(const QuadtreeMesh& mesh) const;
    // Exact original NS domain-definition comparison (maps, faces and BCs).
    // Mesh identity is validated separately; no tolerance-based equivalence.
    bool sameDefinition(const HydroFittedDomain& other) const;

private:
    NativeLeafDomain geometry_;
    std::vector<NativeFace> faces_;
    std::vector<int> velocityRowByNative_, nativeByVelocityRow_;
    Eigen::SparseMatrix<double> divergence_, pressureGradientLoad_;
};
