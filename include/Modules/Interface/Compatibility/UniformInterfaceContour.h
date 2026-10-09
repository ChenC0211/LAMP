#pragma once
#include "Core/Lagrangian/InterfaceNode.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include <vector>

class QuadtreeMesh;
template<typename Mesh> class ImmersedGeometry;

// Owned snapshot of the existing negative-vertex closed-contour construction.
// Fixed uniform mesh only; this is not the fitted/AMR interface reconstruction.
// Nodes, curvature, measures and normalized cell support are read-only. No
// traction, material constants or fluid fields belong to this geometry object.
class UniformInterfaceContour
{
public:
    struct Report
    {
        int nodeCount = 0;
        double surfaceMeasure = 0.0;
        double maxAbsCurvature = 0.0;
    };
    static UniformInterfaceContour build(QuadtreeMesh& mesh,
        const ImmersedGeometry<QuadtreeMesh>& geometry, int levelSetIndex, bool axisymmetric);
    UniformInterfaceContour(const UniformInterfaceContour&) = delete;
    UniformInterfaceContour& operator=(const UniformInterfaceContour&) = delete;
    UniformInterfaceContour(UniformInterfaceContour&&) = default;
    UniformInterfaceContour& operator=(UniformInterfaceContour&&) = default;

    const std::vector<InterfaceNode>& nodes() const { return nodes_; }
    const Report& report() const { return report_; }
    bool axisymmetric() const { return axisymmetric_; }
    int levelSetIndex() const { return levelSetIndex_; }
    // The mesh must stay alive. Geometry/source fields may subsequently change:
    // this snapshot retains its original contour, rather than relabelling it.
    void validateMesh(const QuadtreeMesh& mesh) const;

private:
    UniformInterfaceContour(const QuadtreeMesh& mesh, int levelSetIndex, bool axisymmetric,
        std::vector<InterfaceNode> nodes, Report report);
    const QuadtreeMesh* mesh_;
    NativeMeshIdentity identity_;
    int levelSetIndex_;
    bool axisymmetric_;
    std::vector<InterfaceNode> nodes_;
    Report report_;
};
