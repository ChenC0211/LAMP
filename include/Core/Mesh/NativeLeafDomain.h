#pragma once

#include <memory>
#include <vector>
#include "Core/Base/Vector2D.h"

class FaceSegmentGraph;
class NativeMeshIdentity;
class QuadtreeMesh;

// Geometry of an actual native face, not a velocity/pressure unknown. The
// packed index uses VEdge first, then HEdge, as the existing staggered fields.
// Negative/positive cells are selected native cell indices, or -1 outside
// the selected domain. No wall condition, material or cut fraction is inferred.
struct NativeLeafFace
{
    int nativeRow = -1;
    int edgeLabel = -1;
    bool horizontal = false;
    int negativeCell = -1, positiveCell = -1;
    Vector2D center{0.0, 0.0};
    double area = 0.0;
    double normalSpacing = 0.0, tangentSpacing = 0.0;
    double outwardSign = 0.0; // Selected-domain boundary normal; interior 0.
};

// Selected complete Cartesian leaves on the existing quadtree. Owns compact
// cell indexing, measures and an exact mesh identity, never another mesh.
// XY measures are per unit depth; RZ measures include the full 2*pi*r factor.
// This is not a general cut-cell or free-surface representation.
class NativeLeafDomain
{
public:
    NativeLeafDomain() = default; // Empty until assigned build().
    static NativeLeafDomain build(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, bool axisymmetric);

    const std::vector<int>& cellRowByIndex() const noexcept { return cellRowByIndex_; }
    const std::vector<int>& cellIndexByRow() const noexcept { return cellIndexByRow_; }
    const std::vector<double>& cellVolumes() const noexcept { return cellVolumes_; }
    bool axisymmetric() const noexcept { return axisymmetric_; }
    bool uniformLeaves() const noexcept { return uniformLeaves_; }
    int nativeFaceCount() const noexcept { return verticalFaces_ + horizontalFaces_; }
    Vector2D spacing() const; // Mixed leaves have no single stencil width.

    // Constant-size query; no all-mesh scan/copy per face. Immediately after
    // build the mesh is verified. After ANY possible geometry/topology change,
    // first call validateMeshIdentity, including before a batch of queries.
    // The domain retains no borrowed mesh pointer and no duplicate face array.
    NativeLeafFace face(const QuadtreeMesh& mesh, int nativeRow) const;
    // Borrowed native-row view for geometry consumers. Iteration calls the
    // existing constant-size face query; no face array or mesh is copied.
    // The domain and mesh must outlive the view, remain unchanged, and satisfy
    // the same validateMeshIdentity precondition as face().
    class FaceView
    {
    public:
        struct Iterator
        {
            const FaceView* view;
            int row;
            NativeLeafFace operator*() const { return (*view)[row]; }
            Iterator& operator++() { ++row; return *this; }
            bool operator!=(const Iterator& other) const { return row != other.row; }
        };
        NativeLeafFace operator[](int row) const { return domain_->face(*mesh_, row); }
        std::size_t size() const { return static_cast<std::size_t>(domain_->nativeFaceCount()); }
        Iterator begin() const { return {this, 0}; }
        Iterator end() const { return {this, domain_->nativeFaceCount()}; }
    private:
        friend class NativeLeafDomain;
        FaceView(const NativeLeafDomain& domain, const QuadtreeMesh& mesh) : domain_(&domain), mesh_(&mesh) {}
        const NativeLeafDomain* domain_;
        const QuadtreeMesh* mesh_;
    };
    FaceView faces(const QuadtreeMesh& mesh) const { return FaceView(*this, mesh); }
    void validateMeshIdentity(const QuadtreeMesh& mesh) const;
    // Shared immutable construction identity for retained native-field sources.
    // Domain copies share it; rebuilding/adapting a domain creates a new one.
    std::shared_ptr<const NativeMeshIdentity> meshIdentity() const noexcept { return identity_; }

private:
    int activeCellIndex(const QuadtreeMesh& mesh, int label) const;
    std::shared_ptr<const NativeMeshIdentity> identity_;
    bool axisymmetric_ = false;
    bool uniformLeaves_ = true;
    Vector2D spacing_{0.0, 0.0};
    int verticalFaces_ = 0, horizontalFaces_ = 0;
    std::vector<int> cellRowByIndex_, cellIndexByRow_;
    std::vector<double> cellVolumes_;
};
