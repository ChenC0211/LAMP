#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <cstdint>
#include "Core/Mesh/BackgroundMesh.h"

class QuadtreeMesh : public BackgroundMesh
{
public:    
    QuadtreeMesh(const DictionaryMesh& dictMesh);
    ~QuadtreeMesh();

    // Mesh arrays own their storage. An implicit shallow copy is never valid.
    QuadtreeMesh(const QuadtreeMesh&) = delete;
    QuadtreeMesh& operator=(const QuadtreeMesh&) = delete;

    // In-memory candidate for the EXISTING adapt() implementation. Copies
    // live leaves, parents, cut/derivative data, labels and refinement settings;
    // no background mesh construction, file IO, geometry reinitialization or
    // field interpolation is performed. Only clone-owned storage is allocated.
    // Fields/physical histories are caller-owned and must be staged separately.
    // Mutating/discarding the candidate leaves this accepted mesh unchanged.
    std::unique_ptr<QuadtreeMesh> clone() const;

    // Expand a native leaf-indexed mask by face-neighbor layers, including
    // all coarse/fine neighbors. Each layer reads only the previous layer.
    // Zero layers preserves the mask. Does not clone or change the mesh.
    void expandCellMask(std::vector<std::uint8_t>& mask, int layers) const;

    // Exchange owned storage only within the same background domain/lattice
    // and level-set count. All compatibility checks precede the non-allocating
    // exchange. References into BOTH meshes become stale; the caller must
    // commit matching fields/geometry/domains and invalidate solver caches.
    // This is not, by itself, acceptance of a coupled physical AMR step.
    void swapStorage(QuadtreeMesh& candidate);

    // Apply a list of static refine-box directives (if any) to create a mesh
    // with coarse--fine interfaces at construction time. This is used by demo
    // cases that need a guaranteed coarse--fine interface without running a
    // full adaptive refinement loop.
    void applyRefineBoxes(const std::vector<DictionaryMesh::RefineBox>& boxes);

    template<typename Field, typename VertField, typename HEdgeField, typename VEdgeField,
             typename CellFieldT, typename LevelCellFieldT = CellFieldT>
    void adapt(std::vector<Field>& fields, CellFieldT& magGradU, std::vector<VertField*>& vertFields, 
        std::vector<HEdgeField*>& hEdgeFields, std::vector<VEdgeField*>& vEdgeFields,
        std::vector<CellFieldT*>& cellFields, const LevelCellFieldT* targetLevel = nullptr);

    // Read-only preflight for adapt() with explicit targets and EMPTY level-
    // set inputs. False proves no topology change can start. Reuses the old
    // target clamping, 2:1 neighbor and four-sibling gates. True is conservative
    // (allocation limits may still veto refinement); never replaces adapt().
    // Requires finite targets in the current native cell-label ordering.
    template<typename LevelCellFieldT>
    bool mayAdaptToTargets(const LevelCellFieldT& targetLevel) const;

    // Mesh-only diagnostic: output/fields/mesh_<name>.dat.
    void output(const std::string& name) const;
    template<typename VertField>
    void outputFields(const std::string& name, const std::vector<std::string>& valNames, const std::vector<VertField>& fields) const;

    void writeCache(const std::string& filePath) const;
    // Read filePath + "mesh_cache.bin" (QMC2). The entire cache is staged and
    // validated, including non-leaf parent geometry, before replacing storage.
    // Invalid/truncated/incompatible caches throw std::runtime_error and leave
    // this mesh and the file untouched. Accepted geometry/cut/derivative values
    // are preserved; old inconsistent metrics are rejected, never repaired.
    void readCache(const std::string& filePath);

    inline int maxInitialIterations() const { return maxInitialIterations_; }
    inline int maxLevel() const { return maxLevel_; }

private:
    QuadtreeMesh(const QuadtreeMesh& source, CopyMetadataTag);

    int maxLevel_;

    int maxVertNum_;
    int maxHEdgeNum_;
    int maxVEdgeNum_;
    int maxCellNum_;

    int vertPtr_ = -1;
    int hEdgePtr_ = -1;
    int vEdgePtr_ = -1;
    int cellPtr_ = -1;

    int* vertLabelSpace_ = nullptr;
    int* hEdgeLabelSpace_ = nullptr;
    int* vEdgeLabelSpace_ = nullptr;
    int* cellLabelSpace_ = nullptr;

    int maxInitialIterations_ = 3;
    std::vector<int> minLevel_neg;
    std::vector<int> minLevel_pos;
    std::vector<int> maxLevel_local;

    double C_Lipschitz = 5.0;

    std::vector<int> nCellsPerLevel;

    int minimumAdaptationLevel() const;
    int boundedTargetLevel(double value, int minimum) const;
    bool targetRefinementAllowed(const Cell& cell) const;
    bool targetCoarseningAllowed(const Cell& cell) const;
    template<typename Eligible>
    bool coarseningFamilyEligible(const Cell& parent, const Eligible& eligible) const;

    // 内存扩容
    void resizeMemory();

    // 判定细化或粗化，以levelSetField为例
    template<typename VertField, typename CellFieldT, typename LevelCellFieldT = CellFieldT>
    void refineOrCoarsen0(std::vector<VertField>& LSs, CellFieldT& magGradU, const int currentLevel,
                          const LevelCellFieldT* targetLevel = nullptr);

    // Capture pre-transfer indices before each topology-changing level.
    // Shared by adaptive and static refinement; newly created children retain
    // -1 and resolve their source through their parent's saved index.
    // Interpolation stencils belong to this ONE pre-transfer index space:
    // discard older stencils before rebuilding, including surviving edges.
    void prepareTransferIndices();
    void subdivideCell(const int currentLevel);
    void submergeCell(const int currentLevel);

    void findSubdivideRelatedCells(Cell& cell);
    void findSubmergeRelatedCells(Cell& cell);

    void findSubdivideRelatedHEdges(const Cell& cell, Edge& hEdge0, Edge& hEdge1, Edge& hEdge2, const int direction, const bool isRightSide);
    void findSubdivideRelatedVEdges(const Cell& cell, Edge& vEdge0, Edge& vEdge1, Edge& vEdge2, const int direction, const bool isUpSide);
    void findSubmergeRelatedHEdges(Edge& hEdge, const Cell& sonCellA, const Cell& sonCellB);
    void findSubmergeRelatedVEdges(Edge& vEdge, const Cell& sonCellA, const Cell& sonCellB);

    void statistics();

    template<typename T>
    void updateLabelSpace(T*& list, int*& space, int& ptr, const int maxNum);

    template<typename T>
    void count(T*& list, const int ptr, int& indexNum, int*& indexToLabel);

};

#include "Core/Mesh/QuadtreeMesh_impl.h"
