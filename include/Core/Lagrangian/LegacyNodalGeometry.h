#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include "Core/Lagrangian/InterfaceNode.h"
#include "Core/Geometry/LevelSet.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Base/Utility.h"

class LegacyNodalGeometry
{
public:
    std::vector<std::vector<int>> nodeIndexToLabel;
    std::vector<std::vector<int>> preNodeIndexToLabel;
    std::vector<std::vector<int>> prePreNodeIndexToLabel;
    int noSlipNodeNum = 0;

    LegacyNodalGeometry() = default;
    ~LegacyNodalGeometry() = default;

    // Legacy node construction, labels, geometry history and interpolation stencils.
    template<typename Mesh>
    void initializeGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry);
    // Advance one level set; the Hydro owner controls velocity/history ordering.
    template<typename Mesh>
    void updateGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI);


    InterfaceNode& interfaceNode(const int levelSetI, const int nodeLabel) { return allNodes_[levelSetI][nodeLabel]; }
    InterfaceNode& preInterfaceNode(const int levelSetI, const int nodeLabel) { return allPreNodes_[levelSetI][nodeLabel]; }
    InterfaceNode& prePreInterfaceNode(const int levelSetI, const int nodeLabel) { return allPrePreNodes_[levelSetI][nodeLabel]; }
    InterfaceNode& noSlipNode(const int nodeLabel) { return noSlipNodes_[nodeLabel]; }
    std::vector<InterfaceNode>& allNodes(const int levelSetI) { return allNodes_[levelSetI]; }
    std::vector<InterfaceNode>& allPreNodes(const int levelSetI) { return allPreNodes_[levelSetI]; }
    std::vector<InterfaceNode>& allPrePreNodes(const int levelSetI) { return allPrePreNodes_[levelSetI]; }
    std::vector<InterfaceNode>& noSlipNodes() { return noSlipNodes_; }


    void resize(int n)
    {
        const int count = std::max(0, n);
        nodeIndexToLabel.resize(count);
        preNodeIndexToLabel.resize(count);
        prePreNodeIndexToLabel.resize(count);
        allNodes_.resize(count);
        allPreNodes_.resize(count);
        allPrePreNodes_.resize(count);
    }

    const InterfaceNode& interfaceNode(const int levelSetI, const int nodeLabel) const { return allNodes_[levelSetI][nodeLabel]; }
    const InterfaceNode& preInterfaceNode(const int levelSetI, const int nodeLabel) const { return allPreNodes_[levelSetI][nodeLabel]; }
    const InterfaceNode& prePreInterfaceNode(const int levelSetI, const int nodeLabel) const { return allPrePreNodes_[levelSetI][nodeLabel]; }
    const InterfaceNode& noSlipNode(const int nodeLabel) const { return noSlipNodes_[nodeLabel]; }
    const std::vector<InterfaceNode>& allNodes(const int levelSetI) const { return allNodes_[levelSetI]; }
    const std::vector<InterfaceNode>& allPreNodes(const int levelSetI) const { return allPreNodes_[levelSetI]; }
    const std::vector<InterfaceNode>& allPrePreNodes(const int levelSetI) const { return allPrePreNodes_[levelSetI]; }
    const std::vector<InterfaceNode>& noSlipNodes() const { return noSlipNodes_; }

private:
    std::vector<std::vector<InterfaceNode>> allNodes_;
    std::vector<std::vector<InterfaceNode>> allPreNodes_;
    std::vector<std::vector<InterfaceNode>> allPrePreNodes_;
    std::vector<InterfaceNode> noSlipNodes_;

    template<typename Mesh>
    void createNodes(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI);
    template<typename Mesh>
	void findNeighborNodes(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI);
    template<typename Mesh>
	void calculateGeometricFeatures(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI);
    template<typename Mesh>
    void createNoSlipNodes(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry);

    template<typename Mesh>
    void setIgnore(ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, int& _noSlipNodeNum);
    template<typename Mesh>
    void checkLabel(Mesh& mesh, const int levelSetI);

    template<typename Mesh>
    void findStencil(Mesh& mesh, InterfaceNode& node);
    template<typename Mesh>
    void buildStencil(Mesh& mesh, const int levelSetI);
    template<typename Mesh>
    void buildNoSlipStencil(Mesh& mesh);
};

#include "Core/Lagrangian/LegacyNodalGeometry_impl.h"

