#pragma once

#include "Core/Base/Enumeration.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Lagrangian/LegacyNodalGeometry.h"

using namespace std;

template<typename Mesh>
void LegacyNodalGeometry::initializeGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry)
{

    resize(static_cast<int>(iGeometry.levelSets.size()));
    noSlipNodeNum = 0;
    for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); levelSetI++)
    {
        createNodes(mesh, iGeometry, levelSetI);
        // setIgnore(iGeometry, levelSetI, noSlipNodeNum);
        checkLabel(mesh, levelSetI);
        findNeighborNodes(mesh, iGeometry, levelSetI);
        calculateGeometricFeatures(mesh, iGeometry, levelSetI);
        if (iGeometry.levelSets[levelSetI].moving())
        {
            buildStencil(mesh, levelSetI);
        }
    }
    // createNoSlipNodes(mesh, iGeometry);
    // buildNoSlipStencil(mesh);

}

template<typename Mesh>
void LegacyNodalGeometry::updateGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI)
{
    allPrePreNodes_[levelSetI] = allPreNodes_[levelSetI];
    prePreNodeIndexToLabel[levelSetI] = preNodeIndexToLabel[levelSetI];

    allPreNodes_[levelSetI] = allNodes_[levelSetI];
    preNodeIndexToLabel[levelSetI] = nodeIndexToLabel[levelSetI];

    if (iGeometry.levelSets[levelSetI].moving())
    {
        createNodes(mesh, iGeometry, levelSetI);
        checkLabel(mesh, levelSetI);
        findNeighborNodes(mesh, iGeometry, levelSetI);
        calculateGeometricFeatures(mesh, iGeometry, levelSetI);
        buildStencil(mesh, levelSetI);
    }
}


template<typename Mesh>
void LegacyNodalGeometry::createNodes(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI)
{
    const LevelSet& levelSet = iGeometry.levelSets[levelSetI];

    allNodes_[levelSetI].clear();
    allNodes_[levelSetI].resize(mesh.negFaceVerts[levelSetI].size());
    // 负界面节点
    for (int i = 0; i < mesh.negFaceVerts[levelSetI].size(); i++)
    {
        const int vertL = mesh.negFaceVerts[levelSetI][i];
        Vert& vert = mesh.verts(vertL);

        InterfaceNode& node = allNodes_[levelSetI][i];
        node.setCoord(vert.nearestNodeCoord[levelSetI]);
        node.setIndex(i);
        node.setLabel(i);

        node.setRelatedNegVertLabel(vert.label());
        vert.nearestNodeLabel[levelSetI] = i;
    }

    // 正界面节点
    for (int i = 0; i < mesh.posFaceVerts[levelSetI].size(); i++)
    {
        const int vertI = mesh.posFaceVerts[levelSetI][i];
        Vert& vert = mesh.verts(vertI);

        for (int j = 0; j < 2; j++)
        {
            if (vert.vertToHEdge[j] != -1)
            {
                Edge& hEdge = mesh.hEdges(vert.vertToHEdge[j]);
                if (hEdge.isCut[levelSetI])
                {
                    Vert& neiVert = mesh.verts(hEdge.edgeToVert[j]);
                    if (vert.nearestNodeLabel[levelSetI] == -1)
                    {
                        vert.nearestNodeLabel[levelSetI] = neiVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = neiVert.nearestNodeCoord[levelSetI];
                    }
                    else if ((vert.coord() - neiVert.nearestNodeCoord[levelSetI]).mod2() < (vert.coord() - vert.nearestNodeCoord[levelSetI]).mod2())
                    {
                        vert.nearestNodeLabel[levelSetI] = neiVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = neiVert.nearestNodeCoord[levelSetI];
                    }
                }
            }
            else if (vert.supCell != -1)
            {
                vector<int> supVertLabels;
                Vert& supVert0 = mesh.verts(vert.supVert[0]);
                if (levelSet.value[supVert0.index()] <= 0.0 && supVert0.nearestNodeLabel[levelSetI] != -1)
                {
                    supVertLabels.push_back(vert.supVert[0]);
                }
                Vert& supVert1 = mesh.verts(vert.supVert[1]);
                if (levelSet.value[supVert1.index()] <= 0.0 && supVert1.nearestNodeLabel[levelSetI] != -1)
                {
                    supVertLabels.push_back(vert.supVert[1]);
                }
                if (mesh.cells(vert.supCell).cellToHEdge[0][1] != -1)
                {
                    Vert& southVert = mesh.verts(mesh.hEdges(mesh.cells(vert.supCell).cellToHEdge[0][1]).edgeToVert[0]);
                    if (levelSet.value[southVert.index()] <= 0.0 && southVert.nearestNodeLabel[levelSetI] != -1)
                    {
                        supVertLabels.push_back(southVert.label());
                    }
                }
                if (mesh.cells(vert.supCell).cellToHEdge[1][1] != -1)
                {
                    Vert& northVert = mesh.verts(mesh.hEdges(mesh.cells(vert.supCell).cellToHEdge[1][1]).edgeToVert[0]);
                    if (levelSet.value[northVert.index()] <= 0.0 && northVert.nearestNodeLabel[levelSetI] != -1)
                    {
                        supVertLabels.push_back(northVert.label());
                    }
                }

                if (vert.type() == TVertType::west)
                {
                    if (mesh.cells(vert.supCell).cellToVEdge[0][1] != -1)
                    {
                        Vert& westVert = mesh.verts(mesh.vEdges(mesh.cells(vert.supCell).cellToVEdge[0][1]).edgeToVert[0]);
                        if (levelSet.value[westVert.index()] <= 0.0 && westVert.nearestNodeLabel[levelSetI] != -1)
                        {
                            supVertLabels.push_back(westVert.label());
                        }
                    } 
                }
                else if (vert.type() == TVertType::east)
                {
                    if (mesh.cells(vert.supCell).cellToVEdge[1][1] != -1)
                    {
                        Vert& eastVert = mesh.verts(mesh.vEdges(mesh.cells(vert.supCell).cellToVEdge[1][1]).edgeToVert[0]);
                        if (levelSet.value[eastVert.index()] <= 0.0 && eastVert.nearestNodeLabel[levelSetI] != -1)
                        {
                            supVertLabels.push_back(eastVert.label());
                        }
                    }
                }

                for (int v = 0; v < supVertLabels.size(); v++)
                {
                    Vert& supVert = mesh.verts(supVertLabels[v]);
                    if (vert.nearestNodeLabel[levelSetI] == -1)
                    {
                        vert.nearestNodeLabel[levelSetI] = supVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = supVert.nearestNodeCoord[levelSetI];
                    }
                    else if ((vert.coord() - supVert.nearestNodeCoord[levelSetI]).mod2() < (vert.coord() - vert.nearestNodeCoord[levelSetI]).mod2())
                    {
                        vert.nearestNodeLabel[levelSetI] = supVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = supVert.nearestNodeCoord[levelSetI];
                    }
                }
            }
            
            if (vert.vertToVEdge[j] != -1)
            {
                Edge& vEdge = mesh.vEdges(vert.vertToVEdge[j]);
                if (vEdge.isCut[levelSetI])
                {
                    Vert& neiVert = mesh.verts(vEdge.edgeToVert[j]);
                    if (vert.nearestNodeLabel[levelSetI] == -1)
                    {
                        vert.nearestNodeLabel[levelSetI] = neiVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = neiVert.nearestNodeCoord[levelSetI];
                    }
                    else if ((vert.coord() - neiVert.nearestNodeCoord[levelSetI]).mod2() < (vert.coord() - vert.nearestNodeCoord[levelSetI]).mod2())
                    {
                        vert.nearestNodeLabel[levelSetI] = neiVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = neiVert.nearestNodeCoord[levelSetI];
                    }
                }
            }
            else if (vert.supCell != -1)
            {
                vector<int> supVertLabels;
                Vert& supVert0 = mesh.verts(vert.supVert[0]);
                if (levelSet.value[supVert0.index()] <= 0.0 && supVert0.nearestNodeLabel[levelSetI] != -1)
                {
                    supVertLabels.push_back(vert.supVert[0]);
                }
                Vert& supVert1 = mesh.verts(vert.supVert[1]);
                if (levelSet.value[supVert1.index()] <= 0.0 && supVert1.nearestNodeLabel[levelSetI] != -1)
                {
                    supVertLabels.push_back(vert.supVert[1]);
                }
                if (mesh.cells(vert.supCell).cellToVEdge[0][1] != -1)
                {
                    Vert& westVert = mesh.verts(mesh.vEdges(mesh.cells(vert.supCell).cellToVEdge[0][1]).edgeToVert[0]);
                    if (levelSet.value[westVert.index()] <= 0.0 && westVert.nearestNodeLabel[levelSetI] != -1)
                    {
                        supVertLabels.push_back(westVert.label());
                    }
                }
                if (mesh.cells(vert.supCell).cellToVEdge[1][1] != -1)
                {
                    Vert& eastVert = mesh.verts(mesh.vEdges(mesh.cells(vert.supCell).cellToVEdge[1][1]).edgeToVert[0]);
                    if (levelSet.value[eastVert.index()] <= 0.0 && eastVert.nearestNodeLabel[levelSetI] != -1)
                    {
                        supVertLabels.push_back(eastVert.label());
                    }
                }

                if (vert.type() == TVertType::south)
                {
                    if (mesh.cells(vert.supCell).cellToHEdge[0][1] != -1)
                    {
                        Vert& southVert = mesh.verts(mesh.hEdges(mesh.cells(vert.supCell).cellToHEdge[0][1]).edgeToVert[0]);
                        if (levelSet.value[southVert.index()] <= 0.0 && southVert.nearestNodeLabel[levelSetI] != -1)
                        {
                            supVertLabels.push_back(southVert.label());
                        }
                    }
                }
                else if (vert.type() == TVertType::north)
                {
                    if (mesh.cells(vert.supCell).cellToHEdge[1][1] != -1)
                    {
                        Vert& northVert = mesh.verts(mesh.hEdges(mesh.cells(vert.supCell).cellToHEdge[1][1]).edgeToVert[0]);
                        if (levelSet.value[northVert.index()] <= 0.0 && northVert.nearestNodeLabel[levelSetI] != -1)
                        {
                            supVertLabels.push_back(northVert.label());
                        }
                    }
                }

                for (int v = 0; v < supVertLabels.size(); v++)
                {
                    Vert& supVert = mesh.verts(supVertLabels[v]);
                    if (vert.nearestNodeLabel[levelSetI] == -1)
                    {
                        vert.nearestNodeLabel[levelSetI] = supVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = supVert.nearestNodeCoord[levelSetI];
                    }
                    else if ((vert.coord() - supVert.nearestNodeCoord[levelSetI]).mod2() < (vert.coord() - vert.nearestNodeCoord[levelSetI]).mod2())
                    {
                        vert.nearestNodeLabel[levelSetI] = supVert.nearestNodeLabel[levelSetI];
                        vert.nearestNodeCoord[levelSetI] = supVert.nearestNodeCoord[levelSetI];
                    }
                }
            }
        }
    }
}

template<typename Mesh>
void LegacyNodalGeometry::findNeighborNodes(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI)
{
    const LevelSet& levelSet = iGeometry.levelSets[levelSetI];
    for (int cellI = 0; cellI < mesh.cellIndexNum(); cellI++)
    {
        Cell& cell = mesh.cells(mesh.cellIndexToLabel[cellI]);
        if (cell.isCut[levelSetI])
        {
            vector<int> cutVertLabels;
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    if (cell.cellToVEdge[i][j] != -1)
                    {
                        Edge& vEdge = mesh.vEdges(cell.cellToVEdge[i][j]);
                        if (vEdge.isCut[levelSetI])
                        {
                            Vert& vert0 = mesh.verts(vEdge.edgeToVert[0]);
                            Vert& vert1 = mesh.verts(vEdge.edgeToVert[1]);
                            if (levelSet.value[vert0.index()] <= 0.0)
                            {
                                if (vert0.nearestNodeLabel[levelSetI] != -1)
                                {
                                    cutVertLabels.push_back(vert0.label());
                                }
                            }
                            if (levelSet.value[vert1.index()] <= 0.0)
                            {
                                if (vert1.nearestNodeLabel[levelSetI] != -1)
                                {
                                    cutVertLabels.push_back(vert1.label());
                                }
                            }
                        }
                    }

                    if (cell.cellToHEdge[i][j] != -1)
                    {
                        Edge& hEdge = mesh.hEdges(cell.cellToHEdge[i][j]);
                        if (hEdge.isCut[levelSetI])
                        {
                            Vert& vert0 = mesh.verts(hEdge.edgeToVert[0]);
                            Vert& vert1 = mesh.verts(hEdge.edgeToVert[1]);
                            if (levelSet.value[vert0.index()] <= 0.0)
                            {
                                if (vert0.nearestNodeLabel[levelSetI] != -1)
                                {
                                    cutVertLabels.push_back(vert0.label());
                                }
                            }
                            if (levelSet.value[vert1.index()] <= 0.0)
                            {
                                if (vert1.nearestNodeLabel[levelSetI] != -1)
                                {
                                    cutVertLabels.push_back(vert1.label());
                                }
                            }
                        }
                    }
                }
            }

            if (cutVertLabels.size() < 2)
            {
                continue;
            }
            else if (cutVertLabels.size() > 2)
            {
                // 双界面，后面处理
            }
            else if (cutVertLabels[0] == cutVertLabels[1])
            {
                // 0 1 2 3 (无T型点)
                // jump
                continue;
            }

            if (levelSet.value[mesh.verts(cell.cellToVert[0]).index()] <= 0.0)
            {
                if (levelSet.value[mesh.verts(cell.cellToVert[1]).index()] <= 0.0)
                {
                    if (levelSet.value[mesh.verts(cell.cellToVert[2]).index()] <= 0.0)
                    {
                        // 0-1-2
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                        node0.setUpNodeLabel(node1.label());
                        node1.setDownNodeLabel(node0.label());

                        if (cell.cellToHEdge[0][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.hEdges(cell.cellToHEdge[0][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setUpNodeLabel(node1.label());
                        }
                        if (cell.cellToVEdge[1][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.vEdges(cell.cellToVEdge[1][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setDownNodeLabel(node0.label());
                        }
                    }
                    else if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] <= 0.0)
                    {
                        // 0-1-3
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());

                        if (cell.cellToVEdge[0][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.vEdges(cell.cellToVEdge[0][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setUpNodeLabel(node0.label());
                        }
                        if (cell.cellToHEdge[0][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.hEdges(cell.cellToHEdge[0][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setDownNodeLabel(node1.label());
                        }
                    }
                    else
                    {
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                        // 0-1
                        if (cell.cellToHEdge[0][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.hEdges(cell.cellToHEdge[0][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node0.setUpNodeLabel(node2.label());
                            node2.setDownNodeLabel(node0.label());

                            node1.setDownNodeLabel(node2.label());
                            node2.setUpNodeLabel(node1.label());
                        }
                        else
                        {
                            node0.setUpNodeLabel(node1.label());
                            node1.setDownNodeLabel(node0.label());
                        }
                    }
                }
                else if (levelSet.value[mesh.verts(cell.cellToVert[2]).index()] <= 0.0)
                {
                    if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] < 0.0)
                    {
                        // 0-2-3
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());

                        if (cell.cellToVEdge[0][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.vEdges(cell.cellToVEdge[0][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setDownNodeLabel(node1.label());
                        }
                        if (cell.cellToHEdge[1][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.hEdges(cell.cellToHEdge[1][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setUpNodeLabel(node0.label());
                        }
                    }
                    else
                    {
                        // 0-2
                        // 双界面处理
                        if (cutVertLabels.size() != 4)
                        {
                            cerr << "[SharpInterface] ERROR: invalid interface info for case 0-2." << endl;
                            exit(1);
                        }

                        if (cutVertLabels[0] == cutVertLabels[1])
                        {
                            if (cutVertLabels[2] == cutVertLabels[3])
                            {
                                continue;
                            }
                            else
                            {
                                InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[2]).nearestNodeLabel[levelSetI]];
                                InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[3]).nearestNodeLabel[levelSetI]];

                                node0.setUpNodeLabel(node1.label());
                                node1.setDownNodeLabel(node0.label());
                            }
                        }
                        else if (cutVertLabels[2] == cutVertLabels[3])
                        {
                            InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                            InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                            node0.setUpNodeLabel(node1.label());
                            node1.setDownNodeLabel(node0.label());
                        }
                        else
                        {
                            InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                            InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                            node0.setUpNodeLabel(node1.label());
                            node1.setDownNodeLabel(node0.label());

                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabels[2]).nearestNodeLabel[levelSetI]];
                            InterfaceNode& node3 = allNodes_[levelSetI][mesh.verts(cutVertLabels[3]).nearestNodeLabel[levelSetI]];

                            node2.setUpNodeLabel(node3.label());
                            node3.setDownNodeLabel(node2.label());
                        }
                    }
                }
                else if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] <= 0.0)
                {
                    // 0-3
                    InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                    InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                    if (cell.cellToVEdge[0][1] != -1)// 存在T型点
                    {
                        int cutVertLabel2 = mesh.vEdges(cell.cellToVEdge[0][1]).edgeToVert[0];
                        InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node2.label());
                        node2.setUpNodeLabel(node0.label());

                        node1.setUpNodeLabel(node2.label());
                        node2.setDownNodeLabel(node1.label());
                    }
                    else
                    {
                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());
                    }
                }
                else
                {
                    // 0 存在T型点
                    InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                    InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                    node0.setUpNodeLabel(node1.label());
                    node1.setDownNodeLabel(node0.label());
                }
            }
            else if (levelSet.value[mesh.verts(cell.cellToVert[1]).index()] <= 0.0)
            {
                if (levelSet.value[mesh.verts(cell.cellToVert[2]).index()] <= 0.0)
                {
                    if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] <= 0.0)
                    {
                        // 1-2-3
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());

                        if (cell.cellToVEdge[1][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.vEdges(cell.cellToVEdge[1][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setUpNodeLabel(node0.label());
                        }
                        if (cell.cellToHEdge[1][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.hEdges(cell.cellToHEdge[1][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node2.setDownNodeLabel(node1.label());
                        }
                    }
                    else
                    {
                        // 1-2
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                        if (cell.cellToVEdge[1][1] != -1)// 存在T型点
                        {
                            int cutVertLabel2 = mesh.vEdges(cell.cellToVEdge[1][1]).edgeToVert[0];
                            InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                            node0.setUpNodeLabel(node2.label());
                            node2.setDownNodeLabel(node0.label());

                            node1.setDownNodeLabel(node2.label());
                            node2.setUpNodeLabel(node1.label());
                        }
                        else
                        {
                            node0.setUpNodeLabel(node1.label());
                            node1.setDownNodeLabel(node0.label());
                        }
                    }
                }
                else if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] <= 0.0)
                {
                    // 1-3
                    // 双界面处理
                    if (cutVertLabels.size() != 4)
                    {
                        cerr << "[SharpInterface] ERROR: invalid interface info for case 1-3." << endl;
                        exit(1);
                    }

                    if (cutVertLabels[0] == cutVertLabels[3])
                    {
                        if (cutVertLabels[1] == cutVertLabels[2])
                        {
                            continue;
                        }
                        else
                        {
                            InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];
                            InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[2]).nearestNodeLabel[levelSetI]];

                            node0.setUpNodeLabel(node1.label());
                            node1.setDownNodeLabel(node0.label());
                        }
                    }
                    else if (cutVertLabels[1] == cutVertLabels[2])
                    {
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[3]).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());
                    }
                    else
                    {
                        InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[3]).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());

                        InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];
                        InterfaceNode& node3 = allNodes_[levelSetI][mesh.verts(cutVertLabels[2]).nearestNodeLabel[levelSetI]];

                        node2.setUpNodeLabel(node3.label());
                        node3.setDownNodeLabel(node2.label());
                    }
                }
                else
                {
                    // 1 存在T型点
                    InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                    InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                    node0.setUpNodeLabel(node1.label());
                    node1.setDownNodeLabel(node0.label());
                }
            }
            else if (levelSet.value[mesh.verts(cell.cellToVert[2]).index()] <= 0.0)
            {
                if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] <= 0.0)
                {
                    // 2-3
                    InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                    InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                    if (cell.cellToHEdge[1][1] != -1)// 存在T型点
                    {
                        int cutVertLabel2 = mesh.hEdges(cell.cellToHEdge[1][1]).edgeToVert[0];
                        InterfaceNode& node2 = allNodes_[levelSetI][mesh.verts(cutVertLabel2).nearestNodeLabel[levelSetI]];

                        node0.setDownNodeLabel(node2.label());
                        node2.setUpNodeLabel(node0.label());

                        node1.setUpNodeLabel(node2.label());
                        node2.setDownNodeLabel(node1.label());
                    }
                    else
                    {
                        node0.setDownNodeLabel(node1.label());
                        node1.setUpNodeLabel(node0.label());
                    }
                }
                else
                {
                    // 2 存在T型点
                    InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                    InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                    node0.setUpNodeLabel(node1.label());
                    node1.setDownNodeLabel(node0.label());
                }
            }
            else if (levelSet.value[mesh.verts(cell.cellToVert[3]).index()] <= 0.0)
            {
                // 3 存在T型点
                InterfaceNode& node0 = allNodes_[levelSetI][mesh.verts(cutVertLabels[0]).nearestNodeLabel[levelSetI]];
                InterfaceNode& node1 = allNodes_[levelSetI][mesh.verts(cutVertLabels[1]).nearestNodeLabel[levelSetI]];

                node0.setDownNodeLabel(node1.label());
                node1.setUpNodeLabel(node0.label());
            }
        }
    }

    for (int i = 0; i < nodeIndexToLabel[levelSetI].size(); i++)
    {
        const int nodeLalbe = nodeIndexToLabel[levelSetI][i];
        InterfaceNode& node = allNodes_[levelSetI][nodeLalbe];

        // 修正T型点对应节点的相邻节点
        if (node.upNodeLabel() == -1)
        {
            InterfaceNode& downNode = allNodes_[levelSetI][node.downNodeLabel()];
            node.setUpNodeLabel(allNodes_[levelSetI][downNode.upNodeLabel()].upNodeLabel());
        }
        if (node.downNodeLabel() == -1)
        {
            InterfaceNode& upNode = allNodes_[levelSetI][node.upNodeLabel()];
            node.setDownNodeLabel(allNodes_[levelSetI][upNode.downNodeLabel()].downNodeLabel());
        }

        // 查找相关节点，起始点和终点特殊处理
        node.relatedNodes[0] = node.label();
        if (node.upNodeLabel() == -1)
        {
            node.relatedNodes[1] = node.downNodeLabel();
            node.relatedNodes[2] = allNodes_[levelSetI][node.downNodeLabel()].downNodeLabel();
        }
        else if (node.downNodeLabel() == -1)
        {
            node.relatedNodes[1] = node.upNodeLabel();
            node.relatedNodes[2] = allNodes_[levelSetI][node.upNodeLabel()].upNodeLabel();
        }
        else
        {
            node.relatedNodes[1] = node.upNodeLabel();
            node.relatedNodes[2] = node.downNodeLabel();
        }
    }
}

template<typename Mesh>
void LegacyNodalGeometry::calculateGeometricFeatures(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI)
{
    const LevelSet& levelSet = iGeometry.levelSets[levelSetI];
    // ds
    for (int i = 0; i < nodeIndexToLabel[levelSetI].size(); i++)
    {
        const int nodeLabel = nodeIndexToLabel[levelSetI][i];
        InterfaceNode& node = allNodes_[levelSetI][nodeLabel];

        if (node.upNodeLabel() != -1)
        {
            InterfaceNode& upNode = allNodes_[levelSetI][node.upNodeLabel()];
            node.setDs((upNode.coord() - node.coord()).mod());
        }
        else
        {
            node.setDs(0.0);
        }
    }

    for (int i = 0; i < nodeIndexToLabel[levelSetI].size(); i++)
    {
        const int nodeLabel = nodeIndexToLabel[levelSetI][i];
        InterfaceNode& node = allNodes_[levelSetI][nodeLabel];

        double _downDS = 0.0;
        if (node.downNodeLabel() != -1)
        {
            InterfaceNode& downNode = allNodes_[levelSetI][node.downNodeLabel()];
            _downDS = downNode.ds();
        }
        node.setDs(0.5 * (node.ds() + _downDS));
    }

    // normal
    for (int i = 0; i < nodeIndexToLabel[levelSetI].size(); i++)
    {
        const int nodeLabel = nodeIndexToLabel[levelSetI][i];
        InterfaceNode& node = allNodes_[levelSetI][nodeLabel];

        Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));

        const double dx = Interpolation::biQuadInterpolation(mesh, iGeometry.LS_dx(levelSetI), cell, node.coord());
        const double dy = Interpolation::biQuadInterpolation(mesh, iGeometry.LS_dy(levelSetI), cell, node.coord());
        const double dxx = Interpolation::biLinearInterpolation(mesh, iGeometry.LS_dxx(levelSetI), cell, node.coord());
        const double dyy = Interpolation::biLinearInterpolation(mesh, iGeometry.LS_dyy(levelSetI), cell, node.coord());
        const double dxy = Interpolation::biLinearInterpolation(mesh, iGeometry.LS_dxy(levelSetI), cell, node.coord());

        const double _mod2 = mod2(dx, dy);
        const double mod_inv = 1.0 / sqrt(_mod2);
        node.setNormal(Vector2D(dx * mod_inv, dy * mod_inv));

        const double dx_n1 = dxx * mod_inv + dx * (-0.5 * pow(_mod2, -1.5)) * (2.0 * (dx * dxx + dy * dxy));
        const double dy_n2 = dyy * mod_inv + dy * (-0.5 * pow(_mod2, -1.5)) * (2.0 * (dx * dxy + dy * dyy));

        node.setCurvatureInPlane(-(dx_n1 + dy_n2));
        double systemCurvature = -node.normal().x / node.coord().x;
        if (fabs(node.coord().x) < 0.1 * cell.size().x)// 对称轴附近特殊处理
        {
            systemCurvature = node.curvatureInPlane();
        }

        // systemCurvature = 0.0;
        node.setCurvature(systemCurvature + node.curvatureInPlane());
    }
}



template<typename Mesh>
void LegacyNodalGeometry::setIgnore(ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, int& _noSlipNodeNum)
{
    // 找出节点对应的geometry.dict中对应的Line
    const LevelSet& levelSet = iGeometry.levelSets[levelSetI];
    for (int i = 0; i < allNodes_[levelSetI].size(); i++)
    {
        InterfaceNode& node = allNodes_[levelSetI][i];

        int nearestLineI = -1;
        double minDist2 = std::numeric_limits<double>::max();
        double bestT = 0.0;
        for (int localLineI = 0; localLineI < levelSet.lines().size(); localLineI++)
        {
            const int lineI = levelSet.lines()[localLineI].value;
            vector<Vector2D> discretePoints = iGeometry.discreteLines(levelSetI, localLineI);

            BSplineCurve bspline(discretePoints, 3);

            for (double t = 0; t <= 1; t += 0.00025)
            {
                Vector2D splinePoint = bspline.evaluate(t);
                double dist2 = (splinePoint-node.coord()).mod2();
                if (dist2 < minDist2)
                {
                    minDist2 = dist2;
                    bestT = t;
                    nearestLineI = lineI;
                }
            }

            // Refine the closest point search
            double refinedMinDist2 = minDist2;
            double tStart = (bestT - 0.001 > 0) ? bestT - 0.001 : 0;
            double tEnd = (bestT + 0.001 < 1) ? bestT + 0.001 : 1;
            for (double t = tStart; t <= tEnd; t += 0.0001)
            {
                Vector2D splinePoint = bspline.evaluate(t);
                double dist2 = (splinePoint-node.coord()).mod2();
                if (dist2 < refinedMinDist2)
                {
                    refinedMinDist2 = dist2;
                    bestT = t;
                    nearestLineI = lineI;
                }
            }
        }

        if (nearestLineI == -1)
        {
            cerr << "[SharpInterface] ERROR: failed to find nearest line." << endl;
            exit(1);
        }

        node.setLineIndex(nearestLineI);

        if (iGeometry.lines[node.lineIndex()].ignored())
        {
            node.setIndex(-1);
        }
        else if (iGeometry.lines[node.lineIndex()].geoBoundary() == GeometricBoundary::fixed)
        {
            _noSlipNodeNum++;
        }
    }
}

template<typename Mesh>
void LegacyNodalGeometry::checkLabel(Mesh& mesh, const int levelSetI)
{
    // 筛选掉被忽略的节点
    nodeIndexToLabel[levelSetI].clear();
    for (int i = 0; i < allNodes_[levelSetI].size(); i++)
    {
        InterfaceNode& node = allNodes_[levelSetI][i];

        if (node.index() != -1)
        {
            nodeIndexToLabel[levelSetI].push_back(i);
        }
    }

    for (int i = 0; i < mesh.posFaceVerts[levelSetI].size(); i++)
    {
        Vert& vert = mesh.verts(mesh.posFaceVerts[levelSetI][i]);
        InterfaceNode& node = allNodes_[levelSetI][vert.nearestNodeLabel[levelSetI]];
        if (node.index() == -1)
        {
            vert.nearestNodeLabel[levelSetI] = -1;
        }
    }
    for (int i = 0; i < mesh.negFaceVerts[levelSetI].size(); i++)
    {
        Vert& vert = mesh.verts(mesh.negFaceVerts[levelSetI][i]);
        InterfaceNode& node = allNodes_[levelSetI][i];
        if (node.index() == -1)
        {
            vert.nearestNodeLabel[levelSetI] = -1;
        }
    }
}

template<typename Mesh>
void LegacyNodalGeometry::createNoSlipNodes(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry)
{
    noSlipNodes_.resize(noSlipNodeNum);
    int nodeI = 0;
    for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); levelSetI++)
    {
        for (int i = 0; i < nodeIndexToLabel[levelSetI].size(); i++)
        {
            const int nodeLabel = nodeIndexToLabel[levelSetI][i];
            InterfaceNode node = allNodes_[levelSetI][nodeLabel];

            if (iGeometry.lines[node.lineIndex()].geoBoundary() == GeometricBoundary::fixed)
            {
                noSlipNodes_[nodeI] = node;
                noSlipNodes_[nodeI].setIndex(nodeI);
                noSlipNodes_[nodeI].setLabel(nodeI);

                nodeI++;
            }
        }
    }
    if (nodeI != noSlipNodeNum)
    {
        cerr << "[SharpInterface] ERROR: noSlipNodeNum computation error." << endl;
        exit(1);
    }
}

template<typename Mesh>
void LegacyNodalGeometry::findStencil(Mesh& mesh, InterfaceNode& node)
{
    for (int i = 0; i < 25; i++)
    {
        node.cellStencil[i] = -1;
        node.cellStencilCoeff[i] = 0.0;
    }

    Cell* cells[5][5];
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            cells[i][j] = nullptr;
        }
    }

    const int centerCellLabel = mesh.findCellLabel(node.coord());
    if (centerCellLabel == -1)
    {
        return;
    }

    cells[2][2] = &mesh.cells(centerCellLabel);

    if (cells[2][2]->neiCell[0][0] != -1)
    {
        cells[2][1] = &mesh.cells(cells[2][2]->neiCell[0][0]);

        if (cells[2][1]->neiCell[0][0] != -1)
        {
            cells[2][0] = &mesh.cells(cells[2][1]->neiCell[0][0]);
        }
        if (cells[2][1]->neiCell[1][0] != -1)
        {
            cells[1][1] = &mesh.cells(cells[2][1]->neiCell[1][0]);

            if (cells[1][1]->neiCell[0][0] != -1)
            {
                cells[1][0] = &mesh.cells(cells[1][1]->neiCell[0][0]);

                if (cells[1][0]->neiCell[1][0] != -1)
                {
                    cells[0][0] = &mesh.cells(cells[1][0]->neiCell[1][0]);
                }
            }
            if (cells[1][1]->neiCell[1][0] != -1)
            {
                cells[0][1] = &mesh.cells(cells[1][1]->neiCell[1][0]);
            }
        }
        if (cells[2][1]->neiCell[3][0] != -1)
        {
            cells[3][1] = &mesh.cells(cells[2][1]->neiCell[3][0]);

            if (cells[3][1]->neiCell[0][0] != -1)
            {
                cells[3][0] = &mesh.cells(cells[3][1]->neiCell[0][0]);

                if (cells[3][0]->neiCell[3][0] != -1)
                {
                    cells[4][0] = &mesh.cells(cells[3][0]->neiCell[3][0]);
                }
            }
            if (cells[3][1]->neiCell[3][0] != -1)
            {
                cells[4][1] = &mesh.cells(cells[3][1]->neiCell[3][0]);
            }
        }
    }

    if (cells[2][2]->neiCell[1][0] != -1)
    {
        cells[1][2] = &mesh.cells(cells[2][2]->neiCell[1][0]);

        if (cells[1][2]->neiCell[1][0] != -1)
        {
            cells[0][2] = &mesh.cells(cells[1][2]->neiCell[1][0]);
        }
    }

    if (cells[2][2]->neiCell[2][0] != -1)
    {
        cells[2][3] = &mesh.cells(cells[2][2]->neiCell[2][0]);

        if (cells[2][3]->neiCell[2][0] != -1)
        {
            cells[2][4] = &mesh.cells(cells[2][3]->neiCell[2][0]);
        }
        if (cells[2][3]->neiCell[1][0] != -1)
        {
            cells[1][3] = &mesh.cells(cells[2][3]->neiCell[1][0]);

            if (cells[1][3]->neiCell[2][0] != -1)
            {
                cells[1][4] = &mesh.cells(cells[1][3]->neiCell[2][0]);

                if (cells[1][4]->neiCell[1][0] != -1)
                {
                    cells[0][4] = &mesh.cells(cells[1][4]->neiCell[1][0]);
                }
            }
            if (cells[1][3]->neiCell[1][0] != -1)
            {
                cells[0][3] = &mesh.cells(cells[1][3]->neiCell[1][0]);
            }
        }
        if (cells[2][3]->neiCell[3][0] != -1)
        {
            cells[3][3] = &mesh.cells(cells[2][3]->neiCell[3][0]);

            if (cells[3][3]->neiCell[2][0] != -1)
            {
                cells[3][4] = &mesh.cells(cells[3][3]->neiCell[2][0]);

                if (cells[3][4]->neiCell[3][0] != -1)
                {
                    cells[4][4] = &mesh.cells(cells[3][4]->neiCell[3][0]);
                }
            }
            if (cells[3][3]->neiCell[3][0] != -1)
            {
                cells[4][3] = &mesh.cells(cells[3][3]->neiCell[3][0]);
            }
        }
    }

    if (cells[2][2]->neiCell[3][0] != -1)
    {
        cells[3][2] = &mesh.cells(cells[2][2]->neiCell[3][0]);

        if (cells[3][2]->neiCell[3][0] != -1)
        {
            cells[4][2] = &mesh.cells(cells[3][2]->neiCell[3][0]);
        }
    }

    // 合并检查
    bool virt[5][5];
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            if (cells[i][j] == nullptr)
            {
                virt[i][j] = true;
            }
            else
            {
                virt[i][j] = false;
            }
        }
    }
    // // west
    if (cells[0][1] == nullptr && cells[4][1] == nullptr)// 第1列
    {
        cells[0][1] = cells[0][2];
        cells[1][1] = cells[1][2];
        cells[2][1] = cells[2][2];
        cells[3][1] = cells[3][2];
        cells[4][1] = cells[4][2];

        cells[0][0] = cells[0][3];
        cells[1][0] = cells[1][3];
        cells[2][0] = cells[2][3];
        cells[3][0] = cells[3][3];
        cells[4][0] = cells[4][3];
    }
    else if (cells[0][0] == nullptr && cells[4][0] == nullptr)// 第0列
    {
        cells[0][0] = cells[0][1];
        cells[1][0] = cells[1][1];
        cells[2][0] = cells[2][1];
        cells[3][0] = cells[3][1];
        cells[4][0] = cells[4][1];
    }
    // // south
    if (cells[1][0] == nullptr && cells[1][4] == nullptr)// 第1行
    {
        cells[1][0] = cells[2][0];
        cells[1][1] = cells[2][1];
        cells[1][2] = cells[2][2];
        cells[1][3] = cells[2][3];
        cells[1][4] = cells[2][4];

        cells[0][0] = cells[3][0];
        cells[0][1] = cells[3][1];
        cells[0][2] = cells[3][2];
        cells[0][3] = cells[3][3];
        cells[0][4] = cells[3][4];
    }
    else if (cells[0][0] == nullptr && cells[0][4] == nullptr)// 第0行
    {
        cells[0][0] = cells[1][0];
        cells[0][1] = cells[1][1];
        cells[0][2] = cells[1][2];
        cells[0][3] = cells[1][3];
        cells[0][4] = cells[1][4];
    }
    // // east
    if (cells[0][3] == nullptr && cells[4][3] == nullptr)// 第3列
    {
        cells[0][3] = cells[0][2];
        cells[1][3] = cells[1][2];
        cells[2][3] = cells[2][2];
        cells[3][3] = cells[3][2];
        cells[4][3] = cells[4][2];

        cells[0][4] = cells[0][1];
        cells[1][4] = cells[1][1];
        cells[2][4] = cells[2][1];
        cells[3][4] = cells[3][1];
        cells[4][4] = cells[4][1];
    }
    else if (cells[0][4] == nullptr && cells[4][4] == nullptr)// 第4列
    {
        cells[0][4] = cells[0][3];
        cells[1][4] = cells[1][3];
        cells[2][4] = cells[2][3];
        cells[3][4] = cells[3][3];
        cells[4][4] = cells[4][3];
    }
    // // north
    if (cells[3][0] == nullptr && cells[3][4] == nullptr)// 第3行
    {
        cells[3][0] = cells[2][0];
        cells[3][1] = cells[2][1];
        cells[3][2] = cells[2][2];
        cells[3][3] = cells[2][3];
        cells[3][4] = cells[2][4];

        cells[4][0] = cells[1][0];
        cells[4][1] = cells[1][1];
        cells[4][2] = cells[1][2];
        cells[4][3] = cells[1][3];
        cells[4][4] = cells[1][4];
    }
    else if (cells[4][0] == nullptr && cells[4][4] == nullptr)// 第4行
    {
        cells[4][0] = cells[3][0];
        cells[4][1] = cells[3][1];
        cells[4][2] = cells[3][2];
        cells[4][3] = cells[3][3];
        cells[4][4] = cells[3][4];
    }
    // 角点
    // // 左下角
    if (cells[1][1] == nullptr)
    {
        cells[1][1] = cells[1][2];
        cells[0][1] = cells[0][2];
        cells[1][0] = cells[1][3];
        cells[0][0] = cells[0][3];
    }
    else if (cells[0][0] == nullptr)
    {
        cells[0][0] = cells[1][0];
    }
    // // 右下角
    if (cells[1][3] == nullptr)
    {
        cells[1][3] = cells[1][2];
        cells[0][3] = cells[0][2];
        cells[1][4] = cells[1][1];
        cells[0][4] = cells[0][1];
    }
    else if (cells[0][4] == nullptr)
    {
        cells[0][4] = cells[1][4];
    }
    // // 右上角
    if (cells[3][3] == nullptr)
    {
        cells[3][3] = cells[3][2];
        cells[4][3] = cells[4][2];
        cells[3][4] = cells[3][1];
        cells[4][4] = cells[4][1];
    }
    else if (cells[4][4] == nullptr)
    {
        cells[4][4] = cells[3][4];
    }
    // // 左上角
    if (cells[3][1] == nullptr)
    {
        cells[3][1] = cells[3][2];
        cells[4][1] = cells[4][2];
        cells[3][0] = cells[3][3];
        cells[4][0] = cells[4][3];
    }
    else if (cells[4][0] == nullptr)
    {
        cells[4][0] = cells[3][0];
    }

    // 填充坐标
    int n = 0;
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            double x = cells[2][2]->centerCoord().x + static_cast<double>(j - 2) * cells[2][2]->size().x;
            double y = cells[2][2]->centerCoord().y + static_cast<double>(i - 2) * cells[2][2]->size().y;
            Vector2D coord(x, y);

            node.cellStencil[n] = cells[i][j]->index();
            node.cellStencilCoeff[n] = Interpolation::diracInterpolation(node.coord(), coord, cells[2][2]->size());
            // node.cellStencilCoeff[n] = Interpolation::diracInterpolation(node.coord(), coord, cells[2][2]->size()) * (virt[i][j] ? -1.0 : 1.0);
            n++;
        }
    }
}

template<typename Mesh>
void LegacyNodalGeometry::buildStencil(Mesh& mesh, const int levelSetI)
{
    for (int i = 0; i < nodeIndexToLabel[levelSetI].size(); i++)
    {
        InterfaceNode& node = allNodes_[levelSetI][nodeIndexToLabel[levelSetI][i]];
        findStencil(mesh, node);
    }
}

template<typename Mesh>
void LegacyNodalGeometry::buildNoSlipStencil(Mesh& mesh)
{
    for (int i = 0; i < noSlipNodeNum; i++)
    {
        InterfaceNode& node = noSlipNodes_[i];
        findStencil(mesh, node);
    }
}


