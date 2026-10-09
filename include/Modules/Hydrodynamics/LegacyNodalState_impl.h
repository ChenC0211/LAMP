#pragma once

#include "Modules/Hydrodynamics/LegacyNodalState.h"

using namespace std;

template<typename Mesh>
void LegacyNodalState::initializeGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry)
{
    cout << endl << "[SharpInterface] Initializing interface nodes..." << endl;
    ensureLevelSetCount(static_cast<int>(iGeometry.levelSets.size()));
    geometry_.initializeGeometry(mesh, iGeometry);
    for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); ++levelSetI)
        evaporationVelocity_[levelSetI].assign(geometry_.allNodes(levelSetI).size(), 0.0);
    noSlipFx.initialize(geometry_.noSlipNodeNum, 0.0);
    noSlipFy.initialize(geometry_.noSlipNodeNum, 0.0);
    noSlipUx.initialize(geometry_.noSlipNodeNum, 0.0);
    noSlipUy.initialize(geometry_.noSlipNodeNum, 0.0);
    cout << "[SharpInterface] Interface nodes initialized." << endl;
}

template<typename Mesh>
void LegacyNodalState::updateGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, 
    const VertField<double>& ux, const VertField<double>& uy, const VertField<double>& preUx,
     const VertField<double>& preUy, const double& dt, const double& preDt, const int loop)
{
    cout << endl << "[SharpInterface] Updating interface nodes..." << endl;

    ensureLevelSetCount(static_cast<int>(iGeometry.levelSets.size()));
    for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); levelSetI++)
    {
        prePreEvaporationVelocity_[levelSetI] = preEvaporationVelocity_[levelSetI];
        preEvaporationVelocity_[levelSetI] = evaporationVelocity_[levelSetI];
        geometry_.updateGeometry(mesh, iGeometry, levelSetI);
        if (iGeometry.levelSets[levelSetI].moving())
        {
            evaporationVelocity_[levelSetI].assign(geometry_.allNodes(levelSetI).size(), 0.0);
            computeVelocity(mesh, levelSetI, ux, uy, preUx, preUy, dt, preDt, loop);
        }
    }
    
    cout << "[SharpInterface] Interface nodes updated." << endl;
}

template<typename Mesh>
void LegacyNodalState::computeVelocity(Mesh& mesh, const int levelSetI, const VertField<double>& ux, const VertField<double>& uy, const VertField<double>& preUx, const VertField<double>& preUy, const double& dt, const double& preDt, const int loop)
{
    Uxs[levelSetI].initialize(geometry_.nodeIndexToLabel[levelSetI].size(), 0.0);
    Uys[levelSetI].initialize(geometry_.nodeIndexToLabel[levelSetI].size(), 0.0);
    if (loop > 1)
    {
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < geometry_.nodeIndexToLabel[levelSetI].size(); i++)
        {
            InterfaceNode& node = geometry_.allNodes(levelSetI)[geometry_.nodeIndexToLabel[levelSetI][i]];
            Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));

            Vector2D vel = Vector2D(Interpolation::biQuadInterpolation(mesh, ux, cell, node.coord()), Interpolation::biQuadInterpolation(mesh, uy, cell, node.coord()));

            Vector2D halfCoord = node.coord() - vel * 0.5 * dt;
            Cell& halfCell = mesh.cells(mesh.findCellLabel(halfCoord));

            Vector2D preVel(Interpolation::biQuadInterpolation(mesh, ux, halfCell, halfCoord), Interpolation::biQuadInterpolation(mesh, uy, halfCell, halfCoord));
            Vector2D prePreVel(Interpolation::biQuadInterpolation(mesh, preUx, halfCell, halfCoord), Interpolation::biQuadInterpolation(mesh, preUy, halfCell, halfCoord));

            Vector2D halfVel = preVel * ((2.0 * preDt + dt) / (2.0 * preDt)) - prePreVel * (dt / (2.0 * preDt));

            Uxs[levelSetI].value[i] = halfVel.x;
            Uys[levelSetI].value[i] = halfVel.y;
        }
    }
}

template<typename Field>
double LegacyNodalState::interfaceInterpolation(const Vector2D& coord, Field& field, const int levelSetI, const bool isTimeDiff, const bool isVt, const Time& timeStep, const Vector2D& normal)
{
    double eta[3] = { 0.0 };
	double value[3] = { 0.0 };

    if (timeStep == Time::current)
    {
        int nearestNodeLabel = geometry_.nodeIndexToLabel[levelSetI][0];
        double minDist2 = (geometry_.allNodes(levelSetI)[geometry_.nodeIndexToLabel[levelSetI][0]].coord() - coord).mod2();
        for (int i = 1; i < geometry_.nodeIndexToLabel[levelSetI].size(); i++)
        {
            double dist2 = (geometry_.allNodes(levelSetI)[geometry_.nodeIndexToLabel[levelSetI][i]].coord() - coord).mod2();
            if (dist2 < minDist2)
            {
                minDist2 = dist2;
                nearestNodeLabel = geometry_.nodeIndexToLabel[levelSetI][i];
            }
        }

        if (isVt)
        {
            InterfaceNode& node = geometry_.allNodes(levelSetI)[nearestNodeLabel];
            for (int i = 0; i < 3; i++)
            {
                InterfaceNode& relatedNode = geometry_.allNodes(levelSetI)[node.relatedNodes[i]];

                eta[i] = -(relatedNode.coord().x - coord.x) * normal.y + (relatedNode.coord().y - coord.y) * normal.x;
                value[i] = -Uxs[levelSetI].value[relatedNode.index()] * normal.y + Uys[levelSetI].value[relatedNode.index()] * normal.x;
            }
        }
        else
        {
            InterfaceNode& node = geometry_.allNodes(levelSetI)[nearestNodeLabel];
            for (int i = 0; i < 3; i++)
            {
                InterfaceNode& relatedNode = geometry_.allNodes(levelSetI)[node.relatedNodes[i]];

                eta[i] = -(relatedNode.coord().x - coord.x) * normal.y + (relatedNode.coord().y - coord.y) * normal.x;
                value[i] = field.value[relatedNode.index()];
            }
        }
    }
    else if (timeStep == Time::previous)
    {
        int nearestNodeLabel = geometry_.preNodeIndexToLabel[levelSetI][0];
        double minDist2 = (geometry_.allPreNodes(levelSetI)[geometry_.preNodeIndexToLabel[levelSetI][0]].coord() - coord).mod2();
        for (int i = 1; i < geometry_.preNodeIndexToLabel[levelSetI].size(); i++)
        {
            double dist2 = (geometry_.allPreNodes(levelSetI)[geometry_.preNodeIndexToLabel[levelSetI][i]].coord() - coord).mod2();
            if (dist2 < minDist2)
            {
                minDist2 = dist2;
                nearestNodeLabel = geometry_.preNodeIndexToLabel[levelSetI][i];
            }
        }

        if (isVt)
        {
            InterfaceNode& node = geometry_.allPreNodes(levelSetI)[nearestNodeLabel];
            for (int i = 0; i < 3; i++)
            {
                InterfaceNode& relatedNode = geometry_.allPreNodes(levelSetI)[node.relatedNodes[i]];

                eta[i] = -(relatedNode.coord().x - coord.x) * normal.y + (relatedNode.coord().y - coord.y) * normal.x;
                value[i] = -preUxs[levelSetI].value[relatedNode.index()] * normal.y + preUys[levelSetI].value[relatedNode.index()] * normal.x;
            }
        }
        else
        {
            InterfaceNode& node = geometry_.allPreNodes(levelSetI)[nearestNodeLabel];
            for (int i = 0; i < 3; i++)
            {
                InterfaceNode& relatedNode = geometry_.allPreNodes(levelSetI)[node.relatedNodes[i]];

                eta[i] = -(relatedNode.coord().x - coord.x) * normal.y + (relatedNode.coord().y - coord.y) * normal.x;
                value[i] = field.value[relatedNode.index()];
            }
        }
    }
    else
    {
        int nearestNodeLabel = geometry_.prePreNodeIndexToLabel[levelSetI][0];
        double minDist2 = (geometry_.allPrePreNodes(levelSetI)[geometry_.prePreNodeIndexToLabel[levelSetI][0]].coord() - coord).mod2();
        for (int i = 1; i < geometry_.prePreNodeIndexToLabel[levelSetI].size(); i++)
        {
            double dist2 = (geometry_.allPrePreNodes(levelSetI)[geometry_.prePreNodeIndexToLabel[levelSetI][i]].coord() - coord).mod2();
            if (dist2 < minDist2)
            {
                minDist2 = dist2;
                nearestNodeLabel = geometry_.prePreNodeIndexToLabel[levelSetI][i];
            }
        }

        if (isVt)
        {
            InterfaceNode& node = geometry_.allPrePreNodes(levelSetI)[nearestNodeLabel];
            for (int i = 0; i < 3; i++)
            {
                InterfaceNode& relatedNode = geometry_.allPrePreNodes(levelSetI)[node.relatedNodes[i]];

                eta[i] = -(relatedNode.coord().x - coord.x) * normal.y + (relatedNode.coord().y - coord.y) * normal.x;
                value[i] = -prePreUxs[levelSetI].value[relatedNode.index()] * normal.y + prePreUys[levelSetI].value[relatedNode.index()] * normal.x;
            }
        }
        else
        {
            InterfaceNode& node = geometry_.allPrePreNodes(levelSetI)[nearestNodeLabel];
            for (int i = 0; i < 3; i++)
            {
                InterfaceNode& relatedNode = geometry_.allPrePreNodes(levelSetI)[node.relatedNodes[i]];

                eta[i] = -(relatedNode.coord().x - coord.x) * normal.y + (relatedNode.coord().y - coord.y) * normal.x;
                value[i] = field.value[relatedNode.index()];
            }
        }
    }

    if (isTimeDiff)
    {
        return -(eta[1] + eta[2]) / ((eta[0] - eta[1]) * (eta[0] - eta[2])) * value[0]
            - (eta[0] + eta[2]) / ((eta[1] - eta[0]) * (eta[1] - eta[2])) * value[1]
            - (eta[0] + eta[1]) / ((eta[2] - eta[0]) * (eta[2] - eta[1])) * value[2];
    }
    else
    {
        return (eta[1] * eta[2] / ((eta[0] - eta[1]) * (eta[0] - eta[2]))) * value[0]
            + (eta[0] * eta[2] / ((eta[1] - eta[0]) * (eta[1] - eta[2]))) * value[1]
            + (eta[0] * eta[1] / ((eta[2] - eta[0]) * (eta[2] - eta[1]))) * value[2];
    }
}
