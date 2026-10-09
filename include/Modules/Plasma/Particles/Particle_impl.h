#pragma once

#include "Modules/Plasma/Particles/Particle.h"

using namespace std;

template<typename Mesh>
void Particle::updateWeight(Mesh& mesh, Cell& interpolationCell)
{
    const double hx = interpolationCell.size().x;
    const double hy = interpolationCell.size().y;

	double x = (position.x - mesh.verts(interpolationCell.cellToVert[0]).coord().x) / hx;
	double y = (position.y - mesh.verts(interpolationCell.cellToVert[0]).coord().y) / hy;

	weight[0] = (1.0 - x) * (1.0 - y);
	weight[1] = x * (1.0 - y);
	weight[2] = x * y;
	weight[3] = (1.0 - x) * y;
}

template<typename Mesh, typename Field>
void Particle::getDiscreteVert(Mesh& mesh, Field& LSField, const int levelSetI)
{
    discreteVertL = -1;
    Cell& inCell = mesh.cells(inCellL);
	if(inCell.isCut[levelSetI])
    {
        for(int i = 0; i < 4; i++)
        {
            Vert& vert = mesh.verts(inCell.cellToVert[i]);
            if(LSField.value[vert.index()] > 0.0
                && vert.onCutEdge[levelSetI])
            {
                discreteVertL = vert.label();
                return;
            }
        }

        if(discreteVertL == -1)
        {
            cout << "[Particle] ERROR: failed to find discrete vertex. (Particle.cpp) Line: "
            << __LINE__ << endl;
            exit(-1);
        }
    }
}