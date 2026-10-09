#include "Modules/Plasma/Particles/Particle.h"

using namespace std;

Particle::Particle()
:
position(0.0, 0.0),
velocity(0.0, 0.0)
{
    vz = 0.0;
    inCellL = -1;
    inBackgroundCellL = -1;
    discreteVertL = -1;
    w = 1.0;
    for (int i = 0; i < 4; ++i)
        weight[i] = 0.0;
    emitChance = 0.0;
    levelSetValues.clear();
    rngState = 0;
    trackId = 0;
}

Particle::Particle(double& px, double& py, double& vx, double& vy, int inCellL)
:
position(px, py),
velocity(vx, vy)
{
    this->inCellL = inCellL;
    this->vz = 0.0;
    this->w = 1.0;
    this->rngState = 0;
    this->trackId = 0;
}

Particle::Particle(const Particle& particle)
:
position(particle.position),
velocity(particle.velocity)
{
    this->vz = particle.vz;
    this->inCellL = particle.inCellL;
    this->inBackgroundCellL = particle.inBackgroundCellL;
    this->discreteVertL = particle.discreteVertL;
    this->w = particle.w;

    for(int i = 0; i < 4; i++)
    {
        this->weight[i] = particle.weight[i];
    }

    this->emitChance = particle.emitChance;
    this->levelSetValues = particle.levelSetValues;

    this->rngState = particle.rngState;
    this->trackId = particle.trackId;
}

Particle::~Particle() {}
