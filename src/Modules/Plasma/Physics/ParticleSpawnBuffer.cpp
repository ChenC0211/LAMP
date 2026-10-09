#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"

#include <iomanip>
#include <ostream>
#include <sstream>

void ParticleSpawnBuffer::Push(const std::string& species_name, const Particle& p)
{
    SpawnedParticle sp;
    sp.species_name = species_name;
    sp.particle = p;
    items_.push_back(std::move(sp));
}

void ParticleSpawnBuffer::PrintSummary(std::ostream& os, std::size_t max_items) const
{
    os << "ParticleSpawnBuffer: n=" << items_.size() << "\n";
    const std::size_t nshow = std::min(max_items, items_.size());
    os << std::scientific << std::setprecision(6);
    for (std::size_t i = 0; i < nshow; ++i)
    {
        const auto& it = items_[i];
        os << "  [" << i << "] " << it.species_name
           << "  x=" << it.particle.position.x
           << " y=" << it.particle.position.y
           << "  vx=" << it.particle.velocity.x
           << " vy=" << it.particle.velocity.y
           << "  cell=" << it.particle.inCellL
           << "\n";
    }
    if (items_.size() > nshow)
        os << "  ... (" << (items_.size() - nshow) << " more)\n";
}

std::string ParticleSpawnBuffer::SummaryString(std::size_t max_items) const
{
    std::ostringstream ss;
    PrintSummary(ss, max_items);
    return ss.str();
}
