#pragma once
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/Physics/PicSpawnManagement.h"
#include <array>
#include <random>
#include <unordered_map>

class ParticleSpawnBuffer;
class ParticleType;

// Owns the prepared boundary segments and species source settings. Geometry is
// sampled at setup, as in the original drivers; emission borrows the current
// mesh/solid view and caller's RNG. Buffer clearing, thinning and merge stay at
// the integrator's existing stages.
class PicParticleSources
{
public:
    enum class Mode { Explicit, Implicit };
    using Geometry = PicEmbeddedGeometry::Geometry;
    struct Step { int index; double dtSi, velocityScale, slabThickness; };
    PicParticleSources(Mode mode, const std::vector<SpeciesConfig>& species,
        const std::unordered_map<std::string,std::size_t>& speciesIndex,
        const DictionaryGeometry* input, const QuadtreeMesh& mesh,
        const PicEmbeddedGeometry& embedded, const Geometry* geometry, bool useEmbedded = true);
    bool empty() const { return emitters_.empty(); }
    // An empty pending buffer starts a new batch. Explicit appends retain
    // boundary provenance; filter them through thinExplicit, not externally.
    void emit(const Step& step, const QuadtreeMesh& mesh,
        const PicEmbeddedGeometry& embedded, const Geometry* geometry,
        std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform,
        ParticleSpawnBuffer& pending);
    PicSpawnManagement::Stats thinExplicit(const PicSpawnManagement& manager,
        ParticleSpawnBuffer& pending, bool enabled, double probability,
        std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform);
    struct Report {
        double energyJ = 0.;
        std::array<double,4> domainChargeC{};
    };
    // Explicit retained-particle accounting, before the existing commit/cull
    // stages. Embedded sources contribute energy but not domain-side charge.
    Report accountRetained(const ParticleSpawnBuffer& pending,
        const std::vector<ParticleType>& species, const std::vector<double>& weights,
        const std::unordered_map<std::string,std::size_t>& speciesIndex,
        double velocityScale, const std::array<int,4>& domainBoundaryIds) const;

private:
    void requireExplicit(const ParticleSpawnBuffer& pending) const;
    struct Segment { Vector2D p0, p1, normal; double length = 0.; };
    struct Emitter {
        std::string species_name;
        int boundary_id = 0;
        double number_density_m3 = 0., temperature_eV = -1.;
        Vector2D drift_mps{0.,0.};
        double current_A = 0.;
        int every = 1, max_emit_macro = 0;
        double weight = 1., mass_kg = 0., charge_C = 0.;
        std::vector<Segment> segments;
        double total_length = 0.;
    };
    const Mode mode_;
    const double x0, x1, y0, y1, mx, my;
    std::vector<Emitter> emitters_;
    std::vector<int> pendingBoundaryIds_;
};
