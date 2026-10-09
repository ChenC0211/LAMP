#include "Modules/Plasma/Physics/PicParticleSources.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Base/Interpolation.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
double probeOffset(const QuadtreeMesh& mesh, const Vector2D& p, double length)
{
    double eps = 1e-6 * std::max(1.0, length);
    const int label = mesh.findCellLabel(p);
    if (label >= 0) {
        const Cell& cell = mesh.cells(label);
        const double h = std::min(cell.size().x, cell.size().y);
        if (h > 0.0) eps = std::max(eps, 0.05 * h);
    }
    return eps;
}
}

PicParticleSources::PicParticleSources(Mode mode, const std::vector<SpeciesConfig>& species,
    const std::unordered_map<std::string,std::size_t>& speciesIndexByName,
    const DictionaryGeometry* input, const QuadtreeMesh& mesh,
    const PicEmbeddedGeometry& embeddedGeometry, const Geometry* iGeometry, bool useEmbedded)
    : mode_(mode), x0(mesh.minCoord.x), x1(mesh.maxCoord.x),
      y0(mesh.minCoord.y), y1(mesh.maxCoord.y),
      mx(1e-12 + 1e-6 * (x1 - x0)), my(1e-12 + 1e-6 * (y1 - y0))
{
    if (!input) return;
    const auto& dictGeometry = *input;
    const char* prefix = mode_ == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]";
    auto pickNormal = [&](const Vector2D& p0, const Vector2D& p1, int lsI) -> Vector2D {
        const Vector2D t(p1.x - p0.x, p1.y - p0.y);
        const double len = std::sqrt(t.x * t.x + t.y * t.y);
        if (!(len > 0.0)) return Vector2D(1.0, 0.0);
        const Vector2D n1(-t.y / len, t.x / len), n2(-n1.x, -n1.y);
        if (mode_ == Mode::Explicit && !embeddedGeometry.ready()) return n1;
        const Vector2D mid(0.5 * (p0.x + p1.x), 0.5 * (p0.y + p1.y));
        const double eps = probeOffset(mesh, mid, len);
        const Vector2D pPlus(mid.x + n1.x * eps, mid.y + n1.y * eps);
        const Vector2D pMinus(mid.x + n2.x * eps, mid.y + n2.y * eps);
        if (mode_ == Mode::Explicit)
        {
            const bool in1 = embeddedGeometry.containsPoint(mesh, iGeometry, pPlus, -1);
            const bool in2 = embeddedGeometry.containsPoint(mesh, iGeometry, pMinus, -1);
            return (in1 && !in2) ? n2 : n1;
        }
        if (lsI >= 0)
        {
            auto plasmaSide = [&](const Vector2D& p) {
                const int cell = mesh.findCellLabel(p);
                if (cell < 0) return false;
                const double psi = Interpolation::biQuadInterpolation(mesh,
                    iGeometry->levelSets[static_cast<std::size_t>(lsI)], mesh.cells(cell), p);
                // Implicit sources historically select psi>0 regardless of the
                // embedded geometry's inside-sign policy.
                return std::isfinite(psi) && psi > 0.0;
            };
            const bool in1 = plasmaSide(pPlus), in2 = plasmaSide(pMinus);
            if (in1 && !in2) return n1;
            if (!in1 && in2) return n2;
        }
        auto insideDomain = [&](const Vector2D& p) {
            return p.x > mesh.minCoord.x + mx && p.x < mesh.maxCoord.x - mx &&
                   p.y > mesh.minCoord.y + my && p.y < mesh.maxCoord.y - my;
        };
        const bool in1 = insideDomain(pPlus), in2 = insideDomain(pMinus);
        return (!in1 && in2) ? n2 : n1;
    };
    for (const auto& sp_cfg : species)
    {
        if (!sp_cfg.source.enabled)
            continue;
        if (sp_cfg.source.boundaryId <= 0)
        {
            std::cerr << prefix << "[Source] WARNING: species '" << sp_cfg.name
                      << "' missing sourceBoundaryId; skipping source.\n";
            continue;
        }
        if (!(sp_cfg.source.numberDensity_m3 > 0.0) && !(sp_cfg.source.current_A > 0.0))
        {
            std::cerr << prefix << "[Source] WARNING: species '" << sp_cfg.name
                      << "' has no source density/current; skipping source.\n";
            continue;
        }

        auto sit = speciesIndexByName.find(sp_cfg.name);
        if (sit == speciesIndexByName.end())
        {
            std::cerr << prefix << "[Source] WARNING: species '" << sp_cfg.name
                      << "' not found; skipping source.\n";
            continue;
        }

        const auto bit = dictGeometry.boundaryById.find(sp_cfg.source.boundaryId);
        if (bit == dictGeometry.boundaryById.end())
        {
            std::cerr << prefix << "[Source] WARNING: boundaryId " << sp_cfg.source.boundaryId
                      << " not found in geometry; skipping source for '" << sp_cfg.name << "'.\n";
            continue;
        }

        Emitter src;
        src.species_name = sp_cfg.name;
        src.boundary_id = sp_cfg.source.boundaryId;
        src.number_density_m3 = sp_cfg.source.numberDensity_m3;
        src.temperature_eV = (sp_cfg.source.temperature_eV >= 0.0) ? sp_cfg.source.temperature_eV : sp_cfg.Te_eV;
        src.drift_mps = sp_cfg.source.driftSet ? sp_cfg.source.drift : sp_cfg.drift;
        src.current_A = sp_cfg.source.current_A;
        src.every = (sp_cfg.source.every > 0) ? sp_cfg.source.every : 1;
        src.max_emit_macro = std::max(0, sp_cfg.source.maxEmitMacro);
        src.weight = sp_cfg.weight;
        src.mass_kg = sp_cfg.m;
        src.charge_C = sp_cfg.q;

        int lsForNormal = -1;
        if (mode_ == Mode::Implicit && useEmbedded && embeddedGeometry.ready() && iGeometry &&
            bit->second.kind == GeometricBoundary::fixed)
        {
            lsForNormal = bit->second.levelSetIndex;
            if (lsForNormal < 0 || lsForNormal >= static_cast<int>(iGeometry->levelSets.size()))
                lsForNormal = -1;
        }
        for (int lid : bit->second.lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(dictGeometry.lines.size()))
                continue;
            const Line& ln = dictGeometry.lines[static_cast<std::size_t>(lid)];
            if (ln.ignored())
                continue;
            const Vector2D p0 = ln.points[0];
            const Vector2D p1 = ln.points[1];
            const double len = std::sqrt((p1.x - p0.x) * (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y));
            if (!(len > 0.0))
                continue;
            Segment seg;
            seg.p0 = p0;
            seg.p1 = p1;
            seg.length = len;
            seg.normal = pickNormal(p0, p1, lsForNormal);
            src.segments.push_back(seg);
            src.total_length += len;
        }

        if (src.total_length <= 0.0)
        {
            std::cerr << prefix << "[Source] WARNING: boundaryId " << src.boundary_id
                      << " has no valid segments; skipping source for '" << src.species_name << "'.\n";
            continue;
        }

        emitters_.push_back(std::move(src));
    }
}

void PicParticleSources::emit(const Step& step, const QuadtreeMesh& mesh,
    const PicEmbeddedGeometry& embeddedGeometry, const Geometry* iGeometry,
    std::mt19937_64& rng, std::uniform_real_distribution<double>& unif01,
    ParticleSpawnBuffer& sourceSpawnPending)
{
    if (mode_ == Mode::Explicit) {
        if (sourceSpawnPending.empty()) pendingBoundaryIds_.clear();
        else requireExplicit(sourceSpawnPending);
    }
    const double dt_si_src = step.dtSi;
    const double v0_mps_src = step.velocityScale;
    const double slab_m = step.slabThickness > 0.0 ? step.slabThickness : 1.0;
    constexpr double kPi = 3.141592653589793238462643383279502884;
    for (const auto& src : emitters_)
    {
        if (src.every > 1 && ((step.index + 1) % src.every) != 0)
            continue;
        if (!(src.total_length > 0.0))
            continue;
        if (!(src.weight > 0.0))
            continue;

        double nphys = 0.0;
        const double abs_q = std::abs(src.charge_C);
        if (src.current_A > 0.0 && abs_q > 0.0)
        {
            nphys = src.current_A * dt_si_src / abs_q;
        }
        else
        {
            const double kT_J = std::max(0.0, src.temperature_eV) * PhysicalConstants::eV_J;
            const double vth_flux = (src.mass_kg > 0.0 && kT_J > 0.0)
                                        ? std::sqrt(kT_J / (2.0 * kPi * src.mass_kg))
                                        : 0.0;
            double v_n = 0.0;
            for (const auto& seg : src.segments)
            {
                const double vn_seg = src.drift_mps.x * seg.normal.x + src.drift_mps.y * seg.normal.y;
                v_n += std::max(0.0, vn_seg) * seg.length;
            }
            if (src.total_length > 0.0)
                v_n /= src.total_length;
            const double v_eff = v_n + vth_flux;
            nphys = src.number_density_m3 * (src.total_length * slab_m) * v_eff * dt_si_src;
        }

        if (!(nphys > 0.0))
            continue;

        const double n_macro = nphys / src.weight;
        int n_emit = static_cast<int>(std::floor(n_macro));
        const double frac = n_macro - static_cast<double>(n_emit);
        if (unif01(rng) < frac)
            n_emit += 1;
        if (src.max_emit_macro > 0)
            n_emit = std::min(n_emit, src.max_emit_macro);
        if (n_emit <= 0)
            continue;

        const double kT_J = std::max(0.0, src.temperature_eV) * PhysicalConstants::eV_J;
        const double vStd_code = (src.mass_kg > 0.0 && kT_J > 0.0 && v0_mps_src > 0.0)
                                     ? std::sqrt(kT_J / src.mass_kg) / v0_mps_src
                                     : 0.0;
        std::normal_distribution<double> gauss(0.0, vStd_code);
        const Vector2D drift_code(src.drift_mps.x / v0_mps_src, src.drift_mps.y / v0_mps_src);

        for (int i = 0; i < n_emit; ++i)
        {
            const double r = unif01(rng) * src.total_length;
            const Segment* seg = nullptr;
            double acc = 0.0;
            for (const auto& s : src.segments)
            {
                acc += s.length;
                if (r <= acc)
                {
                    seg = &s;
                    break;
                }
            }
            if (!seg)
                seg = &src.segments.back();

            const double t = unif01(rng);
            Vector2D pos(seg->p0.x + t * (seg->p1.x - seg->p0.x),
                         seg->p0.y + t * (seg->p1.y - seg->p0.y));
            double eps = 1e-6 * std::max(1.0, seg->length);
            const int midCell = mesh.findCellLabel(pos);
            if (midCell >= 0)
            {
                const Cell& c = mesh.cells(midCell);
                const double h = std::min(c.size().x, c.size().y);
                if (h > 0.0)
                    eps = std::max(eps, 0.05 * h);
            }
            Vector2D n = seg->normal;
            Vector2D pos_out(pos.x + n.x * eps, pos.y + n.y * eps);
            if (mode_ == Mode::Explicit && embeddedGeometry.ready() && embeddedGeometry.containsPoint(mesh, iGeometry, pos_out, -1))
                pos_out = Vector2D(pos.x - n.x * eps, pos.y - n.y * eps);
            if (mode_ == Mode::Explicit && embeddedGeometry.ready() && embeddedGeometry.containsPoint(mesh, iGeometry, pos_out, -1))
                continue;
            pos_out.x = std::min(std::max(pos_out.x, x0 + mx), x1 - mx);
            pos_out.y = std::min(std::max(pos_out.y, y0 + my), y1 - my);

            Vector2D v = drift_code;
            if (vStd_code > 0.0)
            {
                v.x += gauss(rng);
                v.y += gauss(rng);
            }
            const double vdot = v.x * n.x + v.y * n.y;
            if (vdot < 0.0)
            {
                v.x -= 2.0 * vdot * n.x;
                v.y -= 2.0 * vdot * n.y;
            }
            const double vz = (vStd_code > 0.0) ? gauss(rng) : 0.0;

            Particle np;
            np.position = pos_out;
            np.velocity = v;
            np.vz = vz;
            np.w = src.weight;
            int inCellL = -1;
            if (midCell >= 0)
                inCellL = mesh.findCellLabelFromHint(midCell, np.position);
            if (inCellL < 0)
                inCellL = mesh.findCellLabel(np.position);
            if (inCellL < 0)
                continue;
            np.inCellL = inCellL;
            np.inBackgroundCellL = mesh.findBackgroundCellLabel(np.position);

            sourceSpawnPending.Push(src.species_name, np);
            if (mode_ == Mode::Explicit) pendingBoundaryIds_.push_back(src.boundary_id);
        }
    }
}

void PicParticleSources::requireExplicit(const ParticleSpawnBuffer& pending) const
{
    if (mode_ != Mode::Explicit || pending.size() != pendingBoundaryIds_.size())
        throw std::logic_error("Explicit source accounting requires retained boundary provenance");
}

PicSpawnManagement::Stats PicParticleSources::thinExplicit(const PicSpawnManagement& manager,
    ParticleSpawnBuffer& pending, bool enabled, double probability,
    std::mt19937_64& rng, std::uniform_real_distribution<double>& uniform)
{
    requireExplicit(pending);
    return manager.thin(pending,enabled,probability,rng,uniform,&pendingBoundaryIds_);
}

PicParticleSources::Report PicParticleSources::accountRetained(const ParticleSpawnBuffer& pending,
    const std::vector<ParticleType>& species, const std::vector<double>& weights,
    const std::unordered_map<std::string,std::size_t>& speciesIndex,
    double velocityScale, const std::array<int,4>& domainBoundaryIds) const
{
    requireExplicit(pending);
    Report report;
    for (std::size_t i = 0; i < pending.size(); ++i) {
        const auto& item = pending.items()[i];
        const auto found = speciesIndex.find(item.species_name);
        if (found == speciesIndex.end()) continue;
        const auto index = found->second;
        const auto& sp = species[index];
        const auto& p = item.particle;
        const double defaultWeight = index < weights.size() ? weights[index] : 1.0;
        const double weight = p.w > 0.0 ? p.w : defaultWeight;
        double energy = 0.;
        if (sp.mass() > 0.0) {
            const double vx = p.velocity.x * velocityScale, vy = p.velocity.y * velocityScale, vz = p.vz * velocityScale;
            const double v2 = vx * vx + vy * vy + vz * vz;
            energy = 0.5 * sp.mass() * v2;
        }
        report.energyJ += energy * weight;
        for (std::size_t side = 0; side < domainBoundaryIds.size(); ++side) {
            if (domainBoundaryIds[side] > 0 && pendingBoundaryIds_[i] == domainBoundaryIds[side]) {
                report.domainChargeC[side] += sp.charge() * weight;
                break;
            }
        }
    }
    return report;
}
