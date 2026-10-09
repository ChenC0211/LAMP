#include "Modules/Plasma/Physics/PostCollisionModel.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <array>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

#include "Modules/Plasma/Physics/PhysicalConstants.h"

namespace
{
std::string ToLowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

double ClampNonNegative(double x)
{
    return (x < 0.0) ? 0.0 : x;
}

double KineticEnergyEV3D(const Vector2D& v, double vz, double m_kg)
{
    const double v2 = v.mod2() + vz * vz;
    const double E_J = 0.5 * m_kg * v2;
    return E_J / PhysicalConstants::eV_J;
}

struct Velocity3D
{
    Vector2D v;
    double vz = 0.0;
};

Velocity3D MakeVelocity3D(const Vector2D& v, double vz)
{
    Velocity3D out;
    out.v = v;
    out.vz = vz;
    return out;
}

double Dot3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> NormalizeDirection3D(const Velocity3D& v)
{
    const double s2 = v.v.mod2() + v.vz * v.vz;
    const double s = std::sqrt(std::max(0.0, s2));
    if (!(s > 0.0))
        return {1.0, 0.0, 0.0};
    return {v.v.x / s, v.v.y / s, v.vz / s};
}

template <class Rng>
Velocity3D SampleMaxwellianVelocity3D(double m_kg,
                                      double T_K,
                                      const Vector2D& drift,
                                      Rng& rng)
{
    Velocity3D out;
    out.v = drift;
    out.vz = 0.0;
    if (!(m_kg > 0.0) || !(T_K > 0.0) || !std::isfinite(m_kg) || !std::isfinite(T_K))
        return out;

    const double sigma = std::sqrt(PhysicalConstants::kB_JpK * T_K / m_kg);
    std::normal_distribution<double> nd(0.0, 1.0);
    out.v = Vector2D(drift.x + sigma * nd(rng), drift.y + sigma * nd(rng));
    out.vz = sigma * nd(rng);
    return out;
}

double SpeedFromEnergyEV(double E_eV, double m_kg)
{
    const double E_J = std::max(0.0, E_eV) * PhysicalConstants::eV_J;
    // v = sqrt(2E/m)
    return std::sqrt(2.0 * E_J / m_kg);
}

template <class Rng>
Vector2D RandomIsotropicDirection2D(Rng& rng)
{
    constexpr double kPi = 3.141592653589793238462643383279502884;
    std::uniform_real_distribution<double> unif(0.0, 2.0 * kPi);
    const double th = unif(rng);
    return Vector2D(std::cos(th), std::sin(th));
}

template <class Rng>
std::array<double, 3> RandomIsotropicDirection3D(Rng& rng)
{
    constexpr double kPi = 3.141592653589793238462643383279502884;
    std::uniform_real_distribution<double> unif01(0.0, 1.0);
    const double u = unif01(rng);
    const double v = unif01(rng);
    const double mu = 2.0 * u - 1.0; // cos(theta) uniform in [-1,1]
    const double phi = 2.0 * kPi * v;
    const double sinT = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    return {sinT * std::cos(phi), sinT * std::sin(phi), mu};
}

template <class Rng>
bool ApplyIonNeutralElasticScatterCOM(Particle& projectile,
                                      const PostCollisionModel::TargetState* target_state,
                                      double m_proj,
                                      double m_tgt,
                                      const std::string& scatter,
                                      Rng& rng,
                                      std::string* err)
{
    if (!(m_proj > 0.0) || !std::isfinite(m_proj))
    {
        if (err)
            *err = "Invalid projectile mass for COM elastic scatter";
        return false;
    }
    if (!(m_tgt > 0.0) || !std::isfinite(m_tgt))
    {
        if (err)
            *err = "Invalid target mass for COM elastic scatter";
        return false;
    }

    const Velocity3D vp = MakeVelocity3D(projectile.velocity, projectile.vz);
    const Velocity3D vt = target_state ? MakeVelocity3D(target_state->velocity, target_state->vz)
                                       : MakeVelocity3D(Vector2D(0.0, 0.0), 0.0);
    const Velocity3D g = MakeVelocity3D(vp.v - vt.v, vp.vz - vt.vz);
    const double g2 = g.v.mod2() + g.vz * g.vz;
    const double gmag = std::sqrt(std::max(0.0, g2));
    if (!(gmag > 0.0))
        return true;

    const double inv_msum = 1.0 / (m_proj + m_tgt);
    const Velocity3D vcm = MakeVelocity3D((vp.v * m_proj + vt.v * m_tgt) * inv_msum,
                                          (m_proj * vp.vz + m_tgt * vt.vz) * inv_msum);

    std::array<double, 3> dir = RandomIsotropicDirection3D(rng);
    if (scatter == "backscatter" || scatter == "backscat")
    {
        const std::array<double, 3> ghat = NormalizeDirection3D(g);
        if (Dot3(dir, ghat) > 0.0)
        {
            dir[0] = -dir[0];
            dir[1] = -dir[1];
            dir[2] = -dir[2];
        }
    }

    // For elastic two-body scattering, the COM relative-speed magnitude is preserved.
    const double alpha = m_tgt * inv_msum;
    projectile.velocity = Vector2D(vcm.v.x + alpha * gmag * dir[0],
                                   vcm.v.y + alpha * gmag * dir[1]);
    projectile.vz = vcm.vz + alpha * gmag * dir[2];
    return true;
}

// Best-effort channel kind inference to keep older dicts compatible.
// We use the channel name, cross-section key, and product list.
std::string InferKind(const ReactionSet& set, const ReactionChannel& ch)
{
    const std::string name = ToLowerCopy(ch.name);
    const std::string key = ToLowerCopy(ch.cross_section_key);

    auto contains = [](const std::string& s, const char* needle) {
        return s.find(needle) != std::string::npos;
    };

    if (contains(name, "elastic") || contains(key, "elastic"))
        return "elastic";
    if (contains(name, "excitation") || contains(name, "excite") || contains(key, "excitation"))
        return "excitation";
    if (contains(name, "ionization") || contains(name, "ionize") || contains(key, "ionization"))
        return "ionization";
    if (contains(name, "charge_exchange") || contains(name, "cx") || contains(key, "charge_exchange") || contains(key, ":cx"))
        return "charge_exchange";

    // Heuristic: ionization typically produces an electron + an ion.
    bool has_e = false;
    bool has_ion = false;
    for (const auto& p : ch.products)
    {
        const std::string pl = ToLowerCopy(p);
        if (pl == "e" || pl == "e-" || pl == "electron")
            has_e = true;
        if (pl.find('+') != std::string::npos)
            has_ion = true;
    }
    if (has_e && has_ion)
        return "ionization";

    // Heuristic: if projectile looks like an ion and products include a neutral distinct from target.
    if (ToLowerCopy(set.projectile).find('+') != std::string::npos)
    {
        for (const auto& p : ch.products)
        {
            const std::string pl = ToLowerCopy(p);
            if (pl != ToLowerCopy(set.projectile) && pl != ToLowerCopy(set.target))
                return "charge_exchange";
        }
    }

    // Default.
    return "elastic";
}

std::string ChooseScatterModel(const PostCollisionModel::Options& opt, const ReactionChannel& ch)
{
    if (!ch.scatter_model.empty())
        return ToLowerCopy(ch.scatter_model);
    return ToLowerCopy(opt.default_scatter_model);
}

std::string ChooseIonizationShareModel(const PostCollisionModel::Options& opt, const ReactionChannel& ch)
{
    if (!ch.energy_share_model.empty())
        return ToLowerCopy(ch.energy_share_model);
    return ToLowerCopy(opt.default_ionization_share_model);
}

} // namespace

PostCollisionModel::PostCollisionModel(const SpeciesDatabase& species_db, Options opt)
    : species_db_(species_db), opt_(std::move(opt))
{
}

namespace
{
template <class Rng>
bool ApplyImpl(const SpeciesDatabase& species_db,
               const PostCollisionModel::Options& opt,
               const ReactionSet& set, const ReactionChannel& ch,
               Particle& projectile,
               Rng& rng,
               const PostCollisionModel::TargetState* target_state,
               ParticleSpawnBuffer* spawn,
               PostCollisionModel::Outcome* out,
               std::string* err)
{
    const Species* projSp = species_db.Find(set.projectile);
    if (!projSp)
    {
        if (err)
            *err = "Unknown projectile species: " + set.projectile;
        return false;
    }
    const double m = projSp->mass_kg;
    if (!(m > 0.0) || !std::isfinite(m))
    {
        if (err)
            *err = "Invalid projectile mass for species: " + set.projectile;
        return false;
    }

    const double E0 = KineticEnergyEV3D(projectile.velocity, projectile.vz, m);
    const std::string kind = InferKind(set, ch);
    const std::string scatter = ChooseScatterModel(opt, ch);
    const bool has_explicit_scatter = !ch.scatter_model.empty();

    auto setSpeedWithScatter = [&](double E_eV)
    {
        const double vmag = SpeedFromEnergyEV(E_eV, m);

        std::array<double, 3> dir{1.0, 0.0, 0.0};
        if (scatter == "isotropic")
        {
            dir = RandomIsotropicDirection3D(rng);
        }
        else if (scatter == "backscatter" || scatter == "backscat")
        {
            const double v2 = projectile.velocity.mod2() + projectile.vz * projectile.vz;
            const double s = std::sqrt(std::max(0.0, v2));
            if (s > 0.0)
            {
                const std::array<double, 3> vhat = {
                    projectile.velocity.x / s,
                    projectile.velocity.y / s,
                    projectile.vz / s
                };

                dir = RandomIsotropicDirection3D(rng);
                const double dot = dir[0] * vhat[0] + dir[1] * vhat[1] + dir[2] * vhat[2];
                if (dot > 0.0)
                {
                    dir[0] = -dir[0];
                    dir[1] = -dir[1];
                    dir[2] = -dir[2];
                }
            }
            else
            {
                dir = RandomIsotropicDirection3D(rng);
            }
        }
        else
        {
            // Default: keep direction (forward scattering) in 3D velocity space.
            const double v2 = projectile.velocity.mod2() + projectile.vz * projectile.vz;
            const double s = std::sqrt(std::max(0.0, v2));
            if (s > 0.0)
                dir = {projectile.velocity.x / s, projectile.velocity.y / s, projectile.vz / s};
        }

        projectile.velocity = Vector2D(dir[0] * vmag, dir[1] * vmag);
        projectile.vz = dir[2] * vmag;
    };

    PostCollisionModel::Outcome local;
    local.channel_name = ch.name;
    local.kind = kind;
    local.E_before_eV = E0;
    local.E_after_eV = E0;

    if (kind == "elastic")
    {
        // Preserve energy by default, but for ion-neutral ELASTIC we may apply
        // a minimal slow-down model by thermalizing the ion to the background
        // gas temperature (without explicitly tracking fast neutrals).
        const bool targetIsGas = (!opt.gas_species_name.empty())
            && (ToLowerCopy(set.target) == ToLowerCopy(opt.gas_species_name));
        const bool comScatterSupported = projSp->IsIon() && targetIsGas
            && has_explicit_scatter
            && (scatter == "isotropic" || scatter == "backscatter" || scatter == "backscat");
        if (comScatterSupported)
        {
            if (!(opt.gas_particle_mass_kg > 0.0) || !std::isfinite(opt.gas_particle_mass_kg))
            {
                if (err)
                    *err = "Ion-neutral elastic scatter requires valid background gas mass for target '" + set.target + "'";
                return false;
            }
            if (!ApplyIonNeutralElasticScatterCOM(projectile, target_state, m, opt.gas_particle_mass_kg, scatter, rng, err))
                return false;
        }
        else if (opt.ion_neutral_elastic_thermalize && !has_explicit_scatter && projSp->IsIon() && targetIsGas)
        {
            const Velocity3D v = SampleMaxwellianVelocity3D(m, opt.gas_temperature_K, /*drift=*/Vector2D(0.0, 0.0), rng);
            projectile.velocity = v.v;
            projectile.vz = v.vz;
        }
        else
        {
            setSpeedWithScatter(std::max(opt.min_energy_eV, E0));
        }
        local.E_after_eV = KineticEnergyEV3D(projectile.velocity, projectile.vz, m);
    }
    else if (kind == "excitation")
    {
        const double E1 = std::max(opt.min_energy_eV, E0 - std::max(0.0, ch.threshold_eV));
        setSpeedWithScatter(E1);
        local.E_after_eV = KineticEnergyEV3D(projectile.velocity, projectile.vz, m);
    }
    else if (kind == "ionization")
    {
        const double Eth = std::max(0.0, ch.threshold_eV);
        const double avail = std::max(0.0, E0 - Eth);

        // Count electron products.
        int n_e_prod = 0;
        std::vector<std::string> ion_products;
        std::vector<std::string> other_products;
        for (const auto& p : ch.products)
        {
            const Species* sp = species_db.Find(p);
            if (!sp)
                continue;
            if (sp->IsElectron())
                ++n_e_prod;
            else if (sp->IsIon())
                ion_products.push_back(p);
            else
                other_products.push_back(p);
        }

        const bool projIsElectron = projSp->IsElectron();
        int n_secondary = 0;
        if (projIsElectron)
        {
            // If the products list contains electrons, assume one corresponds to the scattered primary.
            n_secondary = std::max(0, n_e_prod - 1);
        }
        else
        {
            n_secondary = n_e_prod;
        }

        const std::string share_model = ChooseIonizationShareModel(opt, ch);
        double E_primary = std::max(opt.min_energy_eV, avail);
        double E_secondary = 0.0;

        if (projIsElectron)
        {
            const int n_total_e = 1 + n_secondary;
            if (n_total_e > 0)
            {
                if (share_model == "equal")
                {
                    const double each = avail / static_cast<double>(n_total_e);
                    E_primary = std::max(opt.min_energy_eV, each);
                    E_secondary = std::max(opt.min_energy_eV, each);
                }
                else
                {
                    // Minimal fallback: equal share.
                    const double each = avail / static_cast<double>(n_total_e);
                    E_primary = std::max(opt.min_energy_eV, each);
                    E_secondary = std::max(opt.min_energy_eV, each);
                }
            }
        }
        else
        {
            // Minimal: all loss assigned to threshold; secondaries get zero kinetic energy.
            E_primary = std::max(opt.min_energy_eV, E0 - Eth);
            E_secondary = opt.min_energy_eV;
        }

        // Update projectile (primary).
        setSpeedWithScatter(E_primary);

        // Spawn secondary electrons (if any).
        //
        // IMPORTANT: When the projectile is an electron, older dicts specify two
        // electron products (scattered primary + one secondary). In that case we
        // intentionally spawn only the *secondary* electron(s) and keep species/mass
        // equal to the projectile.
        //
        // When the projectile is *not* an electron (e.g., Ar+ impact ionization), any
        // electron listed in products must spawn an actual electron particle (not the
        // projectile species).
        if (spawn && n_secondary > 0)
        {
            const double wp = (projectile.w > 0.0) ? projectile.w : 1.0;
            std::string electron_name;
            double m_e = 0.0;
            if (!projIsElectron)
            {
                // Prefer the exact token provided in products (preserves naming like "e").
                for (const auto& p : ch.products)
                {
                    const Species* sp = species_db.Find(p);
                    if (sp && sp->IsElectron())
                    {
                        electron_name = p;
                        m_e = sp->mass_kg;
                        break;
                    }
                }
                // Fallback: a common canonical name.
                if (electron_name.empty())
                {
                    electron_name = "electron";
                    const Species* sp = species_db.Find(electron_name);
                    if (sp)
                        m_e = sp->mass_kg;
                }
            }

            for (int i = 0; i < n_secondary; ++i)
            {
                Particle sec = projectile;
                // Keep position/cell labels, but randomize direction.
                const double m_spawn = projIsElectron ? m : ((m_e > 0.0) ? m_e : m);
                const double vmag = SpeedFromEnergyEV(E_secondary, m_spawn);
                const auto dir = RandomIsotropicDirection3D(rng);
                sec.velocity = Vector2D(dir[0] * vmag, dir[1] * vmag);
                sec.vz = dir[2] * vmag;
                sec.w = wp;
                spawn->Push(projIsElectron ? set.projectile : electron_name, sec);
                ++local.n_spawned;
            }
        }

        // Spawn ions (at rest).
        if (spawn)
        {
            const double wp = (projectile.w > 0.0) ? projectile.w : 1.0;
            for (const auto& ion : ion_products)
            {
                Particle ip = projectile;
                ip.velocity = Vector2D(0.0, 0.0);
                ip.vz = 0.0;
                ip.w = wp;
                spawn->Push(ion, ip);
                ++local.n_spawned;
            }
        }

        // Other products: spawn at rest (minimal).
        if (spawn)
        {
            const double wp = (projectile.w > 0.0) ? projectile.w : 1.0;
            for (const auto& op : other_products)
            {
                if (ToLowerCopy(op) == ToLowerCopy(set.projectile))
                    continue;
                Particle pp = projectile;
                pp.velocity = Vector2D(0.0, 0.0);
                pp.vz = 0.0;
                pp.w = wp;
                spawn->Push(op, pp);
                ++local.n_spawned;
            }
        }

        local.E_after_eV = KineticEnergyEV3D(projectile.velocity, projectile.vz, m);
    }
    else if (kind == "charge_exchange")
    {
        // Minimal ion-neutral CX: keep projectile as the ion but reset its velocity;
        // spawn a fast neutral carrying the pre-CX velocity when requested by products.
        const Vector2D v_old = projectile.velocity;
        const double vz_old = projectile.vz;
        const bool targetIsGas = (!opt.gas_species_name.empty())
            && (ToLowerCopy(set.target) == ToLowerCopy(opt.gas_species_name));
        if (opt.ion_neutral_cx_thermalize && projSp->IsIon() && targetIsGas)
        {
            const Velocity3D v = SampleMaxwellianVelocity3D(m, opt.gas_temperature_K, /*drift=*/Vector2D(0.0, 0.0), rng);
            projectile.velocity = v.v;
            projectile.vz = v.vz;
        }
        else
        {
            // Fallback: preserve speed but randomize direction.
            setSpeedWithScatter(std::max(opt.min_energy_eV, E0));
        }

        if (spawn && opt.track_fast_neutrals)
        {
            const double wp = (projectile.w > 0.0) ? projectile.w : 1.0;
            for (const auto& p : ch.products)
            {
                if (ToLowerCopy(p) == ToLowerCopy(set.projectile))
                    continue;
                Particle sp = projectile;
                const Species* s = species_db.Find(p);
                if (s && s->IsNeutral())
                {
                    sp.velocity = v_old;
                    sp.vz = vz_old;
                }
                else
                {
                    sp.velocity = Vector2D(0.0, 0.0);
                    sp.vz = 0.0;
                }
                sp.w = wp;
                spawn->Push(p, sp);
                ++local.n_spawned;
            }
        }

        local.E_after_eV = KineticEnergyEV3D(projectile.velocity, projectile.vz, m);
    }
    else
    {
        // Unknown: treat as elastic.
        setSpeedWithScatter(std::max(opt.min_energy_eV, E0));
        local.E_after_eV = KineticEnergyEV3D(projectile.velocity, projectile.vz, m);
    }

    if (out)
        *out = local;
    return true;
}
} // namespace

bool PostCollisionModel::Apply(const ReactionSet& set, const ReactionChannel& ch,
                               Particle& projectile,
                               std::mt19937_64& rng,
                               const TargetState* target_state,
                               ParticleSpawnBuffer* spawn,
                               Outcome* out,
                               std::string* err) const
{
    return ApplyImpl(species_db_, opt_, set, ch, projectile, rng, target_state, spawn, out, err);
}

std::string PostCollisionModel::InferKind(const ReactionSet& set, const ReactionChannel& channel)
{
    return ::InferKind(set, channel);
}

bool PostCollisionModel::Apply(const ReactionSet& set, const ReactionChannel& ch,
                               Particle& projectile,
                               FastRng::SplitMix64Rng& rng,
                               const TargetState* target_state,
                               ParticleSpawnBuffer* spawn,
                               Outcome* out,
                               std::string* err) const
{
    return ApplyImpl(species_db_, opt_, set, ch, projectile, rng, target_state, spawn, out, err);
}

