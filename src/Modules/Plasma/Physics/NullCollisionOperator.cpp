#include "Modules/Plasma/Physics/NullCollisionOperator.h"

#include "Modules/Plasma/Physics/TabulatedCrossSection.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>

using namespace std;

void NullCollisionStats::Clear()
{
    n_particles = 0;
    n_steps = 0;
    n_attempts = 0;
    n_null = 0;
    n_collisions = 0;
    per_channel.clear();
}

void NullCollisionStats::Print(std::ostream& os) const
{
    os << "NullCollisionStats\n";
    os << "  particles    : " << n_particles << "\n";
    os << "  steps        : " << n_steps << "\n";
    os << "  attempts     : " << n_attempts << "\n";
    os << "  null         : " << n_null << "\n";
    os << "  collisions   : " << n_collisions << "\n";
    if (!per_channel.empty())
    {
        os << "  per-channel:\n";
        // Stable ordering by key.
        vector<pair<string, uint64_t>> entries(per_channel.begin(), per_channel.end());
        sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (const auto& kv : entries)
        {
            os << "    - " << kv.first << " : " << kv.second << "\n";
        }
    }
}

std::string NullCollisionStats::ToString() const
{
    ostringstream ss;
    Print(ss);
    return ss.str();
}

namespace
{
inline double safeExpMinus(double x)
{
    // Returns 1-exp(-x) for x>=0 with better accuracy for small x.
    if (x < 1e-5)
        return x - 0.5 * x * x;
    return 1.0 - exp(-x);
}

inline double clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

}

bool NullCollisionOperator::Initialize(const ReactionSet& set, const CrossSectionDB& xs_db,
                                       const Options& opt, std::string* err)
{
    set_ = set;
    opt_ = opt;
    channels_.clear();
    sigma_majorant_m2_ = 0.0;

    // Validate ReactionSet and resolve models.
    string tmp;
    if (!set_.Validate(&xs_db, &tmp))
    {
        if (err) *err = tmp;
        return false;
    }

    for (const auto& ch : set_.channels)
    {
        if (!ch.enabled)
            continue;
        ChannelRuntime rt;
        rt.name = ch.name;
        rt.threshold_eV = ch.threshold_eV;
        rt.model = xs_db.Get(ch.cross_section_key);
        if (!rt.model)
        {
            if (err) *err = "Missing cross section model for key: " + ch.cross_section_key;
            return false;
        }
        channels_.push_back(std::move(rt));
    }

    if (channels_.empty())
    {
        if (err) *err = "All channels are disabled; nothing to simulate";
        return false;
    }

    if (opt_.majorant_Emax_eV <= 0.0)
    {
        if (err) *err = "majorant_Emax_eV must be > 0";
        return false;
    }
    if (opt_.energy_samples < 4)
        opt_.energy_samples = 4;
    if (opt_.safety_factor < 1.0)
        opt_.safety_factor = 1.0;

    // Optional fast lookup: build per-channel hint tables for tabulated cross sections.
    // These tables map uniform energy bins to an initial lower_bound index guess, allowing
    // SigmaM2WithHint() to skip the per-call binary search in hot loops.
    for (auto& rt : channels_)
    {
        rt.tabulated = dynamic_cast<TabulatedCrossSection*>(rt.model.get());
        rt.hint_i1.clear();
        rt.hint_inv_dE = 0.0;
        rt.hint_bins = 0;
        rt.hint_Emax = 0.0;
    }

    if (opt_.use_fast_lookup && opt_.lookup_bins >= 8)
    {
        const double Emax = opt_.majorant_Emax_eV;
        const size_t bins = std::max<std::size_t>(16, opt_.lookup_bins);
        const double inv_dE = (Emax > 0.0) ? (static_cast<double>(bins) / Emax) : 0.0;

        for (auto& rt : channels_)
        {
            if (!rt.tabulated)
                continue;
            const auto& Egrid = rt.tabulated->energies_eV();
            if (Egrid.size() < 2)
                continue;

            rt.hint_bins = bins;
            rt.hint_Emax = Emax;
            rt.hint_inv_dE = inv_dE;
            rt.hint_i1.assign(bins + 1, 1u);

            for (size_t b = 0; b <= bins; ++b)
            {
                const double t = static_cast<double>(b) / static_cast<double>(bins);
                const double E = t * Emax;
                auto it = std::lower_bound(Egrid.begin(), Egrid.end(), E);
                size_t i1 = static_cast<size_t>(it - Egrid.begin());
                if (i1 == 0)
                    i1 = 1;
                if (i1 >= Egrid.size())
                    i1 = Egrid.size() - 1;
                rt.hint_i1[b] = static_cast<std::uint32_t>(i1);
            }
        }
    }

    return RecomputeMajorant(err);
}

double NullCollisionOperator::ClampEnergy(double E_eV) const
{
    if (E_eV < 0.0) return 0.0;
    if (!opt_.clamp_energy) return E_eV;
    if (E_eV > opt_.majorant_Emax_eV) return opt_.majorant_Emax_eV;
    return E_eV;
}

bool NullCollisionOperator::RecomputeMajorant(std::string* err)
{
    // Sample sigma_tot over [0, Emax] and take max.
    const double Emax = opt_.majorant_Emax_eV;
    const size_t N = opt_.energy_samples;

    double sigma_max = 0.0;
    // Always include 0 and Emax.
    for (size_t i = 0; i < N; ++i)
    {
        const double t = static_cast<double>(i) / static_cast<double>(N - 1);
        const double E = t * Emax;
        const double s = SigmaTotalM2(E);
        if (!(s >= 0.0) || !isfinite(s))
        {
            if (err) *err = "Non-finite sigma_tot during majorant scan";
            return false;
        }
        sigma_max = max(sigma_max, s);
    }

    if (sigma_max <= 0.0)
    {
        if (err) *err = "sigma_tot max is non-positive; check cross sections";
        return false;
    }

    sigma_majorant_m2_ = sigma_max * opt_.safety_factor;
    if (!(sigma_majorant_m2_ > 0.0) || !isfinite(sigma_majorant_m2_))
    {
        if (err) *err = "Invalid sigma_majorant computed";
        return false;
    }
    return true;
}

double NullCollisionOperator::SigmaTotalM2(double E_eV) const
{
    const double E = ClampEnergy(E_eV);
    double sum = 0.0;
    for (const auto& ch : channels_)
    {
        // Below threshold, treat as zero probability.
        if (E < ch.threshold_eV)
            continue;
        double s = 0.0;
        if (opt_.use_fast_lookup && ch.tabulated && ch.HasHintTable())
        {
            std::size_t i1 = ch.HintI1ForEnergy(E);
            s = ch.tabulated->SigmaM2WithHint(E, &i1);
        }
        else
        {
            s = ch.model->SigmaM2(E);
        }
        if (s > 0.0)
            sum += s;
    }
    return sum;
}

int NullCollisionOperator::SampleChannelIndex(double E_eV, std::mt19937_64& rng) const
{
    const double E = ClampEnergy(E_eV);
    // Hot path: avoid per-call allocations by reusing a thread-local scratch vector.
    static thread_local vector<double> weights;
    const size_t n = channels_.size();
    if (weights.size() < n)
        weights.resize(n);
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const auto& ch = channels_[i];
        double w = 0.0;
        if (E >= ch.threshold_eV)
        {
            double s = 0.0;
            if (opt_.use_fast_lookup && ch.tabulated && ch.HasHintTable())
            {
                std::size_t i1 = ch.HintI1ForEnergy(E);
                s = ch.tabulated->SigmaM2WithHint(E, &i1);
            }
            else
            {
                s = ch.model->SigmaM2(E);
            }
            if (s > 0.0 && isfinite(s))
                w = s;
        }
        weights[i] = w;
        sum += w;
    }

    if (!(sum > 0.0))
        return -1;

    static thread_local uniform_real_distribution<double> uni(0.0, 1.0);
    const double r = uni(rng) * sum;
    double acc = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        acc += weights[i];
        if (r <= acc)
            return static_cast<int>(i);
    }
    return static_cast<int>(n - 1);
}

int NullCollisionOperator::SampleChannelIndex(double E_eV, FastRng::SplitMix64Rng& rng) const
{
    const double E = ClampEnergy(E_eV);
    // Hot path: avoid per-call allocations by reusing a thread-local scratch vector.
    static thread_local vector<double> weights;
    const size_t n = channels_.size();
    if (weights.size() < n)
        weights.resize(n);
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const auto& ch = channels_[i];
        double w = 0.0;
        if (E >= ch.threshold_eV)
        {
            double s = 0.0;
            if (opt_.use_fast_lookup && ch.tabulated && ch.HasHintTable())
            {
                std::size_t i1 = ch.HintI1ForEnergy(E);
                s = ch.tabulated->SigmaM2WithHint(E, &i1);
            }
            else
            {
                s = ch.model->SigmaM2(E);
            }
            if (s > 0.0 && isfinite(s))
                w = s;
        }
        weights[i] = w;
        sum += w;
    }

    if (!(sum > 0.0))
        return -1;

    static thread_local uniform_real_distribution<double> uni(0.0, 1.0);
    const double r = uni(rng) * sum;
    double acc = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        acc += weights[i];
        if (r <= acc)
            return static_cast<int>(i);
    }
    return static_cast<int>(n - 1);
}

std::uint32_t NullCollisionOperator::ApplyOne(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                              std::mt19937_64& rng, NullCollisionStats* stats) const
{
    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return 0;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        if (stats)
        {
            stats->n_particles += 1;
            stats->n_steps += 1;
        }
        return 0;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
        return 0;

    uniform_real_distribution<double> uni(0.0, 1.0);

    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps; // 1/s
    if (!(nu_max > 0.0) || !isfinite(nu_max))
        return 0;

    if (stats)
    {
        stats->n_particles += 1;
        stats->n_steps += 1;
    }

    // Poisson process under majorant: draw waiting times and allow multiple events per dt.
    double t_rem = dt_s;
    uint32_t accepted = 0;
    while (t_rem > 0.0)
    {
        // Exponential waiting time with rate nu_max.
        const double u = max(uni(rng), numeric_limits<double>::min());
        const double tau = -log(u) / nu_max;
        if (tau > t_rem)
            break;
        t_rem -= tau;

        if (stats)
            stats->n_attempts += 1;

        const double sigma_tot = SigmaTotalM2(E_eV);
        // Robust ratio.
        const double accept_prob = clamp01(sigma_tot / sigmaM);
        if (uni(rng) > accept_prob)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }

        const int idx = SampleChannelIndex(E_eV, rng);
        if (idx < 0)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }
        accepted += 1;
        if (stats)
        {
            stats->n_collisions += 1;
            stats->per_channel[channels_[static_cast<size_t>(idx)].name] += 1;
        }
    }

    return accepted;
}

std::uint32_t NullCollisionOperator::ApplyOne(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                              FastRng::SplitMix64Rng& rng, NullCollisionStats* stats) const
{
    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return 0;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        if (stats)
        {
            stats->n_particles += 1;
            stats->n_steps += 1;
        }
        return 0;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
        return 0;

    uniform_real_distribution<double> uni(0.0, 1.0);

    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps; // 1/s
    if (!(nu_max > 0.0) || !isfinite(nu_max))
        return 0;

    if (stats)
    {
        stats->n_particles += 1;
        stats->n_steps += 1;
    }

    double t_rem = dt_s;
    uint32_t accepted = 0;
    while (t_rem > 0.0)
    {
        const double u = max(uni(rng), numeric_limits<double>::min());
        const double tau = -log(u) / nu_max;
        if (tau > t_rem)
            break;
        t_rem -= tau;

        if (stats)
            stats->n_attempts += 1;

        const double sigma_tot = SigmaTotalM2(E_eV);
        // Robust ratio.
        const double accept_prob = clamp01(sigma_tot / sigmaM);
        if (uni(rng) > accept_prob)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }

        const int idx = SampleChannelIndex(E_eV, rng);
        if (idx < 0)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }
        accepted += 1;
        if (stats)
        {
            stats->n_collisions += 1;
            stats->per_channel[channels_[static_cast<size_t>(idx)].name] += 1;
        }
    }

    return accepted;
}

std::uint32_t NullCollisionOperator::ApplyOneCollect(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                                     std::mt19937_64& rng,
                                                     std::vector<int>* accepted_channel_indices,
                                                     NullCollisionStats* stats) const
{
    if (accepted_channel_indices)
        accepted_channel_indices->clear();

    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return 0;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        if (stats)
        {
            stats->n_particles += 1;
            stats->n_steps += 1;
        }
        return 0;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
        return 0;

    uniform_real_distribution<double> uni(0.0, 1.0);

    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps; // 1/s
    if (!(nu_max > 0.0) || !isfinite(nu_max))
        return 0;

    if (stats)
    {
        stats->n_particles += 1;
        stats->n_steps += 1;
    }

    double t_rem = dt_s;
    uint32_t accepted = 0;
    while (t_rem > 0.0)
    {
        const double u = max(uni(rng), numeric_limits<double>::min());
        const double tau = -log(u) / nu_max;
        if (tau > t_rem)
            break;
        t_rem -= tau;

        if (stats)
            stats->n_attempts += 1;

        const double sigma_tot = SigmaTotalM2(E_eV);
        const double accept_prob = clamp01(sigma_tot / sigmaM);
        if (uni(rng) > accept_prob)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }

        const int idx = SampleChannelIndex(E_eV, rng);
        if (idx < 0)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }

        accepted += 1;
        if (accepted_channel_indices)
            accepted_channel_indices->push_back(idx);
        if (stats)
        {
            stats->n_collisions += 1;
            stats->per_channel[channels_[static_cast<size_t>(idx)].name] += 1;
        }
    }
    return accepted;
}

std::uint32_t NullCollisionOperator::ApplyOneCollect(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                                     FastRng::SplitMix64Rng& rng,
                                                     std::vector<int>* accepted_channel_indices,
                                                     NullCollisionStats* stats) const
{
    if (accepted_channel_indices)
        accepted_channel_indices->clear();

    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return 0;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        if (stats)
        {
            stats->n_particles += 1;
            stats->n_steps += 1;
        }
        return 0;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
        return 0;

    uniform_real_distribution<double> uni(0.0, 1.0);

    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps; // 1/s
    if (!(nu_max > 0.0) || !isfinite(nu_max))
        return 0;

    if (stats)
    {
        stats->n_particles += 1;
        stats->n_steps += 1;
    }

    double t_rem = dt_s;
    uint32_t accepted = 0;
    while (t_rem > 0.0)
    {
        const double u = max(uni(rng), numeric_limits<double>::min());
        const double tau = -log(u) / nu_max;
        if (tau > t_rem)
            break;
        t_rem -= tau;

        if (stats)
            stats->n_attempts += 1;

        const double sigma_tot = SigmaTotalM2(E_eV);
        const double accept_prob = clamp01(sigma_tot / sigmaM);
        if (uni(rng) > accept_prob)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }

        const int idx = SampleChannelIndex(E_eV, rng);
        if (idx < 0)
        {
            if (stats)
                stats->n_null += 1;
            continue;
        }

        accepted += 1;
        if (accepted_channel_indices)
            accepted_channel_indices->push_back(idx);
        if (stats)
        {
            stats->n_collisions += 1;
            stats->per_channel[channels_[static_cast<size_t>(idx)].name] += 1;
        }
    }
    return accepted;
}


std::uint32_t NullCollisionOperator::ApplyOneCollectCounts(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                                           std::mt19937_64& rng,
                                                           std::vector<int>* accepted_channel_indices,
                                                           NullCollisionCounts* counts) const
{
    if (accepted_channel_indices)
        accepted_channel_indices->clear();

    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return 0;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        if (counts)
        {
            counts->n_particles += 1;
            counts->n_steps += 1;
        }
        return 0;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
        return 0;

    uniform_real_distribution<double> uni(0.0, 1.0);

    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps; // 1/s
    if (!(nu_max > 0.0) || !isfinite(nu_max))
        return 0;

    if (counts)
    {
        counts->n_particles += 1;
        counts->n_steps += 1;
    }

    double t_rem = dt_s;
    uint32_t accepted = 0;
    while (t_rem > 0.0)
    {
        const double u = max(uni(rng), numeric_limits<double>::min());
        const double tau = -log(u) / nu_max;
        if (tau > t_rem)
            break;
        t_rem -= tau;

        if (counts)
            counts->n_attempts += 1;

        const double sigma_tot = SigmaTotalM2(E_eV);
        const double accept_prob = clamp01(sigma_tot / sigmaM);
        if (uni(rng) > accept_prob)
        {
            if (counts)
                counts->n_null += 1;
            continue;
        }

        const int idx = SampleChannelIndex(E_eV, rng);
        if (idx < 0)
        {
            if (counts)
                counts->n_null += 1;
            continue;
        }

        accepted += 1;
        if (accepted_channel_indices)
            accepted_channel_indices->push_back(idx);
        if (counts)
            counts->n_collisions += 1;
    }
    return accepted;
}

std::uint32_t NullCollisionOperator::ApplyOneCollectCounts(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                                           FastRng::SplitMix64Rng& rng,
                                                           std::vector<int>* accepted_channel_indices,
                                                           NullCollisionCounts* counts) const
{
    if (accepted_channel_indices)
        accepted_channel_indices->clear();

    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return 0;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        if (counts)
        {
            counts->n_particles += 1;
            counts->n_steps += 1;
        }
        return 0;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
        return 0;

    uniform_real_distribution<double> uni(0.0, 1.0);

    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps; // 1/s
    if (!(nu_max > 0.0) || !isfinite(nu_max))
        return 0;

    if (counts)
    {
        counts->n_particles += 1;
        counts->n_steps += 1;
    }

    double t_rem = dt_s;
    uint32_t accepted = 0;
    while (t_rem > 0.0)
    {
        const double u = max(uni(rng), numeric_limits<double>::min());
        const double tau = -log(u) / nu_max;
        if (tau > t_rem)
            break;
        t_rem -= tau;

        if (counts)
            counts->n_attempts += 1;

        const double sigma_tot = SigmaTotalM2(E_eV);
        const double accept_prob = clamp01(sigma_tot / sigmaM);
        if (uni(rng) > accept_prob)
        {
            if (counts)
                counts->n_null += 1;
            continue;
        }

        const int idx = SampleChannelIndex(E_eV, rng);
        if (idx < 0)
        {
            if (counts)
                counts->n_null += 1;
            continue;
        }

        accepted += 1;
        if (accepted_channel_indices)
            accepted_channel_indices->push_back(idx);
        if (counts)
            counts->n_collisions += 1;
    }
    return accepted;
}

void NullCollisionOperator::SampleNextEvent(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                            std::mt19937_64& rng,
                                            NextEvent* event,
                                            NullCollisionCounts* counts) const
{
    if (!event)
        return;
    *event = NextEvent{};

    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        event->elapsed_s = dt_s;
        return;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
    {
        event->elapsed_s = dt_s;
        return;
    }

    uniform_real_distribution<double> uni(0.0, 1.0);
    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps;
    if (!(nu_max > 0.0) || !isfinite(nu_max))
    {
        event->elapsed_s = dt_s;
        return;
    }

    const double u = max(uni(rng), numeric_limits<double>::min());
    const double tau = -log(u) / nu_max;
    if (tau > dt_s)
    {
        event->elapsed_s = dt_s;
        return;
    }

    event->elapsed_s = tau;
    if (counts)
        counts->n_attempts += 1;

    const double sigma_tot = SigmaTotalM2(E_eV);
    const double accept_prob = clamp01(sigma_tot / sigmaM);
    if (uni(rng) > accept_prob)
    {
        event->kind = NextEvent::Kind::kNull;
        if (counts)
            counts->n_null += 1;
        return;
    }

    const int idx = SampleChannelIndex(E_eV, rng);
    if (idx < 0)
    {
        event->kind = NextEvent::Kind::kNull;
        if (counts)
            counts->n_null += 1;
        return;
    }

    event->kind = NextEvent::Kind::kCollision;
    event->channel_index = idx;
    if (counts)
        counts->n_collisions += 1;
}

void NullCollisionOperator::SampleNextEvent(double n_gas_m3, double v_rel_mps, double E_eV, double dt_s,
                                            FastRng::SplitMix64Rng& rng,
                                            NextEvent* event,
                                            NullCollisionCounts* counts) const
{
    if (!event)
        return;
    *event = NextEvent{};

    if (!(n_gas_m3 >= 0.0) || !(v_rel_mps >= 0.0) || !(dt_s >= 0.0))
        return;
    if (dt_s == 0.0 || n_gas_m3 == 0.0 || v_rel_mps == 0.0)
    {
        event->elapsed_s = dt_s;
        return;
    }

    const double sigmaM = sigma_majorant_m2_;
    if (!(sigmaM > 0.0))
    {
        event->elapsed_s = dt_s;
        return;
    }

    uniform_real_distribution<double> uni(0.0, 1.0);
    const double nu_max = n_gas_m3 * sigmaM * v_rel_mps;
    if (!(nu_max > 0.0) || !isfinite(nu_max))
    {
        event->elapsed_s = dt_s;
        return;
    }

    const double u = max(uni(rng), numeric_limits<double>::min());
    const double tau = -log(u) / nu_max;
    if (tau > dt_s)
    {
        event->elapsed_s = dt_s;
        return;
    }

    event->elapsed_s = tau;
    if (counts)
        counts->n_attempts += 1;

    const double sigma_tot = SigmaTotalM2(E_eV);
    const double accept_prob = clamp01(sigma_tot / sigmaM);
    if (uni(rng) > accept_prob)
    {
        event->kind = NextEvent::Kind::kNull;
        if (counts)
            counts->n_null += 1;
        return;
    }

    const int idx = SampleChannelIndex(E_eV, rng);
    if (idx < 0)
    {
        event->kind = NextEvent::Kind::kNull;
        if (counts)
            counts->n_null += 1;
        return;
    }

    event->kind = NextEvent::Kind::kCollision;
    event->channel_index = idx;
    if (counts)
        counts->n_collisions += 1;
}
