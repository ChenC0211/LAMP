#include "Modules/Plasma/IO/PicCollisionDiagnostics.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/MCCCollisionOperator.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Core/IO/ConsoleControl.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

void writeExplicitPicCollisionInitialization(const MCCCollisionOperator& operation,
    const GasState& gas, const UnitsContext* units, const NormalizationScales& sc,
    const QuadtreeMesh& mesh, const std::vector<SpeciesConfig>& species,
    double dt, const std::string& reportPath)
{
    {
        auto fmtSci = [](double v, int p = 3) {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(p) << v;
            return oss.str();
        };

        const double dt_si = sc.T0_s * dt;
        const double v0_mps = sc.v0_mps;
        const double n_gas_m3 = gas.number_density_m3;

        // CFL estimate: max( |v| * dt / h_min ), using |drift| + 3*v_th as v_max.
        double h_min_code = std::numeric_limits<double>::max();
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (c.son[0] != -1)
                continue;
            h_min_code = std::min(h_min_code, std::min(c.size().x, c.size().y));
        }
        double h_min_si = 0.0;
        bool h_min_valid = false;
        if (std::isfinite(h_min_code) && h_min_code > 0.0)
        {
            if (units)
                h_min_si = units->LengthCodeToSI(h_min_code);
            else if (sc.L0_m > 0.0)
                h_min_si = sc.L0_m * h_min_code;
            if (std::isfinite(h_min_si) && h_min_si > 0.0)
                h_min_valid = true;
        }

        double vmax_si = 0.0;
        std::string vmax_species;
        const double vth_factor = 3.0;
        for (const auto& sp : species)
        {
            const double vdrift = std::sqrt(sp.drift.x * sp.drift.x + sp.drift.y * sp.drift.y);
            double vth = 0.0;
            if (sp.m > 0.0 && sp.Te_eV > 0.0)
            {
                const double kT_J = sp.Te_eV * PhysicalConstants::eV_J;
                vth = std::sqrt(2.0 * kT_J / sp.m);
            }
            const double v_est = vdrift + vth_factor * vth;
            if (v_est > vmax_si)
            {
                vmax_si = v_est;
                vmax_species = sp.name;
            }
        }
        double cfl = 0.0;
        if (h_min_valid && std::isfinite(vmax_si) && vmax_si > 0.0 && std::isfinite(dt_si))
            cfl = vmax_si * dt_si / h_min_si;

        double sigma_max = 0.0;
        double nu_dt_max = 0.0;
        std::string nu_dt_proj;
        for (const auto& info : operation.MajorantTable())
        {
            sigma_max = std::max(sigma_max, info.sigma_majorant_m2);
            const double nu_dt = n_gas_m3 * v0_mps * info.sigma_majorant_m2 * dt_si;
            if (nu_dt > nu_dt_max)
            {
                nu_dt_max = nu_dt;
                nu_dt_proj = info.projectile;
            }
        }

        std::cout << "[ExplicitPIC][MCC] unit-check: dt_s=" << fmtSci(dt_si)
                  << " v0_mps=" << fmtSci(v0_mps)
                  << " n_gas=" << fmtSci(n_gas_m3)
                  << " sigmaM_max=" << fmtSci(sigma_max)
                  << " nu*dt_max=" << fmtSci(nu_dt_max);
        if (h_min_valid && std::isfinite(vmax_si) && vmax_si > 0.0)
        {
            std::cout << " h_min=" << fmtSci(h_min_si)
                      << " vmax_est=" << fmtSci(vmax_si)
                      << " CFL=" << fmtSci(cfl);
        }
        if (!nu_dt_proj.empty())
            std::cout << " projectile=" << nu_dt_proj;
        std::cout << std::endl;

        if (units && units->mode() != UnitsMode::kSI)
        {
            std::cerr << "[ExplicitPIC][MCC] WARNING: Units.mode is not SI; MCC expects SI gas/cross-sections. "
                         "Ensure code units are SI or set Units.mode=SI.\n";
        }
        if (!std::isfinite(dt_si) || dt_si <= 0.0)
        {
            std::cerr << "[ExplicitPIC][MCC] WARNING: dt_s is non-positive/invalid; check Units/Normalization.\n";
        }
        if (!std::isfinite(v0_mps) || v0_mps <= 0.0)
        {
            std::cerr << "[ExplicitPIC][MCC] WARNING: v0_mps is invalid; check Units/Normalization.\n";
        }
        if (!std::isfinite(n_gas_m3) || n_gas_m3 < 0.0)
        {
            std::cerr << "[ExplicitPIC][MCC] WARNING: n_gas is invalid; check Gas settings.\n";
        }
        if (nu_dt_max > 1.0)
        {
            std::cerr << "[ExplicitPIC][MCC] WARNING: max(nu*dt) > 1.0; check dt and unit consistency.\n";
        }
        if (h_min_valid && std::isfinite(vmax_si) && vmax_si > 0.0 && cfl > 1.0)
        {
            std::cerr << "[ExplicitPIC][CFL] WARNING: vmax*dt/h_min = " << fmtSci(cfl)
                      << " (>1). Consider reducing dt or refining mesh.\n";
            if (!vmax_species.empty())
                std::cerr << "[ExplicitPIC][CFL] WARNING: vmax estimated from species '" << vmax_species
                          << "' (|drift| + 3*v_th, v_th=sqrt(2kT/m)).\n";
        }
    }

    // Initialize collisions report (fresh run).
    {
        std::ofstream out(reportPath, std::ios::out | std::ios::trunc);
        if (!out)
        {
            std::cerr << "[ExplicitPIC] WARNING: failed to open collisions report: " << reportPath << std::endl;
        }
        else
        {
            out << "# collisions_report\n";
            out << "# (all SI units unless noted)\n";
            out << "# gas: " << gas.name
                << "  p=" << gas.pressure_Pa << " Pa"
                << "  T=" << gas.temperature_K << " K"
                << "  n=" << gas.number_density_m3 << " 1/m^3\n";
            // Prefix the MCC summary with '# ' so the report stays comment-friendly.
            std::istringstream iss(operation.SummaryString());
            std::string line;
            while (std::getline(iss, line))
                out << "# " << line << "\n";
            out << "\n";
        }
    }
}

namespace
{
std::vector<std::pair<std::string, std::uint64_t>> sortedPairsByCountDesc(const std::unordered_map<std::string, std::uint64_t>& m)
{
    std::vector<std::pair<std::string, std::uint64_t>> v;
    v.reserve(m.size());
    for (const auto& kv : m)
        v.emplace_back(kv.first, kv.second);
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second)
            return a.second > b.second;
        return a.first < b.first;
    });
    return v;
}

void writeMccReportBlock(std::ostream& os,
                               const std::string& title,
                               int stepForPrint,
                               double timeForPrint,
                               const MCCCollisionOperator::StepStats& st)
{
    os << title << " step=" << stepForPrint << " t=" << std::setprecision(16) << timeForPrint << "\n";
    os << "  totals: particles=" << st.n_particles
       << " attempts=" << st.n_attempts
       << " null=" << st.n_null
       << " collisions=" << st.n_collisions
       << " spawned=" << st.n_spawned
       << " dE_macro_eV=" << st.energy_loss_eV_macro
       << " dE_total_J=" << st.energy_loss_J_total << " max_nu_dtM=" << st.max_nu_dt_majorant << "\n";

    if (!st.per_projectile.empty())
    {
        os << "  per_projectile:\n";
        os << "    " << std::left << std::setw(16) << "projectile"
           << std::right << std::setw(12) << "particles"
           << std::setw(12) << "attempts"
           << std::setw(10) << "null"
           << std::setw(12) << "coll"
           << std::setw(12) << "spawn"
           << std::setw(16) << "dE_eV"
           << std::setw(16) << "dE_J" << std::setw(14) << "max_nu_dtM" << "\n";

        std::vector<std::string> keys;
        keys.reserve(st.per_projectile.size());
        for (const auto& kv : st.per_projectile)
            keys.push_back(kv.first);
        std::sort(keys.begin(), keys.end());

        for (const auto& k : keys)
        {
            const auto& ps = st.per_projectile.at(k);
            os << "    " << std::left << std::setw(16) << k
               << std::right << std::setw(12) << ps.n_particles
               << std::setw(12) << ps.n_attempts
               << std::setw(10) << ps.n_null
               << std::setw(12) << ps.n_collisions
               << std::setw(12) << ps.n_spawned
               << std::setw(16) << ps.energy_loss_eV_macro
               << std::setw(16) << ps.energy_loss_J_total << std::setw(14) << ps.max_nu_dt_majorant << "\n";

            if (!ps.per_channel.empty())
            {
                os << "      channels:\n";
                for (const auto& kv : sortedPairsByCountDesc(ps.per_channel))
                    os << "        " << kv.first << " : " << kv.second << "\n";
            }
        }
    }

    if (!st.per_channel.empty())
    {
        os << "  per_channel_total (key=projectile::channel):\n";
        for (const auto& kv : sortedPairsByCountDesc(st.per_channel))
            os << "    " << kv.first << " : " << kv.second << "\n";
    }

    if (!st.spawn_per_species.empty())
    {
        os << "  spawn_total_by_species:\n";
        for (const auto& kv : sortedPairsByCountDesc(st.spawn_per_species))
            os << "    " << kv.first << " : " << kv.second << "\n";
    }
}
}

ExplicitPicCollisionReports::ExplicitPicCollisionReports(const std::filesystem::path& metaDir, int statsEvery)
    : reportPath_((metaDir / "collisions_report.txt").string()), statsEvery_(statsEvery)
{
    const auto ionReportPath = (metaDir / "ion_report.txt").string();
    const auto neutralReportPath = (metaDir / "neutral_report.txt").string();
    // Ion-neutral collision slow-down report.
    ionOut_.open(ionReportPath, std::ios::out | std::ios::trunc);
    if (ionOut_)
    {
        ionOut_ << "# ion_report (ion-neutral ELASTIC+CX slow-down; optional fast-neutral tracking)\n";
        ionOut_ << "# columns: step,time_s,projectile,elastic,cx,collisions,energyLoss_J_total,max_nu_dt_majorant\n";
    }
    else
    {
        std::cerr << "[ExplicitPIC] WARNING: failed to create ion report: " << ionReportPath << std::endl;
    }

    // Fast-neutral tracking report.
    neutralOut_.open(neutralReportPath, std::ios::out | std::ios::trunc);
    if (neutralOut_)
    {
        neutralOut_ << "# neutral_report (lightweight fast-neutral tracking; q=0, excluded from Poisson)\n";
        neutralOut_ << "# columns: step,time_s,species,macroCount,physCount,meanEnergy_eV,totalEnergy_J,spawned_macro\n";
    }
    else
    {
        std::cerr << "[ExplicitPIC] WARNING: failed to create neutral report: " << neutralReportPath << std::endl;
    }

}

void ExplicitPicCollisionReports::writeStep(int step, double timeSI, const PicSpeciesState& particles,
    double velocityScale, const MCCCollisionOperator* operation, const MCCCollisionOperator::StepStats& stats)
{
    if (statsEvery_ <= 0 || step % statsEvery_ != 0) return;
    // Append ion-neutral ELASTIC/CX slow-down diagnostics.
    if (ionOut_.is_open() && ionOut_ && operation)
    {
        for (const auto& kv : stats.per_projectile)
        {
            const std::string& proj = kv.first;
            if (particles.ions().find(proj) == particles.ions().end())
                continue;
            const auto& ps = kv.second;

            uint64_t n_elastic = 0;
            uint64_t n_cx = 0;
            for (const auto& channel : ps.per_channel)
            {
                const auto kind = operation->ChannelKind(proj, channel.first);
                if (kind == "elastic") n_elastic += channel.second;
                else if (kind == "charge_exchange") n_cx += channel.second;
            }

            ionOut_ << step << ','
                   << std::setprecision(16) << timeSI << ','
                   << proj << ','
                   << n_elastic << ','
                   << n_cx << ','
                   << ps.n_collisions << ','
                   << std::setprecision(16) << ps.energy_loss_J_total << ','
                   << std::setprecision(16) << ps.max_nu_dt_majorant
                   << "\n";
        }
        ionOut_.flush();
    }

    // Append fast-neutral diagnostics. Neutrals are tracked as
    // regular particles with q=0 (kappa_s=0) and are excluded from Poisson.
    if (neutralOut_.is_open() && neutralOut_)
    {
        for (const auto& name : particles.neutrals())
        {
            auto it = particles.index().find(name);
            if (it == particles.index().end())
                continue;
            const std::size_t si = it->second;
            const auto& pt = particles.species()[si];
            const double defaultWeight = (si < particles.weights().size()) ? particles.weights()[si] : 1.0;
            const auto [n_macro, n_phys, Etot_J, meanE_eV] = summarizePicNeutralEnergy(pt, defaultWeight, velocityScale);

            uint64_t spawned = 0;
            if (operation)
            {
                auto sit = stats.spawn_per_species.find(name);
                if (sit != stats.spawn_per_species.end())
                    spawned = sit->second;
            }

            neutralOut_ << step << ','
                       << std::setprecision(16) << timeSI << ','
                       << name << ','
                       << n_macro << ','
                       << std::setprecision(16) << n_phys << ','
                       << std::setprecision(16) << meanE_eV << ','
                       << std::setprecision(16) << Etot_J << ','
                       << spawned
                       << "\n";
        }
        neutralOut_.flush();
    }

}

void ExplicitPicCollisionReports::append(int step, double timeSI, const PicSpeciesState& particles,
    double velocityScale, const MCCCollisionOperator::StepStats& stats,
    const MCCCollisionOperator::StepStats& cumulative)
{
    const auto [electronMean_eV, electronTailFrac, electronNmacro, electronNphys] =
        summarizePicCollisionElectronEnergy(particles.species(), particles.weights(), velocityScale);
    // stdout (debug/trace)
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
    {
        std::cout << "[ExplicitPIC][MCC] report interval reached" << std::endl;
        writeMccReportBlock(std::cout, "[MCC] step", step, timeSI, stats);
        writeMccReportBlock(std::cout, "[MCC] cumulative", step, timeSI, cumulative);
        std::cout << "[MCC] electron_mean_energy_eV=" << std::setprecision(10) << electronMean_eV
                  << " electron_tail_frac_E_gt_11.59eV=" << std::setprecision(10) << electronTailFrac
                  << " electron_n_macro=" << electronNmacro
                  << " electron_n_phys=" << std::setprecision(10) << electronNphys << std::endl;
    }

    // file
    std::ofstream out(reportPath_, std::ios::out | std::ios::app);
    if (!out)
    {
        std::cerr << "[ExplicitPIC] WARNING: failed to append collisions report: " << reportPath_ << std::endl;
        return;
    }
    out << "=== MCC report step=" << step << " t=" << std::setprecision(16) << timeSI << " ===\n";
    writeMccReportBlock(out, "step", step, timeSI, stats);
    writeMccReportBlock(out, "cumulative", step, timeSI, cumulative);

    out << "electron_mean_energy_eV=" << std::setprecision(16) << electronMean_eV << "\n";
    out << "electron_tail_frac_E_gt_11.59eV=" << std::setprecision(16) << electronTailFrac << "\n";
    out << "electron_n_macro=" << electronNmacro << "\n";
    out << "electron_n_phys=" << std::setprecision(16) << electronNphys << "\n";
    out << "\n";
}
