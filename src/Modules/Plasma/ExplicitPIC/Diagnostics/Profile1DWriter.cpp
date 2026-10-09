#include "Modules/Plasma/ExplicitPIC/Diagnostics/Profile1DWriter.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <unordered_set>

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"

namespace ExplicitPIC
{
namespace Diagnostics
{

static void depositCSR(const CouplingWeightsCache& cache,
                       const std::vector<double>& wp,
                       CellField<double>& out)
{
    const auto& cellIndices = cache.cellIndices();
    const auto& weights = cache.weights();
    const auto& offsets = cache.offsets();

    const std::size_t np = wp.size();
    if (offsets.size() != np + 1) return;

    for (std::size_t p = 0; p < np; ++p)
    {
        const double a = wp[p];
        const std::size_t begin = offsets[p];
        const std::size_t end = offsets[p + 1];
        for (std::size_t k = begin; k < end; ++k)
        {
            const int ci = cellIndices[k];
            if (ci >= 0 && ci < static_cast<int>(out.size()))
                out.value[ci] += a * weights[k];
        }
    }
}

void Profile1DWriter::open(const std::string& filePath, int startStep)
{
    close();
    cfg_.path = filePath;
    try
    {
        const std::filesystem::path p(filePath);
        if (p.has_parent_path())
            std::filesystem::create_directories(p.parent_path());
    }
    catch (...)
    {
        // best effort: directory creation failures will surface on open
    }

    headerWritten_ = TimeSeriesIO::openPreparedStepCsv(out_, std::filesystem::path(filePath), startStep);
}

void Profile1DWriter::close()
{
    if (out_.is_open())
        out_.close();
    headerWritten_ = false;
}

std::string Profile1DWriter::safeKey(const std::string& name)
{
    std::string k = ConsoleControl::toLowerCopy(name);
    for (char& c : k)
    {
        const bool ok = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || (c == '_');
        if (!ok) c = '_';
    }
    // trim duplicate underscores
    std::string out;
    out.reserve(k.size());
    bool prevUnderscore = false;
    for (char c : k)
    {
        if (c == '_')
        {
            if (!prevUnderscore) out.push_back(c);
            prevUnderscore = true;
        }
        else
        {
            out.push_back(c);
            prevUnderscore = false;
        }
    }
    // avoid leading underscore
    while (!out.empty() && out.front() == '_') out.erase(out.begin());
    while (!out.empty() && out.back() == '_') out.pop_back();
    if (out.empty()) out = "species";
    return out;
}

void Profile1DWriter::writeHeader(const std::vector<ParticleType>& species, bool writeSIOutputs)
{
    if (!out_.is_open()) return;

    std::unordered_set<std::string> used;
    std::vector<std::string> keys;
    keys.reserve(species.size());
    for (const auto& s : species)
    {
        std::string k = safeKey(s.name());
        if (used.count(k))
        {
            int n = 2;
            while (used.count(k + "_" + std::to_string(n))) ++n;
            k += "_" + std::to_string(n);
        }
        used.insert(k);
        keys.push_back(k);
    }

    out_ << "step,time_code";
    if (writeSIOutputs) out_ << ",time_SI_s";

    out_ << ",bin,x_code";
    if (writeSIOutputs) out_ << ",x_m";

    if (cfg_.includePhi)
    {
        out_ << ",phi_code";
        if (writeSIOutputs) out_ << ",phi_V";
    }
    if (cfg_.includeE)
    {
        out_ << ",Ex_code,Ey_code";
        if (writeSIOutputs) out_ << ",Ex_Vpm,Ey_Vpm";
    }
    if (cfg_.includeRho)
    {
        out_ << ",rho_tilde";
        if (writeSIOutputs) out_ << ",rho_SI_Cpm3";
    }

    for (const auto& k : keys)
        out_ << ",n_" << k << "_1pm3";

    out_ << ",Emean_e_eV,Te_e_eV";

    out_ << "\n";
    headerWritten_ = true;
}

void Profile1DWriter::write(int step,
                            double time_code,
                            double time_SI_s,
                            const QuadtreeMesh& mesh,
                            const UnitsContext* units,
                            const NormalizationContext* norm,
                            const CellField<double>& phi,
                            const CellField<Vector2D>& Ecell,
                            const CellField<double>& rho_tilde,
                            const std::vector<ParticleType>& species,
                            const std::vector<double>& macroWeight,
                            const BuildCacheFn& buildCacheForSpecies,
                            const std::vector<double>* alphaByCellIndex,
                            bool writeSIOutputs)
{
    if (!cfg_.enabled) return;
    if (!out_.is_open()) return;
    if (!headerWritten_) writeHeader(species, writeSIOutputs);

    const int nc = mesh.cellIndexNum();
    if (nc <= 0) return;
    const int nb = std::max(1, cfg_.nbinsX);

    // Deposit number amount per species onto current mesh to construct densities
    std::vector<CellField<double>> Namount;
    Namount.reserve(species.size());
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        CellField<double> tmp;
        tmp.initialize(mesh, 0.0);
        Namount.push_back(std::move(tmp));
    }

    for (std::size_t si = 0; si < species.size(); ++si)
    {
        CouplingWeightsCache cache;
        buildCacheForSpecies(si, cache);

        const auto& pts = species[si].particle;
        const std::size_t np = pts.size();
        if (np == 0) continue;

        std::vector<double> wp(np, 0.0);
        const double wDefault = (si < macroWeight.size() ? macroWeight[si] : 1.0);
        for (std::size_t p = 0; p < np; ++p)
        {
            const double w = pts[p].w;
            wp[p] = (w > 0.0 ? w : wDefault);
        }
        depositCSR(cache, wp, Namount[si]);
    }

    // Optional electron energy amounts
    const std::size_t electronIndex = [&]() -> std::size_t {
        for (std::size_t si = 0; si < species.size(); ++si)
            if (species[si].charge() < 0.0) return si;
        return static_cast<std::size_t>(-1);
    }();

    CellField<double> Eamount;
    Eamount.initialize(mesh, 0.0);
    bool haveElectronEnergy = false;
    if (electronIndex != static_cast<std::size_t>(-1) &&
        electronIndex < species.size() &&
        units != nullptr)
    {
        const auto& e = species[electronIndex];
        CouplingWeightsCache cache;
        buildCacheForSpecies(electronIndex, cache);

        const auto& pts = e.particle;
        const std::size_t np = pts.size();
        if (np > 0)
        {
            std::vector<double> wpE(np, 0.0);
            const double wDefault = (electronIndex < macroWeight.size() ? macroWeight[electronIndex] : 1.0);
            const double v0 = units->scales().v0_mps;
            const double qe = PhysicalConstants::qe_C;
            const double m = e.mass();
            std::size_t nBad = 0;
            for (std::size_t p = 0; p < np; ++p)
            {
                double w = (pts[p].w > 0.0 ? pts[p].w : wDefault);
                if (!std::isfinite(w) || !(w > 0.0))
                    w = wDefault;
                const double vx = pts[p].velocity.x * v0;
                const double vy = pts[p].velocity.y * v0;
                const double vz = pts[p].vz * v0;
                if (!(std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) && std::isfinite(v0) && std::isfinite(m) && (m > 0.0) && std::isfinite(w) && (w > 0.0)))
                {
                    nBad += 1;
                    continue;
                }
                const double ke_J = 0.5 * m * (vx * vx + vy * vy + vz * vz);
                const double ke_eV = ke_J / qe;
                if (!std::isfinite(ke_eV))
                {
                    nBad += 1;
                    continue;
                }
                wpE[p] = w * ke_eV;
            }
            depositCSR(cache, wpE, Eamount);
            haveElectronEnergy = true;
        }
    }

    const double x0 = mesh.minCoord.x;
    const double x1 = mesh.maxCoord.x;
    const double L = std::max(1e-300, x1 - x0);
    const double dx = L / static_cast<double>(nb);

    std::vector<double> sumVolCode(nb, 0.0);
    std::vector<double> sumVolSI(nb, 0.0);

    std::vector<double> sumPhiCodeVol(nb, 0.0);
    std::vector<double> sumExCodeVol(nb, 0.0);
    std::vector<double> sumEyCodeVol(nb, 0.0);
    std::vector<double> sumRhoCodeVol(nb, 0.0);

    std::vector<std::vector<double>> sumN(nb, std::vector<double>(species.size(), 0.0));
    std::vector<double> sumEe(nb, 0.0);

    for (int c = 0; c < nc; ++c)
    {
        const int cLbl = mesh.cellIndexToLabel[c];
        const Cell& cell = mesh.cells(cLbl);
        const double xc = cell.centerCoord().x;
        const double sx = cell.size().x;
        const double sy = cell.size().y;
        if (!(sx > 0.0) || !(sy > 0.0))
            continue;

        const double xa = xc - 0.5 * sx;
        const double xb = xc + 0.5 * sx;
        const int b0_raw = static_cast<int>(std::floor((xa - x0) / dx));
        const int b1_raw = static_cast<int>(std::floor((xb - x0) / dx));
        const int b0 = std::max(0, std::min(nb - 1, b0_raw));
        const int b1 = std::max(0, std::min(nb - 1, b1_raw));

        double alpha = 1.0;
        if (alphaByCellIndex && c < static_cast<int>(alphaByCellIndex->size()))
            alpha = std::max(0.0, (*alphaByCellIndex)[c]);
        if (!(alpha > 0.0))
            continue;

        // Split this cell's contribution across all intersecting x-bins.
        // Assume cell-average fields and uniform particle amount density in the cell.
        for (int b = b0; b <= b1; ++b)
        {
            const double xl = x0 + static_cast<double>(b) * dx;
            const double xr = xl + dx;
            const double overlap = std::min(xb, xr) - std::max(xa, xl);
            if (!(overlap > 0.0))
                continue;

            const double fracX = overlap / sx;
            const double vol_code = overlap * sy * alpha;
            if (!(vol_code > 0.0))
                continue;
            const double vol_SI = (norm ? norm->VolumeCodeToSI(vol_code) : vol_code);

            sumVolCode[b] += vol_code;
            sumVolSI[b] += vol_SI;

            if (cfg_.includePhi && c < static_cast<int>(phi.size()))
                sumPhiCodeVol[b] += phi.value[c] * vol_code;

            if (cfg_.includeE && c < static_cast<int>(Ecell.size()))
            {
                sumExCodeVol[b] += Ecell.value[c].x * vol_code;
                sumEyCodeVol[b] += Ecell.value[c].y * vol_code;
            }

            if (cfg_.includeRho && c < static_cast<int>(rho_tilde.size()))
                sumRhoCodeVol[b] += rho_tilde.value[c] * vol_code;

            for (std::size_t si = 0; si < species.size(); ++si)
            {
                if (c < static_cast<int>(Namount[si].size()))
                    sumN[b][si] += Namount[si].value[c] * fracX;
            }

            if (haveElectronEnergy && c < static_cast<int>(Eamount.size()))
                sumEe[b] += Eamount.value[c] * fracX;
        }
    }

    out_ << std::setprecision(17);
    for (int b = 0; b < nb; ++b)
    {
        const double Vc = sumVolCode[b];
        const double Vs = sumVolSI[b];
        const bool haveVol = (Vc > 0.0) && (Vs > 0.0);

        const double x_code = x0 + (static_cast<double>(b) + 0.5) * dx;
        const double x_SI = (units && writeSIOutputs) ? units->LengthCodeToSI(x_code) : 0.0;

        out_ << step << ',' << time_code;
        if (writeSIOutputs) out_ << ',' << time_SI_s;

        out_ << ',' << b << ',' << x_code;
        if (writeSIOutputs) out_ << ',' << x_SI;

        if (cfg_.includePhi)
        {
            const double phi_code = haveVol ? (sumPhiCodeVol[b] / Vc) : 0.0;
            out_ << ',' << phi_code;
            if (writeSIOutputs)
                out_ << ',' << (units ? units->PotentialCodeToSI(phi_code) : phi_code);
        }
        if (cfg_.includeE)
        {
            const double Ex_code = haveVol ? (sumExCodeVol[b] / Vc) : 0.0;
            const double Ey_code = haveVol ? (sumEyCodeVol[b] / Vc) : 0.0;
            out_ << ',' << Ex_code << ',' << Ey_code;
            if (writeSIOutputs)
            {
                out_ << ',' << (units ? units->ElectricFieldCodeToSI(Ex_code) : Ex_code);
                out_ << ',' << (units ? units->ElectricFieldCodeToSI(Ey_code) : Ey_code);
            }
        }
        if (cfg_.includeRho)
        {
            const double rho_code = haveVol ? (sumRhoCodeVol[b] / Vc) : 0.0;
            out_ << ',' << rho_code;
            if (writeSIOutputs)
                out_ << ',' << (norm ? norm->RhoTildeToSI(rho_code) : rho_code);
        }

        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const double n_SI = haveVol ? (sumN[b][si] / Vs) : 0.0;
            out_ << ',' << n_SI;
        }

        const double Ne_bin = (electronIndex != static_cast<std::size_t>(-1) ? sumN[b][electronIndex] : 0.0);
        const double Emean_e = (haveElectronEnergy && Ne_bin > 0.0) ? (sumEe[b] / Ne_bin) : 0.0;
        out_ << ',' << (std::isfinite(Emean_e) ? Emean_e : 0.0);
        out_ << ',' << (std::isfinite(Emean_e) ? (2.0 / 3.0) * Emean_e : 0.0);

        out_ << "\n";
    }
}

} // namespace Diagnostics
} // namespace ExplicitPIC
