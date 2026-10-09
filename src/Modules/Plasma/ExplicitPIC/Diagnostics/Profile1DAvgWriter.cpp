#include "Modules/Plasma/ExplicitPIC/Diagnostics/Profile1DAvgWriter.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>

#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Units/UnitsContext.h"

#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"

namespace ExplicitPIC
{
namespace Diagnostics
{

namespace
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

static int findElectronIndex(const std::vector<ParticleType>& species)
{
    auto lowerCopy = [](std::string s) {
        for (char& c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const std::string nm = lowerCopy(species[si].name());
        if (nm == "e" || nm == "electron" || nm.find("electron") != std::string::npos)
            return static_cast<int>(si);
    }
    for (std::size_t si = 0; si < species.size(); ++si)
        if (species[si].charge() < 0.0) return static_cast<int>(si);
    return -1;
}

} // namespace

void Profile1DAvgWriter::ensureParentDir(const std::string& filePath)
{
    try
    {
        const std::filesystem::path p(filePath);
        if (p.has_parent_path())
            std::filesystem::create_directories(p.parent_path());
    }
    catch (...)
    {
        // best effort
    }
}

void Profile1DAvgWriter::reset(std::size_t nBinsX)
{
    nBinsX_ = std::max<std::size_t>(1, nBinsX);
    acc_.reset(nBinsX_, static_cast<std::size_t>(NUM_FIELDS),
               static_cast<std::size_t>(std::max(0, cfg_.avgWindowSteps)));
}

void Profile1DAvgWriter::accumulate(const QuadtreeMesh& mesh,
                                    const UnitsContext* units,
                                    const NormalizationContext* norm,
                                    const CellField<double>& phi,
                                    const CellField<Vector2D>& Ecell,
                                    const std::vector<ParticleType>& species,
                                    const std::vector<double>& macroWeight,
                                    const CacheRefFn& cacheForSpecies,
                                    const std::vector<double>* alphaByCellIndex)
{
    if (!enabled()) return;
    if (nBinsX_ == 0) reset(static_cast<std::size_t>(std::max(1, 128)));
    if (!acc_.enabled()) return;

    const int nc = mesh.cellIndexNum();
    if (nc <= 0) return;

    // Deposit number amount per species onto current mesh to construct densities.
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
        const CouplingWeightsCache& cache = cacheForSpecies(si);

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

    const int eIndex = findElectronIndex(species);

    // Electron energy amount (eV) deposited to cells.
    CellField<double> Eamount;
    Eamount.initialize(mesh, 0.0);
    bool haveElectronEnergy = false;
    if (eIndex >= 0 && eIndex < static_cast<int>(species.size()) && units != nullptr)
    {
        const auto& e = species[static_cast<std::size_t>(eIndex)];
        const CouplingWeightsCache& cache = cacheForSpecies(static_cast<std::size_t>(eIndex));

        const auto& pts = e.particle;
        const std::size_t np = pts.size();
        if (np > 0)
        {
            std::vector<double> wpE(np, 0.0);
            const double wDefault = (static_cast<std::size_t>(eIndex) < macroWeight.size() ? macroWeight[static_cast<std::size_t>(eIndex)] : 1.0);
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

            if (nBad > 0)
            {
                // Keep this low-noise: only emit in debug builds / when caller enables it.
                // (We don't have console verbosity here, so no stdout by default.)
            }
        }
    }

    const double x0 = mesh.minCoord.x;
    const double x1 = mesh.maxCoord.x;
    const double L = std::max(1e-300, x1 - x0);
    const double dx = L / static_cast<double>(nBinsX_);

    std::vector<double> sumVolCode(nBinsX_, 0.0);
    std::vector<double> sumVolSI(nBinsX_, 0.0);
    std::vector<double> sumPhiCodeVol(nBinsX_, 0.0);
    std::vector<double> sumExCodeVol(nBinsX_, 0.0);
    std::vector<double> sumNeAmount(nBinsX_, 0.0);
    std::vector<double> sumNiAmount(nBinsX_, 0.0);
    std::vector<double> sumEe(nBinsX_, 0.0);

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
        const int nb = static_cast<int>(nBinsX_);
        const int b0_raw = static_cast<int>(std::floor((xa - x0) / dx));
        const int b1_raw = static_cast<int>(std::floor((xb - x0) / dx));
        const int b0 = std::max(0, std::min(nb - 1, b0_raw));
        const int b1 = std::max(0, std::min(nb - 1, b1_raw));

        double alpha = 1.0;
        if (alphaByCellIndex && c < static_cast<int>(alphaByCellIndex->size()))
            alpha = std::max(0.0, (*alphaByCellIndex)[c]);
        if (!(alpha > 0.0))
            continue;

        // Total ion (positive) number amount in this cell.
        double Ni_cell = 0.0;
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            if (!(species[si].charge() > 0.0)) continue;
            if (c < static_cast<int>(Namount[si].size()))
                Ni_cell += Namount[si].value[c];
        }

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

            sumVolCode[static_cast<std::size_t>(b)] += vol_code;
            sumVolSI[static_cast<std::size_t>(b)] += vol_SI;

            if (c < static_cast<int>(phi.size()))
                sumPhiCodeVol[static_cast<std::size_t>(b)] += phi.value[c] * vol_code;

            // 1D CCP benchmark convention: E(x) is taken as the x-component of the
            // cell-centered electric field, y-averaged within each x-bin.
            if (c < static_cast<int>(Ecell.size()))
                sumExCodeVol[static_cast<std::size_t>(b)] += Ecell.value[c].x * vol_code;

            if (eIndex >= 0 && eIndex < static_cast<int>(Namount.size()) && c < static_cast<int>(Namount[static_cast<std::size_t>(eIndex)].size()))
                sumNeAmount[static_cast<std::size_t>(b)] += Namount[static_cast<std::size_t>(eIndex)].value[c] * fracX;

            sumNiAmount[static_cast<std::size_t>(b)] += Ni_cell * fracX;

            if (haveElectronEnergy && c < static_cast<int>(Eamount.size()))
                sumEe[static_cast<std::size_t>(b)] += Eamount.value[c] * fracX;
        }
    }

    // Assemble sample vector (field-major): [NE, NI, TE, PHI, E] each length nBinsX_.
    std::vector<double> sample(static_cast<std::size_t>(NUM_FIELDS) * nBinsX_,
                               std::numeric_limits<double>::quiet_NaN());

    for (std::size_t b = 0; b < nBinsX_; ++b)
    {
        const double Vc = sumVolCode[b];
        const double Vs = sumVolSI[b];
        const bool haveVol = (Vc > 0.0) && (Vs > 0.0);

        const double ne = haveVol ? (sumNeAmount[b] / Vs) : 0.0;
        const double ni = haveVol ? (sumNiAmount[b] / Vs) : 0.0;

        const double phi_code = haveVol ? (sumPhiCodeVol[b] / Vc) : 0.0;
        const double phi_V = (units ? units->PotentialCodeToSI(phi_code) : phi_code);

        const double Ex_code = haveVol ? (sumExCodeVol[b] / Vc) : 0.0;
        const double Ex_Vpm = (units ? units->ElectricFieldCodeToSI(Ex_code) : Ex_code);

        double Te_eV = 0.0;
        if (haveElectronEnergy && sumNeAmount[b] > 0.0)
        {
            const double Emean = sumEe[b] / sumNeAmount[b];
            Te_eV = (2.0 / 3.0) * Emean;
        }

        sample[static_cast<std::size_t>(NE) * nBinsX_ + b] = ne;
        sample[static_cast<std::size_t>(NI) * nBinsX_ + b] = ni;
        sample[static_cast<std::size_t>(TE) * nBinsX_ + b] = Te_eV;
        sample[static_cast<std::size_t>(PHI) * nBinsX_ + b] = phi_V;
        sample[static_cast<std::size_t>(E) * nBinsX_ + b] = Ex_Vpm;
    }

    acc_.pushSample(sample);
}

void Profile1DAvgWriter::writeSnapshotCSV(const QuadtreeMesh& mesh, const UnitsContext* units) const
{
    if (!enabled()) return;
    if (!acc_.enabled()) return;
    if (nBinsX_ == 0) return;

    ensureParentDir(cfg_.path);
    std::ofstream out(cfg_.path, std::ios::out | std::ios::trunc);
    if (!out) return;

    std::vector<double> mean;
    std::vector<double> stdv;
    acc_.meanStd(mean, stdv);

    const double x0 = mesh.minCoord.x;
    const double x1 = mesh.maxCoord.x;
    const double L = std::max(1e-300, x1 - x0);
    const double dx = L / static_cast<double>(nBinsX_);

    // Column names are chosen to be directly consumed by postprocess/compare_turner.py
    // (PB8.2-4): phi_mean_V, E_mean_Vpm, ne_mean_m-3 ...
    out << "bin,x_m,"
        << "phi_mean_V,phi_std_V,"
        << "E_mean_Vpm,E_std_Vpm,"
        << "ne_mean_m-3,ne_std_m-3,"
        << "ni_mean_m-3,ni_std_m-3,"
        << "Te_mean_eV,Te_std_eV\n";

    out << std::setprecision(17);
    auto finiteOrZero = [](double v) -> double { return std::isfinite(v) ? v : 0.0; };
    for (std::size_t b = 0; b < nBinsX_; ++b)
    {
        const double x_code = x0 + (static_cast<double>(b) + 0.5) * dx;
        const double x_m = (units ? units->LengthCodeToSI(x_code) : x_code);

        const double ph_m = mean[static_cast<std::size_t>(PHI) * nBinsX_ + b];
        const double ph_s = stdv[static_cast<std::size_t>(PHI) * nBinsX_ + b];
        const double ex_m = mean[static_cast<std::size_t>(E) * nBinsX_ + b];
        const double ex_s = stdv[static_cast<std::size_t>(E) * nBinsX_ + b];
        const double ne_m = mean[static_cast<std::size_t>(NE) * nBinsX_ + b];
        const double ne_s = stdv[static_cast<std::size_t>(NE) * nBinsX_ + b];
        const double ni_m = mean[static_cast<std::size_t>(NI) * nBinsX_ + b];
        const double ni_s = stdv[static_cast<std::size_t>(NI) * nBinsX_ + b];
        const double te_m = mean[static_cast<std::size_t>(TE) * nBinsX_ + b];
        const double te_s = stdv[static_cast<std::size_t>(TE) * nBinsX_ + b];

        out << b << ',' << x_m << ','
            << finiteOrZero(ph_m) << ',' << finiteOrZero(ph_s) << ','
            << finiteOrZero(ex_m) << ',' << finiteOrZero(ex_s) << ','
            << finiteOrZero(ne_m) << ',' << finiteOrZero(ne_s) << ','
            << finiteOrZero(ni_m) << ',' << finiteOrZero(ni_s) << ','
            << finiteOrZero(te_m) << ',' << finiteOrZero(te_s) << '\n';
    }
}

} // namespace Diagnostics
} // namespace ExplicitPIC
