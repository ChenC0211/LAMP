#include "Modules/Plasma/IO/PicFieldFrame.h"
#include "Modules/Plasma/Diagnostics/PicParticleFields.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/ConsoleControl.h"
#include "Core/Base/Interpolation.h"
#include "Core/Geometry/LevelSet.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include <cctype>
#include <iostream>
#include <stdexcept>

namespace
{
int findElectronIndex(const std::vector<ParticleType>& species)
{
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        std::string name = species[si].name();
        for (char& c : name)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (name == "e" || name == "electron" || name.find("electron") != std::string::npos)
            return static_cast<int>(si);
    }
    for (std::size_t si = 0; si < species.size(); ++si)
        if (species[si].charge() < 0.0)
            return static_cast<int>(si);
    return -1;
}
}

PicFieldFrame::PicFieldFrame(const QuadtreeMesh& mesh, const PicFieldSelection& selection,
    std::vector<unsigned char> mask, const CellField<double>& potential, bool cellConstant)
    : validCells(std::move(mask)), mesh_(mesh), selection_(selection),
      cellConstant_(selection.mode() == PicFieldMode::Implicit && cellConstant)
{
    if (selection_.mode() == PicFieldMode::Implicit)
        markedCells_ = PicOutputReconstruction::markZeroPotentialCells(mesh_, potential, validCells);
}

void PicFieldFrame::reportClamped(std::uint64_t count, const std::string& name)
{
    if (ConsoleControl::enabled(ConsoleControl::Level::kDebug))
        std::cout << "[ExplicitPIC][output] NOTE: clamped " << count
                  << " negative values in field '" << name << "'\n";
}

void PicFieldFrame::addParticles(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, const NormalizationContext& norm,
    const std::vector<double>* alpha, const BuildCache& buildCache,
    bool fastDensity, const std::vector<CellField<double>>* depositedAmounts)
{
    const bool explicitMode = selection_.mode() == PicFieldMode::Explicit;
    const auto& selected = selection_.fields();
    std::vector<CouplingWeightsCache> outputSpeciesCache(explicitMode ? 0 : species.size());
    std::vector<unsigned char> outputSpeciesCacheBuilt(outputSpeciesCache.size(), 0u);
    auto getSpeciesCache = [&](std::size_t si) -> const CouplingWeightsCache& {
        if (si >= outputSpeciesCache.size())
            throw std::runtime_error("[ImplicitPIC][output] species cache index out of range");
        if (outputSpeciesCacheBuilt[si] == 0u)
        {
            buildCache(si, outputSpeciesCache[si]);
            outputSpeciesCacheBuilt[si] = 1u;
        }
        return outputSpeciesCache[si];
    };

    if (selected.n)
    {
        auto fieldKeys = selection_.keys(species.size());
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const std::string key = fieldKeys.next(species[si].name());
            const int n = static_cast<int>(species[si].particle.size());
            CellField<double> Namount(mesh_, 0.0);
            const double weight = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
            if (!explicitMode && fastDensity && depositedAmounts &&
                si < depositedAmounts->size() && (*depositedAmounts)[si].size() == mesh_.cellIndexNum())
                Namount = (*depositedAmounts)[si];
            else if (n > 0)
            {
                if (explicitMode)
                {
                    CouplingWeightsCache cache;
                    buildCache(si, cache);
                    PicParticleFields::depositNumberAmount(species[si], weight, cache, Namount);
                }
                else
                    PicParticleFields::depositNumberAmount(species[si], weight, getSpeciesCache(si), Namount);
            }
            CellField<double> nSI(mesh_, 0.0);
            PicParticleFields::numberDensitySI(mesh_, norm, Namount, alpha, nSI);
            addScalar("n_" + key + "_SI_1pm3", nSI, [](double v) { return v; }, true);
        }
    }

    if (selected.te || selected.meanEnergy)
    {
        const int eIndex = findElectronIndex(species);
        if (eIndex >= 0)
        {
            const std::size_t si = static_cast<std::size_t>(eIndex);
            const int n = static_cast<int>(species[si].particle.size());
            if (n > 0)
            {
                // Explicit used a fresh cache before creating the output fields.
                CouplingWeightsCache explicitCache;
                if (explicitMode)
                    buildCache(si, explicitCache);
                CellField<double> Emean_eV(mesh_, 0.0);
                CellField<double> Te_eV(mesh_, 0.0);
                const double weight = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
                if (explicitMode)
                {
                    const auto nBad = PicParticleFields::electronEnergyCSR(mesh_, species[si], weight, norm,
                        explicitCache, PicParticleFields::ElectronWeightPolicy::FinitePositive, Emean_eV, Te_eV);
                    if (nBad > 0 && ConsoleControl::enabled(ConsoleControl::Level::kDebug))
                        std::cerr << "[ExplicitPIC][output] WARNING: skipped " << nBad
                                  << " non-finite electron particles when computing Te/mean-energy fields." << std::endl;
                }
                else if (fastDensity)
                    PicParticleFields::electronEnergyHostCells(mesh_, species[si], weight, norm, Emean_eV, Te_eV);
                else
                    PicParticleFields::electronEnergyCSR(mesh_, species[si], weight, norm, getSpeciesCache(si),
                        PicParticleFields::ElectronWeightPolicy::Positive, Emean_eV, Te_eV);
                if (selected.meanEnergy)
                    addScalar("electron_mean_energy_eV", Emean_eV, [](double v) { return v; });
                if (selected.te)
                    addScalar("Te_e_eV", Te_eV, [](double v) { return v; });
            }
        }
    }

    if (selected.weightMean || selected.weightSum || selected.weightMax || selected.macroCount)
    {
        CellField<double> Wsum(mesh_, 0.0), Wmean(mesh_, 0.0), Wmax(mesh_, 0.0), Nmacro(mesh_, 0.0);
        PicParticleFields::particleWeights(mesh_, species, macroWeight, Wsum, Wmean, Wmax, Nmacro);
        if (selected.weightMean)
            addScalar("w_mean", Wmean, [](double v) { return v; });
        if (selected.weightSum)
            addScalar("w_sum", Wsum, [](double v) { return v; });
        if (selected.weightMax)
            addScalar("w_max", Wmax, [](double v) { return v; });
        if (selected.macroCount)
            addScalar("n_macro", Nmacro, [](double v) { return v; }, true);
    }
}

void PicFieldFrame::addElectrostatic(const CellField<double>& phi, const CellField<double>& rhoStar,
    const UnitsContext& units, const NormalizationContext& norm)
{
    const auto& selected = selection_.fields();
    if (selected.phi)
        addScalar("phi_SI_V", phi, [&](double v) { return units.PotentialCodeToSI(v); });
    if (selected.rho)
        addScalar("rho_SI_Cpm3", rhoStar, [&](double v) { return norm.RhoTildeToSI(v); });
}

void PicFieldFrame::addElectric(const CellField<Vector2D>& Ecell, const UnitsContext& units)
{
    const auto& selected = selection_.fields();
    // Electric field components.
    if (selected.ex || selected.ey)
    {
        CellField<double> Ex(mesh_, 0.0);
        CellField<double> Ey(mesh_, 0.0);
        for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
        {
            Ex.value[ci] = Ecell.value[ci].x;
            Ey.value[ci] = Ecell.value[ci].y;
        }
        if (selected.ex)
            addScalar("Ex_SI_Vpm", Ex, [&](double v) { return units.ElectricFieldCodeToSI(v); });
        if (selected.ey)
            addScalar("Ey_SI_Vpm", Ey, [&](double v) { return units.ElectricFieldCodeToSI(v); });
    }
}

void PicFieldFrame::addLevelSets(const UnitsContext& units, const VertField<double>* composite,
    const std::vector<LevelSet>* levelSets, const std::vector<int>& indices,
    const std::vector<std::string>& levelSetNames)
{
    const auto& selected = selection_.fields();
    if (selected.psiAll && levelSets && !indices.empty())
    {
        for (int lsI : indices)
        {
            if (lsI < 0 || lsI >= static_cast<int>((*levelSets).size()))
                continue;

            CellField<double> psiLevelSet(mesh_, 0.0);
            for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
            {
                const int cLbl = mesh_.cellIndexToLabel[ci];
                const Cell& c = mesh_.cells(cLbl);
                psiLevelSet.value[ci] = Interpolation::biLinearInterpolation(
                    mesh_, (*levelSets)[lsI], c, c.centerCoord());
            }

            std::string lsName;
            if (lsI >= 0 && lsI < static_cast<int>(levelSetNames.size()))
                lsName = selection_.key(levelSetNames[lsI]);
            if (lsName.empty())
                lsName = std::string("ls_") + std::to_string(lsI + 1);

            addScalar("psi_" + lsName + "_SI_m", psiLevelSet,
                          [&](double v) { return units.LengthCodeToSI(v); });
        }
    }

    if (selected.psiComposite || selected.psi)
    {
        CellField<double> psiC(mesh_, 0.0);
        if (composite)
        {
            for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
            {
                const int cLbl = mesh_.cellIndexToLabel[ci];
                const Cell& c = mesh_.cells(cLbl);
                psiC.value[ci] = Interpolation::biLinearInterpolation(mesh_, *composite, c, c.centerCoord());
            }
        }
        const std::string compositeFieldName = selected.psi ? "psi_SI_m" : "psi_composite_SI_m";
        addScalar(compositeFieldName, psiC, [&](double v) { return units.LengthCodeToSI(v); });
    }
}

void PicFieldFrame::addSusceptibility(const FaceSegmentGraph& faceGraph,
    const std::vector<double>& muCell_Fpm_now, const std::vector<double>& kappaByCell_now,
    const std::vector<double>& kappaBySeg_now)
{
    const auto& selected = selection_.fields();
    CellField<double> muCell_SI(mesh_, 0.0);
    if (selected.muCell)
    {
        for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
        {
            const std::size_t idx = static_cast<std::size_t>(ci);
            const double muC = (idx < muCell_Fpm_now.size()) ? muCell_Fpm_now[idx] : 0.0;
            muCell_SI.value[ci] = muC;
        }
    }

    CellField<double> muFace_SI(mesh_, 0.0);
    if (selected.muFace)
    {
        // Build a cell-centered muFace proxy by averaging incident face-segment kappa.
        std::vector<double> kappaFaceSum(static_cast<std::size_t>(mesh_.cellIndexNum()), 0.0);
        std::vector<int> kappaFaceCount(static_cast<std::size_t>(mesh_.cellIndexNum()), 0);
        for (std::size_t sid = 0; sid < faceGraph.segments().size(); ++sid)
        {
            const FaceSegment& s = faceGraph.segments()[sid];
            const double kf = (sid < kappaBySeg_now.size()) ? kappaBySeg_now[sid] : 1.0;
            if (s.K >= 0)
            {
                const int kIdx = mesh_.cells(s.K).index();
                if (kIdx >= 0)
                {
                    kappaFaceSum[static_cast<std::size_t>(kIdx)] += kf;
                    kappaFaceCount[static_cast<std::size_t>(kIdx)] += 1;
                }
            }
            if (s.L >= 0)
            {
                const int lIdx = mesh_.cells(s.L).index();
                if (lIdx >= 0)
                {
                    kappaFaceSum[static_cast<std::size_t>(lIdx)] += kf;
                    kappaFaceCount[static_cast<std::size_t>(lIdx)] += 1;
                }
            }
        }

        for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
        {
            const std::size_t idx = static_cast<std::size_t>(ci);
            double kFace = (idx < kappaByCell_now.size()) ? kappaByCell_now[idx] : 1.0;
            if (kappaFaceCount[idx] > 0)
                kFace = kappaFaceSum[idx] / static_cast<double>(kappaFaceCount[idx]);
            muFace_SI.value[ci] = (kFace - 1.0) * PhysicalConstants::epsilon0_Fpm;
        }
    }

    // Optional mu outputs (kept for debugging).
    if (selected.muCell)
        addScalar("muCell_SI_Fpm", muCell_SI, [&](double v) { return v; });
    if (selected.muFace)
        addScalar("muFaceProxy_SI_Fpm", muFace_SI, [&](double v) { return v; });
}

void PicFieldFrame::finish(bool periodicX, bool periodicY)
{
    if (selection_.mode() != PicFieldMode::Implicit)
        return;
    if (cellConstant_)
        selection_.order(names, cells);
    else
    {
        selection_.order(names, vertices);
        PicOutputReconstruction::canonicalizePeriodicBoundaryVertices(mesh_, periodicX, periodicY, vertices);
    }
}
