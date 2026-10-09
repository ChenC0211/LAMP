#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

namespace
{
using Clock = std::chrono::steady_clock;
double secondsBetween(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(end - begin).count();
}
}
PicAmrGeometry::PicAmrGeometry(Mode mode, Options options)
    : mode_(mode), options_(std::move(options)) {}

void PicAmrGeometry::captureBeforeAdapt(const Geometry* geometry)
{
    report_ = {};
    preAdaptDiagnostics_.clear();
    if (mode_ == Mode::Implicit && geometry)
        preAdaptDiagnostics_ = geometry->levelSetUpdateDiagnostics();
}

void PicAmrGeometry::refresh(QuadtreeMesh& mesh, Geometry* iGeometry,
    PicEmbeddedGeometry& embeddedGeometry, const MovingEBController* movingEb,
    const Adaptation& adaptation)
{
    if (mode_ == Mode::Explicit)
    {
        if (iGeometry && (embeddedGeometry.ready() || adaptation.refinementActive))
        {
            iGeometry->reinitializeLevelSets(mesh, 0.0);
            if (embeddedGeometry.ready())
            {
                iGeometry->checkCutInfo(mesh, /*isInit=*/true);
                embeddedGeometry.rebuildComposite(mesh, iGeometry);
            }
            if (adaptation.refinementActive && adaptation.phaseActive)
                iGeometry->markPhases(mesh);
        }
        return;
    }
    std::vector<Geometry::LevelSetUpdateDiagnostics> postAdaptDiagnostics;
    if (iGeometry && options_.amrReinitHybridSensitiveInactive)
    {
        iGeometry->refreshLevelSetScalarDiagnostics(mesh);
        postAdaptDiagnostics = iGeometry->levelSetUpdateDiagnostics();
    }
    std::vector<int> amrActiveMovingLevelSets;
    if (iGeometry && (adaptation.embeddedEnabled || adaptation.refinementActive))
    {
        const bool forceFullGeomRefresh =
            ((movingEb != nullptr) && !options_.amrGeomRefreshAllowWithMovingEb) ||
            adaptation.movingUpdated;
        if (!forceFullGeomRefresh)
        {
            if (options_.mode == RefreshMode::Skip)
                report_.geometryReinitSkipped = true;
            else if (options_.mode == RefreshMode::Adaptive)
                report_.geometryReinitSkipped = (adaptation.cellCountChange <= options_.amrGeomRefreshChangeRelTh);
        }
        const auto tGeometryRefresh0 = Clock::now();
        if (!report_.geometryReinitSkipped)
        {
            if (movingEb != nullptr)
                report_.amrRawActiveMovingLevelSets = movingEb->activeMovingLevelSetIndices(adaptation.step);
            const bool hasExactMovingEbGeometry =
                (movingEb != nullptr) &&
                movingEb->linesInCurrentGeometryFrame() &&
                !movingEb->lastUpdateUsedAdvectedPath();
            // If MovingEB already has explicit line geometry in the current
            // geometry frame, any AMR remesh on the same frozen-geometry
            // interval should rebuild level sets directly from those lines
            // instead of reinitializing from the sampled moving LS field.
            const bool useExactMovingEbRefreshOnly = hasExactMovingEbGeometry;
            if (useExactMovingEbRefreshOnly)
            {
                const auto tGeomInitLs0 = Clock::now();
                iGeometry->initializeLevelSets(mesh);
                report_.geomInitializeLevelSetsS = secondsBetween(tGeomInitLs0, Clock::now());

                const auto tGeomRefreshDerived0 = Clock::now();
                iGeometry->refreshLevelSetDerivedData(mesh, 0.0, /*reinitializeMovingLevelSets=*/false);
                report_.geomRefreshDerivedS = secondsBetween(tGeomRefreshDerived0, Clock::now());
                report_.geometryRefreshUsedExactMovingEb = true;
                selectiveStreak_ = 0;
            }
            else
            {
                bool useSelectiveMovingReinit = false;
                if (options_.amrReinitActiveMovingOnly && movingEb != nullptr)
                {
                    amrActiveMovingLevelSets = report_.amrRawActiveMovingLevelSets;
                    std::vector<unsigned char> activeMask(
                        iGeometry->levelSets.size(), static_cast<unsigned char>(0));
                    for (const int idx : report_.amrRawActiveMovingLevelSets)
                    {
                        if (idx >= 0 && idx < static_cast<int>(activeMask.size()))
                            activeMask[static_cast<std::size_t>(idx)] = 1u;
                    }
                    if (options_.amrReinitHybridSensitiveInactive &&
                        postAdaptDiagnostics.size() == iGeometry->levelSets.size())
                    {
                        for (int lsI = 0; lsI < static_cast<int>(iGeometry->levelSets.size()); ++lsI)
                        {
                            const auto& ls = iGeometry->levelSets[lsI];
                            if (!ls.moving())
                                continue;
                            if (activeMask[static_cast<std::size_t>(lsI)] != 0u)
                                continue;
                            const auto& preDiag = preAdaptDiagnostics_[static_cast<std::size_t>(lsI)];
                            const auto& diag = postAdaptDiagnostics[static_cast<std::size_t>(lsI)];
                            const int nearZeroDelta =
                                std::abs(diag.nearZeroCount - preDiag.nearZeroCount);
                            const int narrowBandDelta =
                                std::abs(diag.narrowBandCount - preDiag.narrowBandCount);
                            const bool remapSensitiveMagnitude =
                                (diag.nearZeroCount >= options_.amrReinitHybridNearZeroMin) ||
                                (diag.narrowBandCount >= options_.amrReinitHybridNarrowBandMin);
                            const bool remapSensitiveDelta =
                                (nearZeroDelta >= options_.amrReinitHybridNearZeroDeltaMin) ||
                                (narrowBandDelta >= options_.amrReinitHybridNarrowBandDeltaMin);
                            const bool remapSensitive =
                                ((options_.amrReinitHybridUseDelta ? remapSensitiveDelta : remapSensitiveMagnitude)) ||
                                (options_.amrReinitHybridIncludeMaxIter &&
                                 diag.reinitIterations >= options_.amrReinitHybridMaxIterThreshold);
                            if (!remapSensitive)
                                continue;
                            report_.amrHybridSensitiveInactiveLevelSets.push_back(lsI);
                        }
                        amrActiveMovingLevelSets.insert(
                            amrActiveMovingLevelSets.end(),
                            report_.amrHybridSensitiveInactiveLevelSets.begin(),
                            report_.amrHybridSensitiveInactiveLevelSets.end());
                        std::sort(amrActiveMovingLevelSets.begin(), amrActiveMovingLevelSets.end());
                        amrActiveMovingLevelSets.erase(
                            std::unique(amrActiveMovingLevelSets.begin(), amrActiveMovingLevelSets.end()),
                            amrActiveMovingLevelSets.end());
                    }
                    int movingLevelSetCount = 0;
                    for (const auto& ls : iGeometry->levelSets)
                    {
                        if (ls.moving())
                            ++movingLevelSetCount;
                    }
                    if (amrActiveMovingLevelSets.empty() ||
                        static_cast<int>(amrActiveMovingLevelSets.size()) >= movingLevelSetCount)
                    {
                        amrActiveMovingLevelSets.clear();
                    }
                    if (!amrActiveMovingLevelSets.empty() &&
                        options_.amrReinitActiveMovingMaxStreak > 0 &&
                        selectiveStreak_ >= options_.amrReinitActiveMovingMaxStreak)
                    {
                        amrActiveMovingLevelSets.clear();
                    }
                }
                useSelectiveMovingReinit = !amrActiveMovingLevelSets.empty();
                if (useSelectiveMovingReinit)
                    report_.amrSelectiveMovingLevelSets = amrActiveMovingLevelSets;
                if (options_.amrGeomRefreshVerbose && !report_.amrSelectiveMovingLevelSets.empty())
                {
                    std::cout << "[ImplicitPIC][AMR] selective moving-LS reinit count="
                              << report_.amrSelectiveMovingLevelSets.size() << std::endl;
                    if (!report_.amrHybridSensitiveInactiveLevelSets.empty())
                    {
                        std::cout << "[ImplicitPIC][AMR] hybrid inactive moving-LS added="
                                  << report_.amrHybridSensitiveInactiveLevelSets.size() << std::endl;
                    }
                }

                const auto tGeomPdeReinit0 = Clock::now();
                if (useSelectiveMovingReinit)
                    iGeometry->pdeReinitializeMovingLevelSets(mesh, 0.0, report_.amrSelectiveMovingLevelSets);
                else
                    iGeometry->pdeReinitializeMovingLevelSets(mesh, 0.0);
                report_.geomPdeReinitS = secondsBetween(tGeomPdeReinit0, Clock::now());
                selectiveStreak_ = useSelectiveMovingReinit
                                                     ? (selectiveStreak_ + 1)
                                                     : 0;

                const auto tGeomRefreshDerived0 = Clock::now();
                iGeometry->refreshLevelSetDerivedData(mesh, 0.0, /*reinitializeMovingLevelSets=*/false);
                report_.geomRefreshDerivedS = secondsBetween(tGeomRefreshDerived0, Clock::now());
            }
        }
        else
        {
            selectiveStreak_ = 0;
        }
        if (adaptation.embeddedEnabled && embeddedGeometry.ready())
        {
            const auto tGeomCheckCutInfo0 = Clock::now();
            iGeometry->checkCutInfo(mesh, /*isInit=*/true);
            report_.geomCheckCutInfoS = secondsBetween(tGeomCheckCutInfo0, Clock::now());

            const auto tGeomRebuildComposite0 = Clock::now();
            embeddedGeometry.rebuildComposite(mesh, iGeometry);
            report_.geomRebuildCompositeS = secondsBetween(tGeomRebuildComposite0, Clock::now());
        }
        if (adaptation.refinementActive && adaptation.phaseActive)
        {
            const auto tGeomMarkPhases0 = Clock::now();
            iGeometry->markPhases(mesh);
            report_.geomMarkPhasesS = secondsBetween(tGeomMarkPhases0, Clock::now());
        }
        report_.geometryRefreshS = secondsBetween(tGeometryRefresh0, Clock::now());

        if (options_.amrGeomRefreshVerbose)
        {
            std::cout << "[ImplicitPIC][AMR] geometry refresh "
                      << (report_.geometryReinitSkipped ? "SKIP_REINIT"
                          : (report_.geometryRefreshUsedExactMovingEb ? "REFRESH_ONLY_EXACT_MOVINGEB" : "FULL_REINIT"))
                      << " relChange=" << std::fixed << std::setprecision(4) << adaptation.cellCountChange
                      << " th=" << options_.amrGeomRefreshChangeRelTh
                      << std::defaultfloat << std::setprecision(6)
                      << std::endl;
        }
    }

}
