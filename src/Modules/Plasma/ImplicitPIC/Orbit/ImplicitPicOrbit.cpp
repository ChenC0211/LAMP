#include "Modules/Plasma/ImplicitPIC/Orbit/ImplicitPicOrbit.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/ExplicitPIC/Field/MatchedGather.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicOverrides.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
using Clock = std::chrono::steady_clock;
double secondsBetween(const Clock::time_point& a, const Clock::time_point& b)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(b - a).count();
}
void validateCacheOrThrow(const CouplingWeightsCache& cache, std::size_t expectedNp,
    const char* stageTag, std::size_t speciesIdx)
{
    const std::size_t np = cache.numParticles();
    const std::size_t nnz = cache.nnz();
    const auto& off = cache.offsets();
    const auto& idx = cache.cellIndices();
    const auto& w = cache.weights();

    auto fail = [&](const std::string& reason) {
        std::ostringstream oss;
        oss << "[ImplicitPIC] ERROR: invalid CouplingWeightsCache at " << stageTag
            << " (species=" << speciesIdx << ") : " << reason
            << " | expectedNp=" << expectedNp
            << " np=" << np
            << " off.size=" << off.size()
            << " nnz=" << nnz
            << " idx.size=" << idx.size()
            << " w.size=" << w.size();
        if (!off.empty())
            oss << " off[0]=" << off.front() << " off.back=" << off.back();
        throw std::runtime_error(oss.str());
    };

    if (np != expectedNp)
        fail("numParticles mismatch");
    if (off.size() != np + 1u)
        fail("offsets size must be np+1");
    if (off.empty() || off.front() != 0)
        fail("offsets[0] must be 0");
    if (idx.size() != w.size())
        fail("cellIndices/weights size mismatch");
    if (static_cast<std::size_t>(off.back()) != nnz)
        fail("offsets.back must equal nnz");
    for (std::size_t i = 1; i < off.size(); ++i)
    {
        if (off[i] < off[i - 1])
            fail("offsets must be non-decreasing");
        if (off[i] < 0)
            fail("offsets must be non-negative");
        if (static_cast<std::size_t>(off[i]) > nnz)
            fail("offsets must not exceed nnz");
    }
}
}

ImplicitPicOrbit::ImplicitPicOrbit(const QuadtreeMesh& mesh, std::size_t speciesCount)
    : mesh(mesh), cachesBar(speciesCount), orbitAvgWorkspaces(speciesCount),
      gatherScratch(speciesCount), predictorX0Scratch(speciesCount), predictorX1Scratch(speciesCount),
      predictorHostScratch(speciesCount), segOffsetsScratch(speciesCount), segsScratch(speciesCount)
{}

ImplicitPicOrbit::ImplicitPicOrbit(const QuadtreeMesh& mesh, std::size_t speciesCount,
    const PicConfig& config, const ImplicitPicOverrides& overrides,
    const UnitsContext* units, bool periodicX, bool periodicY)
    : ImplicitPicOrbit(mesh, speciesCount)
{
    auto& averaging = configuredOptions_.averaging;
    averaging.renormalize = true;
    // Patching requires the fused path, regardless of the fused override.
    averaging.useFusedPath = overrides.orbitAvgFused.value_or(config.orbitAvgFused) || config.orbitAvgPatch.enabled;
    averaging.adaptiveEnabled = config.orbitAvgAdaptive.enabled;
    averaging.adaptiveTol = config.orbitAvgAdaptive.tol;
    averaging.adaptiveMaxRefineDepth = config.orbitAvgAdaptive.maxRefineDepth;
    averaging.adaptiveRule = config.orbitAvgAdaptive.rule;
    averaging.adaptiveCheapRefine = config.orbitAvgAdaptive.cheapRefine;
    averaging.adaptiveCrossSplitAtMaxDepth = config.orbitAvgAdaptive.crossSplitAtMaxDepth;
    averaging.adaptiveCrossSplitMaxIters = config.orbitAvgAdaptive.crossSplitMaxIters;
    double sheathPadCode = config.orbitAvgAdaptive.sheathPad_m;
    if (units) sheathPadCode = units->LengthSIToCode(config.orbitAvgAdaptive.sheathPad_m);
    averaging.adaptiveCrossSplitSheathPad = std::max(0.0, sheathPadCode);
    averaging.adaptiveCrossSplitMaxSplits = config.orbitAvgAdaptive.crossSplitMaxSplits;
    averaging.segmentPatchEnabled = config.orbitAvgPatch.enabled;
    averaging.segmentPatchTol = config.orbitAvgPatch.tol;
    averaging.segmentPatchMaxRefineDepth = config.orbitAvgPatch.maxRefineDepth;

    configuredOptions_.parallelSpecies = overrides.orbitAvgSpeciesPar.value_or(false);
    auto& segmenter = configuredOptions_.segmenter;
    segmenter.trace = overrides.segmenterTrace;
    segmenter.ignoreHostHints = overrides.segmenterIgnoreHostHints;
    segmenter.periodicX = periodicX;
    segmenter.periodicY = periodicY;
    segmenter.maxSegmentsPerParticle = std::max(16, overrides.segmenterMaxSegments.value_or(std::max(16, config.segmenterMaxSegmentsPerParticle)));
    segmenter.segmentCellMaxLevel = overrides.segmenterMaxLevel.value_or(config.segmenterMaxLevel);
    auto& weights = configuredOptions_.weights;
    weights.periodicX = periodicX;
    weights.periodicY = periodicY;
    weights.enableCompositeBand = config.deposition.compositeBand;
    weights.enableMomentExactCut = config.deposition.momentExactCut;
}

ImplicitPicOrbit::Options ImplicitPicOrbit::optionsFor(const PicEmbeddedGeometry& geometry,
    bool useCutMetrics, bool enableEmbeddedIntercept) const
{
    auto options = configuredOptions_;
    if (useCutMetrics)
    {
        options.segmenter.enableEmbeddedIntercept = enableEmbeddedIntercept;
        options.segmenter.embeddedLevelSetVert = geometry.levelSet();
        options.segmenter.embeddedInsideIsNegative = geometry.insideIsNegative();
        options.segmenter.alphaByCellIndex = &geometry.alphaByCell();
        options.weights.levelSetI = geometry.levelSetIndex();
        options.weights.levelSetVert = geometry.levelSet();
        options.weights.insideIsNegative = geometry.insideIsNegative();
        options.weights.alphaByCellIndex = &geometry.alphaByCell();
        options.weights.cutByAlpha = true;
    }
    return options;
}

std::vector<ImplicitPicOrbit::SpeciesPerf> ImplicitPicOrbit::build(const std::vector<ParticleType>& species,
    const std::vector<double>& hp, double dt, const PicEmbeddedGeometry& geometry,
    bool useCutMetrics, bool enableEmbeddedIntercept, Trace trace)
{
    return build(species, hp, dt, optionsFor(geometry, useCutMetrics, enableEmbeddedIntercept), trace);
}

void ImplicitPicOrbit::resetScratch()
{
    const auto count = cachesBar.size();
    cachesBar.clear();
    cachesBar.resize(count);
    orbitAvgWorkspaces.assign(count, {});
    for (auto& buf : gatherScratch) buf.clear();
    for (auto& buf : predictorX0Scratch) buf.clear();
    for (auto& buf : predictorX1Scratch) buf.clear();
    for (auto& buf : predictorHostScratch) buf.clear();
    for (auto& buf : segOffsetsScratch) buf.clear();
    for (auto& buf : segsScratch) buf.clear();
}

void ImplicitPicOrbit::validateCurrentCaches(const std::vector<ParticleType>& species,
    const std::vector<CouplingWeightsCache>& caches)
{
    if (caches.size() < species.size())
        throw std::invalid_argument("Implicit PIC current cache species count mismatch");
    for (std::size_t si = 0; si < species.size(); ++si)
        validateCacheOrThrow(caches[si], species[si].particle.size(), "buildCachesAtCurrentPositions", si);
}

void ImplicitPicOrbit::gather(const std::vector<ParticleType>& species,
    const std::vector<CouplingWeightsCache>& cachesN, const CellField<Vector2D>& Ecell, int step)
{
    if (species.size() != cachesBar.size() || cachesN.size() < species.size())
        throw std::invalid_argument("Implicit PIC gather species count mismatch");
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        gatherScratch[si].clear();
        ExplicitPIC::Field::MatchedGather::gatherVector2D(cachesN[si], Ecell, gatherScratch[si], false);
        if (step == 0 && gatherScratch[si].size() != species[si].particle.size())
        {
            std::ostringstream oss;
            oss << "[ImplicitPIC] ERROR: gather size mismatch at gather E^n (species=" << si << ")"
                << " cache.np=" << cachesN[si].numParticles()
                << " particles=" << species[si].particle.size()
                << " Ep=" << gatherScratch[si].size();
            throw std::runtime_error(oss.str());
        }
    }
}

std::vector<ImplicitPicOrbit::SpeciesPerf> ImplicitPicOrbit::build(const std::vector<ParticleType>& species,
    const std::vector<double>& hp, double dt, const Options& options, Trace trace)
{
    if (species.size() != cachesBar.size())
        throw std::invalid_argument("Implicit PIC orbit species count mismatch");
    if (options.segmenter.periodicX != options.weights.periodicX
        || options.segmenter.periodicY != options.weights.periodicY)
        throw std::invalid_argument("Implicit PIC predictor and weights must use the same periodicity");
    const int step = trace.step;
    const bool consoleInfo = trace.info, consoleDebug = trace.debug;
    const bool periodicX = options.weights.periodicX, periodicY = options.weights.periodicY;
    const auto& orbitAvgBuildOpt = options.averaging;
    const bool adaptiveOrbitAvgAllowed = orbitAvgBuildOpt.adaptiveEnabled;
    const bool adaptiveOrbitAvgNeedsEbClip = adaptiveOrbitAvgAllowed && options.segmenter.enableEmbeddedIntercept;
    const bool orbitAvgSpeciesParallel = options.parallelSpecies;
    std::vector<SpeciesPerf> orbitSpeciesPerf(species.size());
    std::exception_ptr orbitBuildError;

#ifdef _OPENMP
#pragma omp parallel for if(orbitAvgSpeciesParallel && species.size() > 2u) schedule(static)
#endif
    for (int sii = 0; sii < static_cast<int>(species.size()); ++sii)
    {
        const std::size_t si = static_cast<std::size_t>(sii);
        auto& perf = orbitSpeciesPerf[si];
        try
        {
        auto& sp = species[si];
        const int n = static_cast<int>(sp.particle.size());
        if (n <= 0)
        {
            cachesBar[si].clear();
            continue;
        }

        auto& x0Pred = predictorX0Scratch[si];
        auto& x1Pred = predictorX1Scratch[si];
        auto& host = predictorHostScratch[si];
        x0Pred.resize(static_cast<std::size_t>(n));
        x1Pred.resize(static_cast<std::size_t>(n));
        host.resize(static_cast<std::size_t>(n));

        // Predictor prep: build linear endpoints consistent with leapfrog kick+drift.
        const bool useAdaptiveOrbitAvg = adaptiveOrbitAvgAllowed;

        const auto tPredictPrep0 = Clock::now();
        const double kappa = sp.kappa_s();

        auto wrap01 = [](double a) -> double {
            a = std::fmod(a, 1.0);
            if (a < 0.0)
                a += 1.0;
            if (a >= 1.0)
                a = std::nextafter(1.0, 0.0);
            return a;
        };
        auto wrapCoord = [&](double& v, double vmin, double vmax, double eps) {
            const double L = vmax - vmin;
            if (L > 0.0)
            {
                const double s = (v - vmin) / L;
                v = vmin + wrap01(s) * L;
                if (v >= vmax)
                    v = vmax - std::abs(eps);
            }
        };

        auto minimalImage = [](double dx, double L) -> double {
            if (!(L > 0.0))
                return dx;
            // Map dx to [-0.5L, 0.5L] (ties toward even via nearbyint).
            const double k = std::nearbyint(dx / L);
            return dx - k * L;
        };

        const double epsWrap = 1e-12;
        const double xmin = mesh.minCoord.x;
        const double xmax = mesh.maxCoord.x;
        const double ymin = mesh.minCoord.y;
        const double ymax = mesh.maxCoord.y;
        const double Lx = xmax - xmin;
        const double Ly = ymax - ymin;

        const auto& Ep = gatherScratch[si];
        if (static_cast<int>(Ep.size()) != n)
            throw std::runtime_error("[ImplicitPIC] predictor Ep size mismatch");

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (int p = 0; p < n; ++p)
        {
            const auto& part = sp.particle[static_cast<std::size_t>(p)];
            const Vector2D x0p = part.position;
            Vector2D vHalf = part.velocity;
            vHalf.x += kappa * Ep[static_cast<std::size_t>(p)].x * dt;
            vHalf.y += kappa * Ep[static_cast<std::size_t>(p)].y * dt;
            Vector2D x1p = Vector2D(x0p.x + vHalf.x * dt, x0p.y + vHalf.y * dt);
            if (useAdaptiveOrbitAvg)
            {
                if (periodicX)
                {
                    const double dx = minimalImage(x1p.x - x0p.x, Lx);
                    x1p.x = x0p.x + dx;
                    wrapCoord(x1p.x, xmin, xmax, epsWrap);
                }
                if (periodicY)
                {
                    const double dy = minimalImage(x1p.y - x0p.y, Ly);
                    x1p.y = x0p.y + dy;
                    wrapCoord(x1p.y, ymin, ymax, epsWrap);
                }
            }
            x0Pred[static_cast<std::size_t>(p)] = x0p;
            x1Pred[static_cast<std::size_t>(p)] = x1p;
            host[static_cast<std::size_t>(p)] = part.inCellL;
        }
        perf.predictorPrep_s += secondsBetween(tPredictPrep0, Clock::now());
        const int predictorSegmentBudget = std::max(16, options.segmenter.maxSegmentsPerParticle);
        auto runPredictorSegmenter = [&](const char* tag,
                                         ImplicitPIC::Orbit::PredictorOrbitSegmenter::Options& sopt)
            -> ImplicitPIC::Orbit::PredictorOrbitSegmenter::Stats {
            sopt.maxSegmentsPerParticle = predictorSegmentBudget;
            const auto tPredictSeg0 = Clock::now();
            ImplicitPIC::Orbit::PredictorOrbitSegmenter::Stats segStats;
            try
            {
                ImplicitPIC::Orbit::PredictorOrbitSegmenter::segmentLinearOrbits(
                    mesh,
                    /*x0=*/x0Pred,
                    /*x1=*/x1Pred,
                    /*hostCellHints=*/&host,
                    dt,
                    segOffsetsScratch[si],
                    segsScratch[si],
                    sopt,
                    &segStats);
            }
            catch (const std::runtime_error& e)
            {
                const std::string msg = e.what();
                const bool hitSegCap =
                    (msg.find("exceeded maxSegmentsPerParticle") != std::string::npos);
                const int retryBudget =
                    std::min(131072,
                             std::max(sopt.maxSegmentsPerParticle * 8, 16384));
                if (!hitSegCap || retryBudget <= sopt.maxSegmentsPerParticle)
                    throw;

                if (consoleInfo)
                {
                    std::cout << "[ImplicitPIC][segmenter] retry " << tag
                              << " at step=" << (step + 1)
                              << " species=" << si
                              << " budget " << sopt.maxSegmentsPerParticle
                              << " -> " << retryBudget
                              << std::endl;
                }

                sopt.maxSegmentsPerParticle = retryBudget;
                segOffsetsScratch[si].clear();
                segsScratch[si].clear();
                segStats = ImplicitPIC::Orbit::PredictorOrbitSegmenter::Stats{};
                ImplicitPIC::Orbit::PredictorOrbitSegmenter::segmentLinearOrbits(
                    mesh,
                    /*x0=*/x0Pred,
                    /*x1=*/x1Pred,
                    /*hostCellHints=*/&host,
                    dt,
                    segOffsetsScratch[si],
                    segsScratch[si],
                    sopt,
                    &segStats);
            }
            perf.predictorSegment_s += secondsBetween(tPredictSeg0, Clock::now());
            perf.predictor_segments_total += static_cast<std::uint64_t>(segsScratch[si].size());
            perf.predictor_fast_exit_hits += segStats.fastExitHits;
            perf.predictor_bisection_fallback_hits += segStats.bisectionFallbackHits;
            perf.predictor_coarsen_hits += segStats.coarsenHits;
            perf.predictor_coarsen_parent_climb += segStats.coarsenParentClimb;
            return segStats;
        };

        if (!useAdaptiveOrbitAvg)
        {
            auto sopt = options.segmenter;
            if (sopt.enableEmbeddedIntercept)
                sopt.segmentCellMaxLevel = -1;

            if (consoleDebug && step == 0)
            {
                std::cout << "[ImplicitPIC] segmenter: begin (species=" << si
                          << ", n=" << n
                          << ", ebIntercept=" << (sopt.enableEmbeddedIntercept ? "on" : "off")
                          << ", maxSeg=" << predictorSegmentBudget
                          << ")" << std::endl;
            }
            runPredictorSegmenter("baseline", sopt);

            if (consoleDebug && step == 0)
            {
                const int nSeg = static_cast<int>(segsScratch[si].size());
                std::cout << "[ImplicitPIC] segmenter: end (species=" << si
                          << ", segments=" << nSeg << ")" << std::endl;
            }
        }
        else if (adaptiveOrbitAvgNeedsEbClip)
        {
            // EB-aware clipping for adaptive whole-step orbit averaging:
            // run interception segmenter once, then collapse each particle
            // predictor endpoint to the last in-plasma segment endpoint.
            auto sopt = options.segmenter;
            sopt.segmentCellMaxLevel = -1;

            if (consoleDebug && step == 0)
            {
                std::cout << "[ImplicitPIC] segmenter(adaptive-eb-clip): begin (species=" << si
                          << ", n=" << n
                          << ", maxSeg=" << predictorSegmentBudget
                          << ")" << std::endl;
            }
            runPredictorSegmenter("adaptive-eb-clip", sopt);

            const auto& segOffsets = segOffsetsScratch[si];
            const auto& segs = segsScratch[si];
            for (int p = 0; p < n; ++p)
            {
                const int s0 = segOffsets[static_cast<std::size_t>(p)];
                const int s1 = segOffsets[static_cast<std::size_t>(p) + 1u];
                if (s1 > s0)
                {
                    const auto& firstSeg = segs[static_cast<std::size_t>(s0)];
                    const auto& lastSeg = segs[static_cast<std::size_t>(s1 - 1)];
                    host[static_cast<std::size_t>(p)] = firstSeg.cellLabel;
                    x1Pred[static_cast<std::size_t>(p)] = lastSeg.x1;
                }
                else
                {
                    // No in-plasma segment survived (e.g. starts in solid):
                    // collapse predictor path to avoid integrating through solid.
                    x1Pred[static_cast<std::size_t>(p)] = x0Pred[static_cast<std::size_t>(p)];
                }
            }

            if (consoleDebug && step == 0)
            {
                const int nSeg = static_cast<int>(segsScratch[si].size());
                std::cout << "[ImplicitPIC] segmenter(adaptive-eb-clip): end (species=" << si
                          << ", segments=" << nSeg << ")" << std::endl;
            }
        }
        else
        {
            // Report as one interval per particle (consistent perf metric).
            perf.predictor_segments_total += static_cast<std::uint64_t>(n);
        }

        // (4) Orbit-averaged overlaps -> cacheBar.
        if (consoleDebug && step == 0)
            std::cout << "[ImplicitPIC] stage: orbitAvgOverlaps (species=" << si << ")" << std::endl;
        auto wopt = options.weights;
        wopt.recordHpStats = false;
        wopt.recordPerParticleHp = false;
        wopt.recordPerParticleCompositeBand = false;

        const std::size_t np = static_cast<std::size_t>(n);
        const auto tOrbitAvg0 = Clock::now();
        ImplicitPIC::Orbit::OrbitAveragedWeightsBuilderTiming orbitAvgTiming;
        bool okW = false;
        if (useAdaptiveOrbitAvg)
        {
            okW = ImplicitPIC::Orbit::buildOrbitAveragedCacheAdaptiveGauss24(
                mesh, hp, wopt,
                np,
                x0Pred,
                x1Pred,
                host,
                dt,
                cachesBar[si],
                &orbitAvgWorkspaces[si],
                orbitAvgBuildOpt,
                &orbitAvgTiming);
        }
        else
        {
            okW = ImplicitPIC::Orbit::buildOrbitAveragedCacheGauss2(
                mesh, hp, wopt,
                np,
                segOffsetsScratch[si],
                segsScratch[si],
                host,
                dt,
                cachesBar[si],
                &orbitAvgWorkspaces[si],
                orbitAvgBuildOpt,
                &orbitAvgTiming);
        }
        perf.orbitAvgWeights_s += secondsBetween(tOrbitAvg0, Clock::now());
        perf.orbitAvg_fused_mode = std::max(perf.orbitAvg_fused_mode, orbitAvgTiming.fusedMode);
        perf.orbitAvg_fused_eval_s += orbitAvgTiming.fusedEval_s;
        perf.orbitAvg_fused_accum_s += orbitAvgTiming.fusedAccum_s;
        perf.orbitAvg_sampleFill_s += orbitAvgTiming.sampleFill_s;
        perf.orbitAvg_sampleCacheBuild_s += orbitAvgTiming.sampleCacheBuild_s;
        perf.orbitAvg_localReduce_s += orbitAvgTiming.localReduce_s;
        perf.orbitAvg_concat_s += orbitAvgTiming.concat_s;
        perf.orbitAvg_sampleCwFallbackCount += static_cast<std::uint64_t>(orbitAvgTiming.sampleCacheFallbackCount);
        perf.orbitAvg_sampleCount += static_cast<std::uint64_t>(orbitAvgTiming.sampleCount);
        perf.orbitAvg_sampleCache_nnz += static_cast<std::uint64_t>(orbitAvgTiming.sampleCacheNnz);
        perf.orbitAvg_out_nnz += static_cast<std::uint64_t>(orbitAvgTiming.outNnz);
        perf.orbitAvg_adaptAcceptCount += static_cast<std::uint64_t>(orbitAvgTiming.adaptAcceptCount);
        perf.orbitAvg_adaptRefineCount += static_cast<std::uint64_t>(orbitAvgTiming.adaptRefineCount);
        perf.orbitAvg_adaptHitMaxDepthCount += static_cast<std::uint64_t>(orbitAvgTiming.adaptHitMaxDepthCount);
        perf.orbitAvg_adaptAcceptedDepthSum += static_cast<std::uint64_t>(orbitAvgTiming.adaptAcceptedDepthSum);
        if (!okW)
            throw std::runtime_error("[ImplicitPIC] buildOrbitAveragedCache failed");

        if (step == 0)
            validateCacheOrThrow(cachesBar[si], np,
                                 useAdaptiveOrbitAvg ? "buildOrbitAveragedCacheAdaptiveGauss24"
                                                     : "buildOrbitAveragedCacheGauss2",
                                 si);
        }
        catch (...)
        {
#ifdef _OPENMP
            #pragma omp critical(quantum_orbit_species_build_error)
            {
                if (!orbitBuildError)
                    orbitBuildError = std::current_exception();
            }
#else
            if (!orbitBuildError)
                orbitBuildError = std::current_exception();
#endif
    }

}


    if (orbitBuildError)
        std::rethrow_exception(orbitBuildError);

    return orbitSpeciesPerf;
}
