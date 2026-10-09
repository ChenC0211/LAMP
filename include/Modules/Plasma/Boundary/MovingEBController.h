#include "Core/IO/FileManagement.h"
#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Core/Base/Interpolation.h"
#include "Core/Base/Vector2D.h"
#include "Core/Field/VertField.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Modules/Plasma/Boundary/MovingEBConfig.h"
#include "Modules/Plasma/Particles/Particle.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Units/UnitsContext.h"
#include "Solvers/Poisson/CutCellMetrics.h"

struct WallImpactSample
{
    int boundaryId = 0;
    int levelSetIndex = -1;
    int speciesId = -1;
    int lineIndex = -1;
    int cellLabel = -1;
    Vector2D pos;
    Vector2D normal;
    Vector2D vImpact;
    double vn = 0.0;
    double energy_eV = 0.0;
    double theta_rad = 0.0;
    double macroWeight = 0.0;
    double time_s = 0.0;
    bool crossing = false;
    bool ionSpecies = false;
    bool usedCutMetrics = false;
    bool resolved = false;
};

class MovingEBController
{
public:
    struct ImpactResolveContext
    {
        struct LineBin
        {
            double minX = 0.0;
            double minY = 0.0;
            double cellSize = 1.0;
            int nx = 0;
            int ny = 0;
            std::vector<std::vector<int>> bins;
            std::unordered_map<int, std::vector<int>> boundaryLineIds;
        };

        bool enabled = false;
        LineBin lineBin;
    };

    struct GeometryDiagnostics
    {
        double gsa_mm = 0.0;
        double gad_mm = 0.0;
        double da_mm = 0.0;
    };

    struct RegionDiagnostics
    {
        std::string name;
        int boundaryId = 0;
        int levelSetIndex = -1;
        bool active = false;
        bool normalErosion = false;
        double sampleCount = 0.0;
        double area_m2 = 0.0;
        double flux_1pm2s = 0.0;
        double meanEnergy_eV = 0.0;
        double meanCos = 1.0;
        double vn_mps = 0.0;
    };

    explicit MovingEBController(const DictionaryGeometry& dictGeometry, MovingEBConfig config,
                                const UnitsContext* units = nullptr,
                                const NormalizationContext* normalization = nullptr)
        : dictGeometry_(&dictGeometry), config_(std::move(config)), units_(units), normalization_(normalization), enabled_(config_.enabled)
    {
        states_.resize(config_.regions.size());
        regionCaches_.resize(config_.regions.size());
        for (std::size_t i = 0; i < config_.regions.size(); ++i)
        {
            states_[i].name = config_.regions[i].name;
            const auto& region = config_.regions[i];
            auto& cache = regionCaches_[i];
            cache.motionMode = toLower_(region.motionMode);
            cache.model = toLower_(region.model);
            cache.energyModel = toLower_(region.energyModel);
            cache.angleModel = toLower_(region.angleModel);
            cache.lineIds = region.selector.lineIds;
            if (cache.lineIds.empty() && region.boundaryId > 0)
            {
                auto it = dictGeometry_->boundaryById.find(region.boundaryId);
                if (it != dictGeometry_->boundaryById.end())
                    cache.lineIds = it->second.lineIds;
            }
        }
    }

    bool enabled() const { return enabled_; }
    int resolvedImpactCount() const { return resolvedImpactsStep_; }
    int unresolvedImpactCount() const { return unresolvedImpactsStep_; }
    int ignoredNonIonImpactCount() const { return ignoredNonIonImpactsStep_; }
    int regionMissImpactCount() const { return regionMissImpactsStep_; }
    bool linesInCurrentGeometryFrame() const { return linesCurrentGeometryFrame_; }
    bool lastUpdateUsedAdvectedPath() const { return lastUpdateUsedAdvectedPath_; }
    bool lastUpdateHadActiveNormalErosion() const { return lastUpdateHadActiveNormalErosion_; }
    GeometryDiagnostics geometryDiagnostics() const
    {
        return GeometryDiagnostics{lastGeometryMetrics_.gsa_mm, lastGeometryMetrics_.gad_mm, lastGeometryMetrics_.da_mm};
    }

    std::vector<RegionDiagnostics> activeRegionDiagnostics(int step) const
    {
        std::vector<RegionDiagnostics> out;
        if (dictGeometry_ == nullptr)
            return out;
        for (std::size_t ri = 0; ri < config_.regions.size() && ri < states_.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            const auto& state = states_[ri];
            if (!regionIsActive_(region, step))
                continue;
            RegionDiagnostics diag;
            diag.name = region.name;
            diag.boundaryId = region.boundaryId;
            diag.levelSetIndex = region.levelSetIndex;
            diag.active = state.activeThisStep;
            diag.normalErosion = (toLower_(region.motionMode) == "normalerosion");
            diag.sampleCount = state.lastSampleCount;
            diag.area_m2 = state.lastArea_m2;
            diag.flux_1pm2s = state.smoothedFlux;
            diag.meanEnergy_eV = state.smoothedEnergy;
            diag.meanCos = state.smoothedCos;
            diag.vn_mps = state.lastNormalSpeed_mps;
            out.push_back(diag);
        }
        return out;
    }

    std::vector<int> activeMovingLevelSetIndices(int step) const
    {
        std::vector<int> out;
        if (dictGeometry_ == nullptr)
            return out;
        out.reserve(config_.regions.size());
        for (const auto& region : config_.regions)
        {
            if (!regionIsActive_(region, step))
                continue;
            if (toLower_(region.motionMode) == "fixed")
                continue;
            if (region.levelSetIndex < 0)
                continue;
            out.push_back(region.levelSetIndex);
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    Vector2D boundaryRigidShiftCode(int boundaryId) const
    {
        Vector2D shift(0.0, 0.0);
        if (!enabled_ || boundaryId <= 0 || dictGeometry_ == nullptr)
            return shift;

        for (std::size_t ri = 0; ri < config_.regions.size() && ri < states_.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            if (region.boundaryId != boundaryId)
                continue;
            if (toLower_(region.motionMode) != "rigidshift")
                continue;
            shift.x += states_[ri].rigidShiftX_code;
            shift.y += states_[ri].rigidShiftY_code;
        }
        return shift;
    }

    template<typename Mesh>
    void initialize(const Mesh& mesh, const ImmersedGeometry<Mesh>& iGeometry, bool restarting = false)
    {
        if (!enabled_)
            return;
        if (!restarting)
        {
            initializeFields_(mesh);
            linesCurrentGeometryFrame_ = true;
            lastOutputStep_ = 0;
        }
        initializeNormalDirections_(mesh, iGeometry);
        lastGeometryMetrics_ = computeGeometryMetrics_(mesh, iGeometry);
        openOutputs_(restarting);
        if (!restarting)
        {
            writeGeometryMetrics_(0, 0.0, false);
            writeConsistencyMetrics_(0, 0.0, false, iGeometry);
            writeRegionMetrics_(0, 0.0, false);
            writeImpactStats_(0, 0.0, false);
        }
    }

    template<typename Mesh>
    void onMeshChanged(const Mesh& mesh)
    {
        if (!enabled_)
            return;
        initializeFields_(mesh);
    }

    void beginStep()
    {
        pendingSamples_.clear();
        resolvedImpactsStep_ = 0;
        unresolvedImpactsStep_ = 0;
        ignoredNonIonImpactsStep_ = 0;
        regionMissImpactsStep_ = 0;
        lastUpdateHadActiveNormalErosion_ = false;
        for (auto& st : states_)
        {
            st.activeThisStep = false;
            st.activeInWindow = false;
            st.activeWindowStepCount = 0;
        }
    }

    void ingest(const WallImpactSample& sample)
    {
        if (!enabled_)
            return;
        if (sample.resolved)
            ++resolvedImpactsStep_;
        else
        {
            ++unresolvedImpactsStep_;
            if (!sample.ionSpecies)
                ++ignoredNonIonImpactsStep_;
        }
        pendingSamples_.push_back(sample);
    }

    template<typename Mesh>
    bool maybeUpdate(int step, double time_s, double dt_s, Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry)
    {
        if (!enabled_)
            return false;

        const auto& cfg = config_;
        const int updateEvery = std::max(1, cfg.updateEvery);
        const bool cadenceHit = (step <= 0) || ((step % updateEvery) == 0);
        const bool outputCadenceHit = ((step % std::max(1, cfg.outputEvery)) == 0);
        bufferStepSamples_(dt_s);
        const int updateWindowStartExclusive = std::max(lastUpdateStep_, 0);
        const int outputWindowStartExclusive = std::max(lastOutputStep_, 0);
        const int updateWindowStepCount = std::max(1, step - updateWindowStartExclusive);
        const int outputWindowStepCount = std::max(1, step - outputWindowStartExclusive);
        const double updateWindowDuration_s = std::max(dt_s, updateWindowDuration_s_);
        const double outputWindowDuration_s = std::max(dt_s, outputWindowDuration_s_);
        const bool useExactRigidShift = cadenceHit && canUseExactRigidShiftPath_(updateWindowStartExclusive, step);
        lastUpdateUsedAdvectedPath_ = false;
        bool shouldUpdate = false;
        bool shouldOutput = false;

        if (cadenceHit || outputCadenceHit)
        {
            initializeNormalDirections_(mesh, iGeometry);
            if (cadenceHit)
            {
                accumulateRegionStats_(updateWindowSamples_, step, updateWindowStartExclusive, updateWindowDuration_s,
                                       iGeometry);
                const MotionSummary motion = updateRegionMotionState_(step, updateWindowStartExclusive);
                shouldUpdate = cadenceHit && motion.anyActive && motion.anyMotion;
                shouldOutput = shouldUpdate || outputCadenceHit;
                if (shouldUpdate)
                    buildVelocityFields_(mesh, iGeometry, updateWindowStepCount);
            }
            else
            {
                accumulateRegionStats_(outputWindowSamples_, step, outputWindowStartExclusive, outputWindowDuration_s,
                                       iGeometry);
                updateRegionMotionState_(step, outputWindowStartExclusive);
                shouldOutput = true;
            }
        }

        if (shouldUpdate)
        {
            for (std::size_t ri = 0; ri < config_.regions.size() && ri < states_.size(); ++ri)
            {
                const auto& region = config_.regions[ri];
                auto& state = states_[ri];
                if (regionCaches_[ri].motionMode != "rigidshift")
                    continue;
                if (!state.activeInWindow || state.activeWindowStepCount <= 0)
                    continue;
                const double activeDuration_s =
                    updateWindowDuration_s *
                    (static_cast<double>(state.activeWindowStepCount) / static_cast<double>(updateWindowStepCount));
                state.rigidShiftX_code += state.lastRigidVx_codeps * activeDuration_s;
                state.rigidShiftY_code += state.lastRigidVy_codeps * activeDuration_s;
            }
            if (useExactRigidShift)
            {
                applyExactRigidShiftGeometry_(iGeometry);
                // Exact rigid-shift updates mutate iGeometry.lines directly.
                // Rebuild level-set values from the moved interface before
                // refreshing cut topology; reinitializeLevelSets() alone only
                // re-distances the existing field and would otherwise keep
                // stale cut topology from the pre-shift geometry.
                iGeometry.initializeLevelSets(mesh);
                // The signed field has already been rebuilt from the exact
                // moved geometry above. Refresh only derivatives/cut-info here;
                // PDE reinitialization is reserved for advected level sets.
                iGeometry.refreshLevelSetDerivedData(mesh, time_s, /*reinitializeMovingLevelSets=*/false);
                linesCurrentGeometryFrame_ = true;
            }
            else
            {
                iGeometry.updateLevelSets(mesh, currentVx_, currentVy_, previousVx_, previousVy_, updateWindowDuration_s);
                iGeometry.markPhases(mesh);
                // Exact redistancing from the advected zero contour replaces
                // PDE reinit to avoid far-field roughness accumulation.
                iGeometry.exactRedistanceMovingLevelSets(mesh);
                iGeometry.refreshLevelSetDerivedData(mesh, time_s, /*reinitializeMovingLevelSets=*/false);
                linesCurrentGeometryFrame_ = false;
                lastUpdateUsedAdvectedPath_ = true;
            }
            iGeometry.checkCutInfo(mesh, true);
            if (useExactRigidShift)
                iGeometry.markPhases(mesh);
            lastUpdateStep_ = step;
            lastUpdateTime_s_ = time_s;
            lastGeometryMetrics_ = computeGeometryMetrics_(mesh, iGeometry);
        }

        if (shouldOutput)
        {
            writeGeometryMetrics_(step, time_s, shouldUpdate);
            writeConsistencyMetrics_(step, time_s, shouldUpdate, iGeometry);
            writeRegionMetrics_(step, time_s, shouldUpdate);
            writeImpactStats_(step, time_s, shouldUpdate);
            writeMatchStats_(step, time_s, shouldUpdate);
            lastOutputStep_ = step;
        }

        if (cadenceHit)
            clearUpdateWindow_();
        if (shouldOutput)
            clearOutputWindow_();
        pendingSamples_.clear();
        return shouldUpdate;
    }

    template<typename Mesh>
    static ImpactResolveContext buildImpactResolveContext(const ImmersedGeometry<Mesh>& iGeometry,
                                                         const std::vector<int>& lineToBoundaryId)
    {
        ImpactResolveContext ctx;
        ctx.enabled = !iGeometry.lines.empty();
        if (!ctx.enabled)
            return ctx;

        double minX = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();
        for (std::size_t li = 0; li < iGeometry.lines.size(); ++li)
        {
            const auto& ln = iGeometry.lines[li];
            minX = std::min(minX, std::min(ln.points[0].x, ln.points[1].x));
            maxX = std::max(maxX, std::max(ln.points[0].x, ln.points[1].x));
            minY = std::min(minY, std::min(ln.points[0].y, ln.points[1].y));
            maxY = std::max(maxY, std::max(ln.points[0].y, ln.points[1].y));
            const int bid = (li < lineToBoundaryId.size()) ? lineToBoundaryId[li] : 0;
            ctx.lineBin.boundaryLineIds[bid].push_back(static_cast<int>(li));
        }

        if (!std::isfinite(minX) || !std::isfinite(maxX) || !std::isfinite(minY) || !std::isfinite(maxY))
            return ctx;

        const double spanX = std::max(maxX - minX, 1.0e-12);
        const double spanY = std::max(maxY - minY, 1.0e-12);
        const int target = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(iGeometry.lines.size()))));
        const double cellSize = std::max(1.0e-12, std::max(spanX, spanY) / static_cast<double>(target));
        ctx.lineBin.minX = minX;
        ctx.lineBin.minY = minY;
        ctx.lineBin.cellSize = cellSize;
        ctx.lineBin.nx = std::max(1, static_cast<int>(std::ceil(spanX / cellSize)));
        ctx.lineBin.ny = std::max(1, static_cast<int>(std::ceil(spanY / cellSize)));
        ctx.lineBin.bins.assign(static_cast<std::size_t>(ctx.lineBin.nx * ctx.lineBin.ny), {});

        for (std::size_t li = 0; li < iGeometry.lines.size(); ++li)
        {
            const auto& ln = iGeometry.lines[li];
            const double lx0 = std::min(ln.points[0].x, ln.points[1].x);
            const double lx1 = std::max(ln.points[0].x, ln.points[1].x);
            const double ly0 = std::min(ln.points[0].y, ln.points[1].y);
            const double ly1 = std::max(ln.points[0].y, ln.points[1].y);
            const int ix0 = clampBinIndex_(static_cast<int>(std::floor((lx0 - ctx.lineBin.minX) / cellSize)), ctx.lineBin.nx);
            const int ix1 = clampBinIndex_(static_cast<int>(std::floor((lx1 - ctx.lineBin.minX) / cellSize)), ctx.lineBin.nx);
            const int iy0 = clampBinIndex_(static_cast<int>(std::floor((ly0 - ctx.lineBin.minY) / cellSize)), ctx.lineBin.ny);
            const int iy1 = clampBinIndex_(static_cast<int>(std::floor((ly1 - ctx.lineBin.minY) / cellSize)), ctx.lineBin.ny);
            for (int iy = iy0; iy <= iy1; ++iy)
            {
                for (int ix = ix0; ix <= ix1; ++ix)
                    ctx.lineBin.bins[static_cast<std::size_t>(iy * ctx.lineBin.nx + ix)].push_back(static_cast<int>(li));
            }
        }

        return ctx;
    }

    std::filesystem::path resolvedCachePath(const std::filesystem::path& baseDir) const
    {
        if (dictGeometry_ == nullptr)
            return baseDir / cacheFileName();

        const std::string& configured = config_.restartPath;
        if (configured.empty())
            return baseDir / cacheFileName();

        const std::filesystem::path configuredPath(configured);
        if (configuredPath.is_absolute())
            return configuredPath;
        return baseDir / configuredPath;
    }

    void writeCache(const std::filesystem::path& baseDir) const
    {
        if (!enabled_)
            return;
        const std::filesystem::path cachePath = resolvedCachePath(baseDir);
        if (!cachePath.parent_path().empty())
            std::filesystem::create_directories(cachePath.parent_path());
        auto out = FileIO::openOutputFile(cachePath, std::ios::binary | std::ios::out | std::ios::trunc);

        const char magic[4] = {'M', 'E', 'B', '5'};
        out.write(magic, 4);
        out.write(reinterpret_cast<const char*>(&lastUpdateStep_), sizeof(int));
        out.write(reinterpret_cast<const char*>(&lastOutputStep_), sizeof(int));
        out.write(reinterpret_cast<const char*>(&lastUpdateTime_s_), sizeof(double));
        const unsigned char linesCurrent = linesCurrentGeometryFrame_ ? 1u : 0u;
        out.write(reinterpret_cast<const char*>(&linesCurrent), sizeof(unsigned char));
        const int stateCount = static_cast<int>(states_.size());
        out.write(reinterpret_cast<const char*>(&stateCount), sizeof(int));
        for (const auto& st : states_)
        {
            writeString_(out, st.name);
            out.write(reinterpret_cast<const char*>(&st.lastArea_m2), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.lastSampleCount), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.smoothedFlux), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.smoothedEnergy), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.smoothedCos), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.lastNormalSpeed_mps), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.lastRigidVx_codeps), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.lastRigidVy_codeps), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.rigidShiftX_code), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.rigidShiftY_code), sizeof(double));
            out.write(reinterpret_cast<const char*>(&st.normalDirectionSign), sizeof(double));
            const unsigned char normalInit = st.normalDirectionInitialized ? 1u : 0u;
            out.write(reinterpret_cast<const char*>(&normalInit), sizeof(unsigned char));
        }
        writeField_(out, currentVx_);
        writeField_(out, currentVy_);
        writeField_(out, previousVx_);
        writeField_(out, previousVy_);
        FileIO::finishOutputFile(out, cachePath);
    }

    template<typename Mesh>
    bool readCache(const Mesh& mesh, const std::filesystem::path& baseDir)
    {
        if (!enabled_)
            return false;
        const std::filesystem::path cachePath = resolvedCachePath(baseDir);
        std::ifstream in(cachePath.string(), std::ios::binary | std::ios::in);
        if (!in)
            return false;

        char magic[4] = {};
        in.read(magic, 4);
        const std::string magicStr(magic, 4);
        if (magicStr != "MEB2" && magicStr != "MEB3" && magicStr != "MEB4" && magicStr != "MEB5")
            return false;

        initializeFields_(mesh);
        in.read(reinterpret_cast<char*>(&lastUpdateStep_), sizeof(int));
        if (magicStr == "MEB5")
            in.read(reinterpret_cast<char*>(&lastOutputStep_), sizeof(int));
        in.read(reinterpret_cast<char*>(&lastUpdateTime_s_), sizeof(double));
        if (magicStr != "MEB5")
            lastOutputStep_ = std::max(lastUpdateStep_, 0);
        if (magicStr == "MEB3" || magicStr == "MEB4" || magicStr == "MEB5")
        {
            unsigned char linesCurrent = 0u;
            in.read(reinterpret_cast<char*>(&linesCurrent), sizeof(unsigned char));
            linesCurrentGeometryFrame_ = (linesCurrent != 0u);
        }
        else
        {
            // Preserve pre-MEB3 behavior for legacy caches.
            linesCurrentGeometryFrame_ = false;
        }
        int stateCount = 0;
        in.read(reinterpret_cast<char*>(&stateCount), sizeof(int));
        if (stateCount != static_cast<int>(states_.size()))
            return false;
        for (auto& st : states_)
        {
            st.name = readString_(in);
            in.read(reinterpret_cast<char*>(&st.lastArea_m2), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.lastSampleCount), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.smoothedFlux), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.smoothedEnergy), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.smoothedCos), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.lastNormalSpeed_mps), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.lastRigidVx_codeps), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.lastRigidVy_codeps), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.rigidShiftX_code), sizeof(double));
            in.read(reinterpret_cast<char*>(&st.rigidShiftY_code), sizeof(double));
            if (magicStr == "MEB4" || magicStr == "MEB5")
            {
                in.read(reinterpret_cast<char*>(&st.normalDirectionSign), sizeof(double));
                unsigned char normalInit = 0u;
                in.read(reinterpret_cast<char*>(&normalInit), sizeof(unsigned char));
                st.normalDirectionInitialized = (normalInit != 0u);
            }
            else
            {
                st.normalDirectionSign = -1.0;
                st.normalDirectionInitialized = false;
            }
        }
        readField_(in, mesh, currentVx_);
        readField_(in, mesh, currentVy_);
        readField_(in, mesh, previousVx_);
        readField_(in, mesh, previousVy_);
        return static_cast<bool>(in);
    }

    static std::string cacheFileName() { return "moving_eb_cache.bin"; }

    template<typename Mesh>
    static bool resolveWallImpact(const Mesh& mesh,
                                  const ImmersedGeometry<Mesh>& iGeometry,
                                  const Solvers::Poisson::CutCellMetrics* cutMetrics,
                                  const std::vector<int>& lineToBoundaryId,
                                  const ImpactResolveContext* queryContext,
                                  const Vector2D& pos,
                                  int cellLabelHint,
                                  const Particle& particle,
                                  int speciesId,
                                  double speciesMass,
                                  double speciesCharge,
                                  double defaultWeight,
                                  double velocityScale_mps,
                                  double time_s,
                                  bool crossing,
                                  WallImpactSample& out)
    {
        auto pointSegDist2 = [](const Vector2D& pt, const Vector2D& a, const Vector2D& b) {
            const Vector2D ab = b - a;
            const double denom = ab.mod2();
            if (!(denom > 0.0))
                return (pt - a).mod2();
            const Vector2D ap = pt - a;
            const double proj = ap.x * ab.x + ap.y * ab.y;
            const double t = std::clamp(proj / denom, 0.0, 1.0);
            return (pt - (a + ab * t)).mod2();
        };
        auto lineDistance2 = [&](int li) {
            if (li < 0 || li >= static_cast<int>(iGeometry.lines.size()))
                return std::numeric_limits<double>::infinity();
            const auto& ln = iGeometry.lines[static_cast<std::size_t>(li)];
            return pointSegDist2(pos, ln.points[0], ln.points[1]);
        };
        auto nearestLineForBoundary = [&](int requiredBoundaryId, double* bestD2Out = nullptr) {
            double bestD2 = std::numeric_limits<double>::infinity();
            int nearest = -1;
            std::vector<int> candidateLines;
            if (queryContext != nullptr)
                collectImpactCandidateLines_(*queryContext, pos, requiredBoundaryId, candidateLines);

            const auto scanLine = [&](int li) {
                const int bid = (li >= 0 && li < static_cast<int>(lineToBoundaryId.size()))
                                    ? lineToBoundaryId[static_cast<std::size_t>(li)]
                                    : 0;
                if (requiredBoundaryId > 0 && bid != requiredBoundaryId)
                    return;
                const double d2 = lineDistance2(li);
                if (d2 < bestD2)
                {
                    bestD2 = d2;
                    nearest = li;
                }
            };

            for (int li : candidateLines)
                scanLine(li);

            if (nearest < 0)
            {
                if (requiredBoundaryId > 0 && queryContext != nullptr)
                {
                    auto it = queryContext->lineBin.boundaryLineIds.find(requiredBoundaryId);
                    if (it != queryContext->lineBin.boundaryLineIds.end())
                    {
                        for (int li : it->second)
                            scanLine(li);
                    }
                }
                else
                {
                    for (int li = 0; li < static_cast<int>(iGeometry.lines.size()); ++li)
                        scanLine(li);
                }
            }
            if (bestD2Out != nullptr)
                *bestD2Out = bestD2;
            return nearest;
        };

        out = WallImpactSample{};
        out.pos = pos;
        out.speciesId = speciesId;
        out.crossing = crossing;
        out.time_s = time_s;
        out.macroWeight = (particle.w > 0.0) ? particle.w : defaultWeight;
        out.ionSpecies = (speciesCharge > 0.0);

        int cellLabel = cellLabelHint;
        if (cellLabel < 0)
        {
            try { cellLabel = mesh.findCellLabel(pos); } catch (...) { cellLabel = -1; }
        }
        out.cellLabel = cellLabel;

        int boundaryId = 0;
        int bestLine = -1;
        Vector2D normal(0.0, 0.0);
        bool usedCutMetrics = false;
        if (cutMetrics != nullptr && cellLabel >= 0)
        {
            const int cellIndex = mesh.cells(cellLabel).index();
            const auto& cm = cutMetrics->cellMetricsByIndex();
            if (cellIndex >= 0 && cellIndex < static_cast<int>(cm.size()))
            {
                boundaryId = cm[static_cast<std::size_t>(cellIndex)].gammaBoundaryId;
                bestLine = cm[static_cast<std::size_t>(cellIndex)].gammaLineIndex;
                normal = cm[static_cast<std::size_t>(cellIndex)].n_gamma;
                usedCutMetrics = true;
            }
        }

        if (usedCutMetrics || !iGeometry.lines.empty())
        {
            double nearestGlobalD2 = std::numeric_limits<double>::infinity();
            const int nearestGlobalLine = nearestLineForBoundary(0, &nearestGlobalD2);
            const int nearestGlobalBoundary =
                (nearestGlobalLine >= 0 && nearestGlobalLine < static_cast<int>(lineToBoundaryId.size()))
                    ? lineToBoundaryId[static_cast<std::size_t>(nearestGlobalLine)]
                    : 0;

            if (usedCutMetrics)
            {
                const int taggedLine = bestLine;
                const int taggedBoundary0 = boundaryId;
                if (boundaryId <= 0 && bestLine >= 0 && bestLine < static_cast<int>(lineToBoundaryId.size()))
                    boundaryId = lineToBoundaryId[static_cast<std::size_t>(bestLine)];

                double nearestBoundaryD2 = std::numeric_limits<double>::infinity();
                const int nearestBoundaryLine =
                    (boundaryId > 0) ? nearestLineForBoundary(boundaryId, &nearestBoundaryD2) : -1;
                const int taggedBoundary =
                    (bestLine >= 0 && bestLine < static_cast<int>(lineToBoundaryId.size()))
                        ? lineToBoundaryId[static_cast<std::size_t>(bestLine)]
                        : 0;
                const double taggedD2 = lineDistance2(bestLine);
                const bool taggedLineValid = std::isfinite(taggedD2);
                const bool taggedBoundaryConsistent =
                    (boundaryId <= 0) || (taggedBoundary == 0) || (taggedBoundary == boundaryId);
                bool trustTaggedBoundary = (boundaryId > 0 && nearestBoundaryLine >= 0 && taggedBoundaryConsistent);

                if (trustTaggedBoundary && nearestGlobalLine >= 0 && nearestGlobalBoundary > 0 &&
                    nearestGlobalBoundary != boundaryId)
                {
                    const double boundarySlack = 1.25;
                    const double tol = 1e-24;
                    trustTaggedBoundary =
                        std::isfinite(nearestBoundaryD2) &&
                        nearestBoundaryD2 <= boundarySlack * nearestGlobalD2 + tol;
                }

                if (trustTaggedBoundary)
                {
                    bestLine = nearestBoundaryLine;
                }
                else
                {
                    bestLine = nearestGlobalLine;
                    if (nearestGlobalBoundary > 0)
                        boundaryId = nearestGlobalBoundary;
                }

                usedCutMetrics =
                    taggedLineValid &&
                    bestLine == taggedLine &&
                    (taggedBoundary0 <= 0 || boundaryId == taggedBoundary0);
            }
            else
            {
                bestLine = nearestGlobalLine;
                if (boundaryId <= 0 && nearestGlobalBoundary > 0)
                    boundaryId = nearestGlobalBoundary;
            }
        }

        if (bestLine >= 0 && normal.mod2() <= 0.0)
        {
            const int lsI = iGeometry.lines[static_cast<std::size_t>(bestLine)].levelSetIndex;
            if (lsI >= 0 && lsI < static_cast<int>(iGeometry.levelSets.size()))
            {
                try
                {
                    const Cell& c = mesh.cells(mesh.findCellLabel(pos));
                    const double gx = Interpolation::biLinearInterpolation(mesh, iGeometry.LS_dx(lsI), c, pos);
                    const double gy = Interpolation::biLinearInterpolation(mesh, iGeometry.LS_dy(lsI), c, pos);
                    const double gmag = std::sqrt(gx * gx + gy * gy);
                    if (gmag > 0.0 && std::isfinite(gmag))
                        normal = Vector2D(-gx / gmag, -gy / gmag);
                }
                catch (...) {}
            }
        }

        out.boundaryId = boundaryId;
        out.lineIndex = bestLine;
        out.levelSetIndex = (bestLine >= 0 && bestLine < static_cast<int>(iGeometry.lines.size()))
                                ? iGeometry.lines[static_cast<std::size_t>(bestLine)].levelSetIndex
                                : -1;
        out.normal = normal;
        out.usedCutMetrics = usedCutMetrics;

        const double vScale = (velocityScale_mps > 0.0) ? velocityScale_mps : 1.0;
        out.vImpact = Vector2D(particle.velocity.x * vScale, particle.velocity.y * vScale);
        const double nmag = std::sqrt(out.normal.mod2());
        if (nmag > 0.0)
        {
            out.vn = std::max(0.0, (out.vImpact.x * out.normal.x + out.vImpact.y * out.normal.y) / nmag);
            const double vz = particle.vz * vScale;
            const double speed = std::sqrt(out.vImpact.mod2() + vz * vz);
            const double cosTheta = std::clamp(out.vn / std::max(speed, 1e-30), -1.0, 1.0);
            out.theta_rad = std::acos(cosTheta);
        }
        if (speciesMass > 0.0)
        {
            const double vx = particle.velocity.x * vScale;
            const double vy = particle.velocity.y * vScale;
            const double vz = particle.vz * vScale;
            out.energy_eV = 0.5 * speciesMass * (vx * vx + vy * vy + vz * vz) / 1.602176634e-19;
        }
        out.resolved = (out.ionSpecies && boundaryId > 0 && out.levelSetIndex >= 0 && nmag > 0.0);
        return out.resolved;
    }

private:
    struct Bounds2D
    {
        double minX = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();

        bool valid() const
        {
            return std::isfinite(minX) && std::isfinite(maxX) && std::isfinite(minY) && std::isfinite(maxY) &&
                   minX <= maxX && minY <= maxY;
        }

        void include(const Vector2D& p)
        {
            minX = std::min(minX, p.x);
            maxX = std::max(maxX, p.x);
            minY = std::min(minY, p.y);
            maxY = std::max(maxY, p.y);
        }

        void expand(double margin)
        {
            if (!valid())
                return;
            minX -= margin;
            maxX += margin;
            minY -= margin;
            maxY += margin;
        }

        void intersectWith(const Bounds2D& other)
        {
            if (!valid())
            {
                *this = other;
                return;
            }
            if (!other.valid())
                return;
            minX = std::max(minX, other.minX);
            maxX = std::min(maxX, other.maxX);
            minY = std::max(minY, other.minY);
            maxY = std::min(maxY, other.maxY);
        }
    };

    struct RegionState
    {
        std::string name;
        double lastArea_m2 = 0.0;
        double lastSampleCount = 0.0;
        int lastResolvedIonCount = 0;
        int lastDirectMatchCount = 0;
        int lastFallbackMatchCount = 0;
        double smoothedFlux = 0.0;
        double smoothedEnergy = 0.0;
        double smoothedCos = 1.0;
        double lastNormalSpeed_mps = 0.0;
        double lastRigidVx_codeps = 0.0;
        double lastRigidVy_codeps = 0.0;
        double rigidShiftX_code = 0.0;
        double rigidShiftY_code = 0.0;
        double normalDirectionSign = -1.0;
        bool normalDirectionInitialized = false;
        bool activeThisStep = false;
        bool activeInWindow = false;
        int activeWindowStepCount = 0;
    };

    struct RegionCache
    {
        std::string motionMode;
        std::string model;
        std::string energyModel;
        std::string angleModel;
        std::vector<int> lineIds;
    };

    struct VertexSpatialIndex
    {
        bool valid = false;
        double minX = 0.0;
        double minY = 0.0;
        double cellSize = 1.0;
        int nx = 0;
        int ny = 0;
        std::vector<std::vector<int>> bins;
    };

    struct MotionSummary
    {
        bool anyActive = false;
        bool anyMotion = false;
    };

    struct GeometryMetrics
    {
        double gsa_mm = 0.0;
        double gad_mm = 0.0;
        double da_mm = 0.0;
    };

    struct MatchStatsRow
    {
        int boundaryId = 0;
        int levelSetIndex = -1;
        int resolvedIonCount = 0;
        int matchedResolvedIonCount = 0;
        int unmatchedResolvedIonCount = 0;
        int directMatchCount = 0;
        int fallbackMatchCount = 0;
    };

    const DictionaryGeometry* dictGeometry_ = nullptr;
    const MovingEBConfig config_;
    const UnitsContext* units_ = nullptr;
    const NormalizationContext* normalization_ = nullptr;
    bool enabled_ = false;
    std::vector<RegionState> states_;
    std::vector<RegionCache> regionCaches_;
    std::vector<WallImpactSample> pendingSamples_;
    std::vector<WallImpactSample> updateWindowSamples_;
    std::vector<WallImpactSample> outputWindowSamples_;
    VertField<double> currentVx_;
    VertField<double> currentVy_;
    VertField<double> previousVx_;
    VertField<double> previousVy_;
    VertexSpatialIndex vertexIndex_;
    int lastUpdateStep_ = -1;
    int lastOutputStep_ = 0;
    double lastUpdateTime_s_ = 0.0;
    double updateWindowDuration_s_ = 0.0;
    double outputWindowDuration_s_ = 0.0;
    int resolvedImpactsStep_ = 0;
    int unresolvedImpactsStep_ = 0;
    int ignoredNonIonImpactsStep_ = 0;
    int regionMissImpactsStep_ = 0;
    GeometryMetrics lastGeometryMetrics_;
    std::vector<MatchStatsRow> lastMatchStatsRows_;
    bool linesCurrentGeometryFrame_ = true;
    bool lastUpdateUsedAdvectedPath_ = false;
    bool lastUpdateHadActiveNormalErosion_ = false;
    std::ofstream geometryOut_;
    std::ofstream consistencyOut_;
    std::ofstream regionOut_;
    std::ofstream impactOut_;
    std::ofstream matchOut_;

    enum class RegionMatchKind : unsigned char
    {
        none = 0u,
        direct = 1u,
        fallback = 2u,
    };

    template<typename Mesh>
    void initializeFields_(const Mesh& mesh)
    {
        currentVx_.initialize(mesh, 0.0);
        currentVy_.initialize(mesh, 0.0);
        previousVx_.initialize(mesh, 0.0);
        previousVy_.initialize(mesh, 0.0);
        buildVertexSpatialIndex_(mesh);
    }

    static int clampBinIndex_(int index, int count)
    {
        if (count <= 0)
            return 0;
        if (index < 0)
            return 0;
        if (index >= count)
            return count - 1;
        return index;
    }

    static std::size_t flatBinIndex_(int ix, int iy, int nx)
    {
        return static_cast<std::size_t>(iy * nx + ix);
    }

    void bufferStepSamples_(double dt_s)
    {
        updateWindowSamples_.insert(updateWindowSamples_.end(), pendingSamples_.begin(), pendingSamples_.end());
        outputWindowSamples_.insert(outputWindowSamples_.end(), pendingSamples_.begin(), pendingSamples_.end());
        updateWindowDuration_s_ += std::max(0.0, dt_s);
        outputWindowDuration_s_ += std::max(0.0, dt_s);
    }

    void clearUpdateWindow_()
    {
        updateWindowSamples_.clear();
        updateWindowDuration_s_ = 0.0;
    }

    void clearOutputWindow_()
    {
        outputWindowSamples_.clear();
        outputWindowDuration_s_ = 0.0;
    }

    template<typename Mesh>
    void buildVertexSpatialIndex_(const Mesh& mesh)
    {
        vertexIndex_ = VertexSpatialIndex{};
        if (mesh.vertIndexNum() <= 0)
            return;

        double minX = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();
        for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
        {
            const Vector2D coord = mesh.verts(mesh.vertIndexToLabel[vi]).coord();
            minX = std::min(minX, coord.x);
            maxX = std::max(maxX, coord.x);
            minY = std::min(minY, coord.y);
            maxY = std::max(maxY, coord.y);
        }
        if (!std::isfinite(minX) || !std::isfinite(maxX) || !std::isfinite(minY) || !std::isfinite(maxY))
            return;

        const double spanX = std::max(maxX - minX, 1.0e-12);
        const double spanY = std::max(maxY - minY, 1.0e-12);
        const int target = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(mesh.vertIndexNum()) / 8.0)));
        const double cellSize = std::max(1.0e-12, std::max(spanX, spanY) / static_cast<double>(target));

        vertexIndex_.valid = true;
        vertexIndex_.minX = minX;
        vertexIndex_.minY = minY;
        vertexIndex_.cellSize = cellSize;
        vertexIndex_.nx = std::max(1, static_cast<int>(std::ceil(spanX / cellSize)));
        vertexIndex_.ny = std::max(1, static_cast<int>(std::ceil(spanY / cellSize)));
        vertexIndex_.bins.assign(static_cast<std::size_t>(vertexIndex_.nx * vertexIndex_.ny), {});

        for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
        {
            const Vector2D coord = mesh.verts(mesh.vertIndexToLabel[vi]).coord();
            const int ix = clampBinIndex_(static_cast<int>(std::floor((coord.x - minX) / cellSize)), vertexIndex_.nx);
            const int iy = clampBinIndex_(static_cast<int>(std::floor((coord.y - minY) / cellSize)), vertexIndex_.ny);
            vertexIndex_.bins[flatBinIndex_(ix, iy, vertexIndex_.nx)].push_back(vi);
        }
    }

    template<typename Mesh>
    std::vector<int> collectCandidateVerts_(const Mesh& mesh, const Bounds2D& bounds) const
    {
        std::vector<int> candidates;
        if (!vertexIndex_.valid || !bounds.valid())
        {
            candidates.resize(static_cast<std::size_t>(mesh.vertIndexNum()));
            for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
                candidates[static_cast<std::size_t>(vi)] = vi;
            return candidates;
        }

        const int ix0 = clampBinIndex_(static_cast<int>(std::floor((bounds.minX - vertexIndex_.minX) / vertexIndex_.cellSize)), vertexIndex_.nx);
        const int ix1 = clampBinIndex_(static_cast<int>(std::floor((bounds.maxX - vertexIndex_.minX) / vertexIndex_.cellSize)), vertexIndex_.nx);
        const int iy0 = clampBinIndex_(static_cast<int>(std::floor((bounds.minY - vertexIndex_.minY) / vertexIndex_.cellSize)), vertexIndex_.ny);
        const int iy1 = clampBinIndex_(static_cast<int>(std::floor((bounds.maxY - vertexIndex_.minY) / vertexIndex_.cellSize)), vertexIndex_.ny);
        for (int iy = iy0; iy <= iy1; ++iy)
        {
            for (int ix = ix0; ix <= ix1; ++ix)
            {
                const auto& bin = vertexIndex_.bins[flatBinIndex_(ix, iy, vertexIndex_.nx)];
                candidates.insert(candidates.end(), bin.begin(), bin.end());
            }
        }

        if (candidates.empty())
        {
            candidates.resize(static_cast<std::size_t>(mesh.vertIndexNum()));
            for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
                candidates[static_cast<std::size_t>(vi)] = vi;
        }
        return candidates;
    }

    Bounds2D selectorBounds_(const MovingEBSelector& selector) const
    {
        Bounds2D bounds;
        const auto lengthToCode = [&](double v) {
            return (units_ != nullptr) ? units_->LengthSIToCode(v) : v;
        };
        if (selector.hasXRange)
        {
            bounds.minX = lengthToCode(selector.xMin_m);
            bounds.maxX = lengthToCode(selector.xMax_m);
        }
        if (selector.hasYRange)
        {
            const double y0 = std::abs(lengthToCode(selector.yMin_m));
            const double y1 = std::abs(lengthToCode(selector.yMax_m));
            bounds.minY = -std::max(y0, y1);
            bounds.maxY = std::max(y0, y1);
        }
        return bounds;
    }

    template<typename Mesh>
    Bounds2D regionLineBounds_(std::size_t regionIndex, const ImmersedGeometry<Mesh>& iGeometry) const
    {
        Bounds2D bounds;
        if (regionIndex >= regionCaches_.size())
            return bounds;
        for (int lid : regionCaches_[regionIndex].lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(iGeometry.lines.size()))
                continue;
            const auto& ln = iGeometry.lines[static_cast<std::size_t>(lid)];
            bounds.include(ln.points[0]);
            bounds.include(ln.points[1]);
        }
        return bounds;
    }

    template<typename Mesh>
    Bounds2D regionCandidateBounds_(const MovingEBRegion& region,
                                    std::size_t regionIndex,
                                    const ImmersedGeometry<Mesh>& iGeometry,
                                    double margin) const
    {
        Bounds2D bounds = regionLineBounds_(regionIndex, iGeometry);
        if (bounds.valid())
            bounds.expand(margin);
        bounds.intersectWith(selectorBounds_(region.selector));
        return bounds;
    }

    static void collectImpactCandidateLines_(const ImpactResolveContext& ctx,
                                             const Vector2D& pos,
                                             int boundaryId,
                                             std::vector<int>& out)
    {
        out.clear();
        if (!ctx.enabled || ctx.lineBin.nx <= 0 || ctx.lineBin.ny <= 0 || ctx.lineBin.bins.empty())
            return;

        const int baseIx = clampBinIndex_(static_cast<int>(std::floor((pos.x - ctx.lineBin.minX) / ctx.lineBin.cellSize)), ctx.lineBin.nx);
        const int baseIy = clampBinIndex_(static_cast<int>(std::floor((pos.y - ctx.lineBin.minY) / ctx.lineBin.cellSize)), ctx.lineBin.ny);
        for (int radius = 0; radius <= 1 && out.empty(); ++radius)
        {
            const int ix0 = std::max(0, baseIx - radius);
            const int ix1 = std::min(ctx.lineBin.nx - 1, baseIx + radius);
            const int iy0 = std::max(0, baseIy - radius);
            const int iy1 = std::min(ctx.lineBin.ny - 1, baseIy + radius);
            for (int iy = iy0; iy <= iy1; ++iy)
            {
                for (int ix = ix0; ix <= ix1; ++ix)
                {
                    const auto& bin = ctx.lineBin.bins[flatBinIndex_(ix, iy, ctx.lineBin.nx)];
                    out.insert(out.end(), bin.begin(), bin.end());
                }
            }
        }

        if (boundaryId > 0 && !out.empty())
        {
            const auto it = ctx.lineBin.boundaryLineIds.find(boundaryId);
            if (it != ctx.lineBin.boundaryLineIds.end())
            {
                std::vector<int> filtered;
                filtered.reserve(out.size());
                for (int li : out)
                {
                    if (std::find(it->second.begin(), it->second.end(), li) != it->second.end())
                        filtered.push_back(li);
                }
                out.swap(filtered);
            }
        }

        if (!out.empty())
        {
            std::sort(out.begin(), out.end());
            out.erase(std::unique(out.begin(), out.end()), out.end());
        }
    }

    template<typename Mesh>
    std::vector<Vector2D> collectBoundaryCutPoints_(const Mesh& mesh,
                                                    const ImmersedGeometry<Mesh>& iGeometry,
                                                    int boundaryId) const
    {
        std::vector<Vector2D> pts;
        if (dictGeometry_ == nullptr)
            return pts;
        const auto bit = dictGeometry_->boundaryById.find(boundaryId);
        if (bit == dictGeometry_->boundaryById.end())
            return pts;
        const int levelSetIndex = bit->second.levelSetIndex;
        if (levelSetIndex < 0 || levelSetIndex >= static_cast<int>(iGeometry.levelSets.size()))
            return pts;

        const auto addEdgePoint = [&](const Edge& edge) {
            if (levelSetIndex >= static_cast<int>(edge.isCut.size()) || !edge.isCut[static_cast<std::size_t>(levelSetIndex)])
                return;
            const Vert& v0 = mesh.verts(edge.edgeToVert[0]);
            const Vert& v1 = mesh.verts(edge.edgeToVert[1]);
            const double psi0 = iGeometry.levelSets[levelSetIndex].value[v0.index()];
            const double psi1 = iGeometry.levelSets[levelSetIndex].value[v1.index()];
            if (!std::isfinite(psi0) || !std::isfinite(psi1))
                return;
            double t = 0.5;
            const double denom = psi0 - psi1;
            if (std::abs(denom) > 1.0e-14)
                t = psi0 / denom;
            t = std::clamp(t, 0.0, 1.0);
            const Vector2D p0 = v0.coord();
            const Vector2D p1 = v1.coord();
            pts.push_back(p0 + (p1 - p0) * t);
        };

        for (int ei = 0; ei < mesh.hEdgeIndexNum(); ++ei)
            addEdgePoint(mesh.hEdges(mesh.hEdgeIndexToLabel[ei]));
        for (int ei = 0; ei < mesh.vEdgeIndexNum(); ++ei)
            addEdgePoint(mesh.vEdges(mesh.vEdgeIndexToLabel[ei]));
        return pts;
    }

    template<typename Mesh>
    std::vector<Vector2D> collectBoundaryPolylinePoints_(const ImmersedGeometry<Mesh>& iGeometry,
                                                         int boundaryId) const
    {
        std::vector<Vector2D> pts;
        if (dictGeometry_ == nullptr)
            return pts;
        auto it = dictGeometry_->boundaryById.find(boundaryId);
        if (it == dictGeometry_->boundaryById.end())
            return pts;
        for (int lid : it->second.lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(iGeometry.lines.size()))
                continue;
            const auto& ln = iGeometry.lines[static_cast<std::size_t>(lid)];
            pts.push_back(ln.points[0]);
            pts.push_back(ln.points[1]);
        }
        return pts;
    }

    template<typename Mesh>
    std::vector<Vector2D> collectBoundaryMetricPoints_(const Mesh& mesh,
                                                       const ImmersedGeometry<Mesh>& iGeometry,
                                                       int boundaryId) const
    {
        std::vector<Vector2D> pts = collectBoundaryCutPoints_(mesh, iGeometry, boundaryId);
        const std::vector<Vector2D> polyPts = collectBoundaryPolylinePoints_(iGeometry, boundaryId);
        if (!pts.empty() && !polyPts.empty())
        {
            double minPx = std::numeric_limits<double>::infinity();
            double maxPx = -std::numeric_limits<double>::infinity();
            double minPy = std::numeric_limits<double>::infinity();
            double maxPy = -std::numeric_limits<double>::infinity();
            for (const auto& p : polyPts)
            {
                minPx = std::min(minPx, p.x);
                maxPx = std::max(maxPx, p.x);
                minPy = std::min(minPy, p.y);
                maxPy = std::max(maxPy, p.y);
            }

            const double span = std::max(maxPx - minPx, maxPy - minPy);
            const double margin = std::max(5.0e-3, 0.5 * span);
            std::vector<Vector2D> filtered;
            filtered.reserve(pts.size());
            for (const auto& p : pts)
            {
                if (p.x < (minPx - margin) || p.x > (maxPx + margin))
                    continue;
                if (p.y < (minPy - margin) || p.y > (maxPy + margin))
                    continue;
                filtered.push_back(p);
            }
            if (!filtered.empty())
                return filtered;
        }
        if (!pts.empty())
            return pts;
        return polyPts;
    }

    template<typename Mesh>
    GeometryMetrics computeGeometryMetrics_(const Mesh& mesh, const ImmersedGeometry<Mesh>& iGeometry) const
    {
        const auto codeLengthToSI = [&](double v) {
            return (units_ != nullptr) ? units_->LengthCodeToSI(v) : v;
        };
        GeometryMetrics gm;
        auto minX = [](const std::vector<Vector2D>& pts) {
            double out = std::numeric_limits<double>::infinity();
            for (const auto& p : pts)
                out = std::min(out, p.x);
            return std::isfinite(out) ? out : 0.0;
        };
        auto maxX = [](const std::vector<Vector2D>& pts) {
            double out = -std::numeric_limits<double>::infinity();
            for (const auto& p : pts)
                out = std::max(out, p.x);
            return std::isfinite(out) ? out : 0.0;
        };
        auto minAbsY = [](const std::vector<Vector2D>& pts) {
            double out = std::numeric_limits<double>::infinity();
            for (const auto& p : pts)
                out = std::min(out, std::abs(p.y));
            return std::isfinite(out) ? out : 0.0;
        };

        auto screen = collectBoundaryMetricPoints_(mesh, iGeometry, 10);
        auto lowerScreen = collectBoundaryMetricPoints_(mesh, iGeometry, 11);
        auto accel = collectBoundaryMetricPoints_(mesh, iGeometry, 12);
        auto lowerAccel = collectBoundaryMetricPoints_(mesh, iGeometry, 13);
        auto decel = collectBoundaryMetricPoints_(mesh, iGeometry, 14);
        auto lowerDecel = collectBoundaryMetricPoints_(mesh, iGeometry, 15);
        screen.insert(screen.end(), lowerScreen.begin(), lowerScreen.end());
        accel.insert(accel.end(), lowerAccel.begin(), lowerAccel.end());
        decel.insert(decel.end(), lowerDecel.begin(), lowerDecel.end());
        gm.gsa_mm = 1.0e3 * codeLengthToSI(minX(accel) - maxX(screen));
        gm.gad_mm = 1.0e3 * codeLengthToSI(minX(decel) - maxX(accel));
        gm.da_mm = 2.0e3 * codeLengthToSI(minAbsY(accel));
        return gm;
    }

    template<typename Mesh>
    void accumulateRegionStats_(const std::vector<WallImpactSample>& samples,
                                int step,
                                int windowStartExclusive,
                                double dt_s,
                                const ImmersedGeometry<Mesh>& iGeometry)
    {
        const auto& cfg = config_;
        const int windowStepCount = std::max(1, step - std::max(windowStartExclusive, 0));
        std::vector<unsigned char> matchedToActiveErosion(samples.size(), 0u);
        std::vector<int> assignedActiveRegion(samples.size(), -1);
        std::vector<RegionMatchKind> assignedMatchKind(samples.size(), RegionMatchKind::none);
        std::map<std::pair<int, int>, MatchStatsRow> matchBuckets;
        bool anyActiveErosion = false;
        for (std::size_t si = 0; si < samples.size(); ++si)
        {
            const auto& sample = samples[si];
            if (!sample.resolved || !sample.ionSpecies)
                continue;
            RegionMatchKind matchKind = RegionMatchKind::none;
            assignedActiveRegion[si] = matchActiveErosionRegion_(sample, step, iGeometry, &matchKind);
            assignedMatchKind[si] = matchKind;
            if (assignedActiveRegion[si] >= 0)
                matchedToActiveErosion[si] = 1u;

            auto& bucket = matchBuckets[std::make_pair(sample.boundaryId, sample.levelSetIndex)];
            bucket.boundaryId = sample.boundaryId;
            bucket.levelSetIndex = sample.levelSetIndex;
            ++bucket.resolvedIonCount;
            if (assignedActiveRegion[si] >= 0)
            {
                ++bucket.matchedResolvedIonCount;
                if (matchKind == RegionMatchKind::direct)
                    ++bucket.directMatchCount;
                else if (matchKind == RegionMatchKind::fallback)
                    ++bucket.fallbackMatchCount;
            }
            else
            {
                ++bucket.unmatchedResolvedIonCount;
            }
        }
        for (std::size_t ri = 0; ri < cfg.regions.size(); ++ri)
        {
            const auto& region = cfg.regions[ri];
            auto& state = states_[ri];
            state.lastArea_m2 = computeRegionArea_(ri, region, iGeometry);
            state.activeThisStep = regionIsActive_(region, step);
            state.activeWindowStepCount = activeStepCountInWindow_(region, windowStartExclusive, step);
            state.activeInWindow = (state.activeWindowStepCount > 0);
            state.lastResolvedIonCount = 0;
            state.lastDirectMatchCount = 0;
            state.lastFallbackMatchCount = 0;
            double count = 0.0;
            double eWeighted = 0.0;
            double cosWeighted = 0.0;
            const bool activeNormalErosion = state.activeInWindow && (regionCaches_[ri].motionMode == "normalerosion");
            anyActiveErosion = anyActiveErosion || activeNormalErosion;
            for (std::size_t si = 0; si < samples.size(); ++si)
            {
                const auto& sample = samples[si];
                if (!sample.resolved)
                    continue;
                if (!state.activeInWindow)
                    continue;
                if (assignedActiveRegion[si] != static_cast<int>(ri))
                    continue;
                ++state.lastResolvedIonCount;
                if (assignedMatchKind[si] == RegionMatchKind::direct)
                    ++state.lastDirectMatchCount;
                else if (assignedMatchKind[si] == RegionMatchKind::fallback)
                    ++state.lastFallbackMatchCount;
                count += std::max(0.0, sample.macroWeight);
                eWeighted += std::max(0.0, sample.macroWeight) * std::max(0.0, sample.energy_eV);
                cosWeighted += std::max(0.0, sample.macroWeight) * std::max(0.0, std::cos(sample.theta_rad));
            }
            state.lastSampleCount = count;
            const double duration =
                std::max(1e-30,
                         dt_s * (static_cast<double>(state.activeWindowStepCount) / static_cast<double>(windowStepCount)));
            const double flux = (state.lastArea_m2 > 0.0) ? (count / (state.lastArea_m2 * duration)) : 0.0;
            const double meanEnergy = (count > 0.0) ? (eWeighted / count) : state.smoothedEnergy;
            const double meanCos = (count > 0.0) ? (cosWeighted / count) : state.smoothedCos;
            if (toLower_(cfg.impactModel.smoothing) != "none" && lastUpdateStep_ >= 0)
            {
                state.smoothedFlux = 0.5 * state.smoothedFlux + 0.5 * flux;
                state.smoothedEnergy = 0.5 * state.smoothedEnergy + 0.5 * meanEnergy;
                state.smoothedCos = 0.5 * state.smoothedCos + 0.5 * meanCos;
            }
            else
            {
                state.smoothedFlux = flux;
                state.smoothedEnergy = meanEnergy;
                state.smoothedCos = meanCos;
            }
        }
        regionMissImpactsStep_ = 0;
        if (anyActiveErosion)
        {
            for (std::size_t si = 0; si < samples.size(); ++si)
            {
                const auto& sample = samples[si];
                if (sample.resolved &&
                    sample.ionSpecies &&
                    !matchedToActiveErosion[si] &&
                    sampleTouchesActiveErosionSupport_(sample, step, iGeometry))
                    ++regionMissImpactsStep_;
            }
        }
        lastMatchStatsRows_.clear();
        for (const auto& item : matchBuckets)
            lastMatchStatsRows_.push_back(item.second);
    }

    template<typename Mesh>
    double computeRegionArea_(std::size_t regionIndex,
                              const MovingEBRegion& region,
                              const ImmersedGeometry<Mesh>& iGeometry) const
    {
        const double slab = (normalization_ != nullptr && normalization_->config().slabThickness_m > 0.0)
                                ? normalization_->config().slabThickness_m
                                : 1.0;
        double length = 0.0;
        for (int lid : regionCaches_[regionIndex].lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(iGeometry.lines.size()))
                continue;
            const auto& ln = iGeometry.lines[static_cast<std::size_t>(lid)];
            const Vector2D mid = (ln.points[0] + ln.points[1]) * 0.5;
            if (!selectorMatchesPoint_(mid, region.selector))
                continue;
            length += (ln.points[1] - ln.points[0]).mod();
        }
        const double lengthSI = (units_ != nullptr) ? units_->LengthCodeToSI(length) : length;
        return std::max(0.0, lengthSI * slab);
    }

    MotionSummary updateRegionMotionState_(int step, int windowStartExclusive)
    {
        const double codePerMeter = (units_ != nullptr) ? units_->LengthSIToCode(1.0) : 1.0;
        MotionSummary summary;
        lastUpdateHadActiveNormalErosion_ = false;
        for (std::size_t ri = 0; ri < config_.regions.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            const auto& cache = regionCaches_[ri];
            auto& state = states_[ri];
            state.activeThisStep = regionIsActive_(region, step) && cache.motionMode != "fixed";
            state.activeWindowStepCount = (cache.motionMode != "fixed")
                                              ? activeStepCountInWindow_(region, windowStartExclusive, step)
                                              : 0;
            state.activeInWindow = (state.activeWindowStepCount > 0);
            if (!state.activeInWindow)
            {
                state.lastNormalSpeed_mps = 0.0;
                state.lastRigidVx_codeps = 0.0;
                state.lastRigidVy_codeps = 0.0;
                continue;
            }
            summary.anyActive = true;
            if (cache.motionMode == "normalerosion")
                lastUpdateHadActiveNormalErosion_ = true;
            Vector2D rigidVel(region.velocityX_mps * codePerMeter, region.velocityY_mps * codePerMeter);
            state.lastRigidVx_codeps = rigidVel.x;
            state.lastRigidVy_codeps = rigidVel.y;
            double normalSpeed = 0.0;
            if (cache.model == "ionerosion")
                normalSpeed = computeIonErosionSpeed_(region, state);
            else if (cache.model == "prescribed")
                normalSpeed = std::max(0.0, region.normalSpeed_mps);
            state.lastNormalSpeed_mps = normalSpeed;
            const double normalSpeedCode = normalSpeed * codePerMeter;
            if (rigidVel.mod2() > 0.0 || normalSpeedCode > 0.0)
                summary.anyMotion = true;
        }
        return summary;
    }

    template<typename Mesh>
    void buildVelocityFields_(const Mesh& mesh, const ImmersedGeometry<Mesh>& iGeometry, int windowStepCount)
    {
        previousVx_ = currentVx_;
        previousVy_ = currentVy_;
        currentVx_.initialize(mesh, 0.0);
        currentVy_.initialize(mesh, 0.0);

        for (std::size_t ri = 0; ri < config_.regions.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            const auto& cache = regionCaches_[ri];
            auto& state = states_[ri];
            if (!state.activeInWindow || state.activeWindowStepCount <= 0)
                continue;

            const std::string& motionMode = cache.motionMode;
            const double activityScale =
                static_cast<double>(state.activeWindowStepCount) / static_cast<double>(std::max(1, windowStepCount));
            const Vector2D rigidVel(state.lastRigidVx_codeps * activityScale, state.lastRigidVy_codeps * activityScale);
            const double normalSpeedCode =
                state.lastNormalSpeed_mps * ((units_ != nullptr) ? units_->LengthSIToCode(1.0) : 1.0) * activityScale;

            const double bandSI = (region.smoothingLength_m > 0.0)
                                      ? region.smoothingLength_m
                                      : std::max(1e-8, config_.impactModel.smoothingLength_m);
            const double band = std::max(1e-8, (units_ != nullptr) ? units_->LengthSIToCode(bandSI) : bandSI);
            const bool useCurrentLevelSetBand =
                (region.levelSetIndex >= 0 && region.levelSetIndex < static_cast<int>(iGeometry.levelSets.size()) &&
                 (region.selector.wholeBoundary || (region.selector.lineIds.empty() && region.boundaryId > 0)));
            const Bounds2D candidateBounds = regionCandidateBounds_(region, ri, iGeometry, 4.0 * band);
            const std::vector<int> candidateVerts = collectCandidateVerts_(mesh, candidateBounds);
            for (int vi : candidateVerts)
            {
                const Vector2D coord = mesh.verts(mesh.vertIndexToLabel[vi]).coord();
                if (!selectorMatchesPoint_(coord, region.selector))
                    continue;
                double weight = 1.0;
                double d2 = 0.0;
                if (motionMode == "rigidshift" && useCurrentLevelSetBand)
                {
                    const double phi = std::abs(iGeometry.levelSets[region.levelSetIndex].value[vi]);
                    d2 = phi * phi;
                }
                else
                {
                    d2 = distanceToRegionLines2_(coord, ri, iGeometry);
                }
                if (d2 > 0.0)
                {
                    const double d = std::sqrt(d2);
                    if (d > 4.0 * band)
                        continue;
                    weight = std::exp(-(d * d) / std::max(1e-30, band * band));
                }

                currentVx_.value[vi] += weight * rigidVel.x;
                currentVy_.value[vi] += weight * rigidVel.y;
                if (motionMode == "normalerosion" &&
                    normalSpeedCode > 0.0 &&
                    region.levelSetIndex >= 0 &&
                    region.levelSetIndex < static_cast<int>(iGeometry.levelSets.size()))
                {
                    const double nx = state.normalDirectionSign * iGeometry.LS_dx(region.levelSetIndex).value[vi];
                    const double ny = state.normalDirectionSign * iGeometry.LS_dy(region.levelSetIndex).value[vi];
                    const double nmag = std::sqrt(nx * nx + ny * ny);
                    if (nmag > 0.0)
                    {
                        currentVx_.value[vi] += weight * normalSpeedCode * nx / nmag;
                        currentVy_.value[vi] += weight * normalSpeedCode * ny / nmag;
                    }
                }
            }
        }
    }

    bool canUseExactRigidShiftPath_(int windowStartExclusive, int windowEndInclusive) const
    {
        if (!dictGeometry_)
            return false;

        bool anyActive = false;
        for (std::size_t ri = 0; ri < config_.regions.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            const auto& cache = regionCaches_[ri];
            if (cache.motionMode == "fixed")
                continue;
            if (activeStepCountInWindow_(region, windowStartExclusive, windowEndInclusive) <= 0)
                continue;
            anyActive = true;
            if (cache.motionMode != "rigidshift")
                return false;
            if (!region.selector.wholeBoundary)
                return false;
            if (region.boundaryId <= 0 || region.levelSetIndex < 0)
                return false;
        }
        return anyActive;
    }

    template<typename Mesh>
    void applyExactRigidShiftGeometry_(ImmersedGeometry<Mesh>& iGeometry) const
    {
        if (!dictGeometry_)
            return;

        iGeometry.lines = dictGeometry_->lines;
        for (std::size_t lsI = 0; lsI < std::min(iGeometry.levelSets.size(), dictGeometry_->innerPoints.size()); ++lsI)
            iGeometry.levelSets[lsI].setInnerPoint(dictGeometry_->innerPoints[lsI]);

        for (std::size_t ri = 0; ri < config_.regions.size() && ri < states_.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            if (toLower_(region.motionMode) != "rigidshift" || region.boundaryId <= 0)
                continue;

            auto bit = dictGeometry_->boundaryById.find(region.boundaryId);
            if (bit == dictGeometry_->boundaryById.end())
                continue;

            const Vector2D shift(states_[ri].rigidShiftX_code, states_[ri].rigidShiftY_code);
            for (int lid : bit->second.lineIds)
            {
                if (lid < 0 || lid >= static_cast<int>(iGeometry.lines.size()))
                    continue;
                auto& ln = iGeometry.lines[static_cast<std::size_t>(lid)];
                ln.points[0] += shift;
                ln.points[1] += shift;
                if (ln.type() == LineType::arc)
                    ln.setArcCenter(ln.arcCenter() + shift);
            }

            if (region.levelSetIndex >= 0 && region.levelSetIndex < static_cast<int>(iGeometry.levelSets.size()) &&
                region.levelSetIndex < static_cast<int>(dictGeometry_->innerPoints.size()))
            {
                iGeometry.levelSets[region.levelSetIndex].setInnerPoint(
                    dictGeometry_->innerPoints[static_cast<std::size_t>(region.levelSetIndex)] + shift);
            }
        }
    }

    template<typename Mesh>
    void initializeNormalDirections_(const Mesh& mesh, const ImmersedGeometry<Mesh>& iGeometry)
    {
        for (std::size_t ri = 0; ri < config_.regions.size() && ri < states_.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            auto& state = states_[ri];
            if (regionCaches_[ri].motionMode != "normalerosion")
                continue;
            if (region.levelSetIndex < 0 || region.levelSetIndex >= static_cast<int>(iGeometry.levelSets.size()))
                continue;
            if (state.normalDirectionInitialized)
                continue;

            double sign = -1.0;
            const Vector2D innerPoint = iGeometry.levelSets[region.levelSetIndex].innerPoint();
            try
            {
                int cellLabel = mesh.findCellLabel(innerPoint);
                if (cellLabel < 0)
                    cellLabel = 0;
                const Cell& cell = mesh.cells(cellLabel);
                const double phi = Interpolation::biLinearInterpolation(mesh, iGeometry.levelSets[region.levelSetIndex], cell, innerPoint);
                if (std::isfinite(phi) && std::abs(phi) > 1e-30)
                    sign = (phi < 0.0) ? -1.0 : 1.0;
            }
            catch (...) {}
            state.normalDirectionSign = sign;
            state.normalDirectionInitialized = true;
        }
    }

    double computeIonErosionSpeed_(const MovingEBRegion& region, const RegionState& state) const
    {
        // Gate erosion on the number of resolved ion hits, not on macro-weight.
        // Using weighted count here makes a single macro-particle look like
        // millions of samples and defeats the intended anti-noise threshold.
        if (state.lastResolvedIonCount < std::max(1, config_.impactModel.minSamples))
            return 0.0;
        double energyFactor = 1.0;
        const std::string& energyModel = regionCaches_[&region - config_.regions.data()].energyModel;
        if (energyModel == "linear")
            energyFactor = std::max(0.0, state.smoothedEnergy / std::max(region.energyRef_eV, 1e-12));
        else if (energyModel == "sqrt")
            energyFactor = std::sqrt(std::max(0.0, state.smoothedEnergy / std::max(region.energyRef_eV, 1e-12)));

        double angleFactor = 1.0;
        if (regionCaches_[&region - config_.regions.data()].angleModel == "cosine")
            angleFactor = std::pow(std::max(0.0, state.smoothedCos), std::max(0.0, region.angleExponent));

        double speed = region.C0 * std::pow(std::max(0.0, state.smoothedFlux), std::max(0.0, region.alpha)) * energyFactor * angleFactor;
        if (region.vmax_mps > 0.0)
            speed = std::min(speed, region.vmax_mps);
        return std::max(0.0, speed);
    }

    bool sampleMatchesRegion_(const WallImpactSample& sample, const MovingEBRegion& region) const
    {
        if (region.boundaryId > 0 && sample.boundaryId != region.boundaryId)
            return false;
        if (!region.selector.lineIds.empty())
        {
            if (sample.lineIndex < 0)
                return false;
            if (std::find(region.selector.lineIds.begin(), region.selector.lineIds.end(), sample.lineIndex) == region.selector.lineIds.end())
                return false;
        }
        if (region.levelSetIndex >= 0 && sample.levelSetIndex != region.levelSetIndex)
            return false;
        return selectorMatchesPoint_(sample.pos, region.selector);
    }

    double erosionFallbackLineTol2_() const
    {
        if (dictGeometry_ == nullptr)
            return 0.0;
        const double tolSI =
            std::max(5.0e-6, 0.25 * std::max(0.0, config_.impactModel.smoothingLength_m));
        const double tolCode = (units_ != nullptr) ? units_->LengthSIToCode(tolSI) : tolSI;
        return tolCode * tolCode;
    }

    template<typename Mesh>
    int matchActiveErosionRegion_(const WallImpactSample& sample,
                                  int step,
                                  const ImmersedGeometry<Mesh>& iGeometry,
                                  RegionMatchKind* matchKindOut = nullptr) const
    {
        if (matchKindOut != nullptr)
            *matchKindOut = RegionMatchKind::none;
        int bestRegion = -1;
        double bestD2 = std::numeric_limits<double>::infinity();
        const double fallbackTol2 = erosionFallbackLineTol2_();

        for (std::size_t ri = 0; ri < config_.regions.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            if (!regionIsActive_(region, step))
                continue;
            if (regionCaches_[ri].motionMode != "normalerosion")
                continue;
            if (region.boundaryId > 0 && sample.boundaryId != region.boundaryId)
                continue;
            if (region.levelSetIndex >= 0 && sample.levelSetIndex != region.levelSetIndex)
                continue;
            if (!selectorMatchesPoint_(sample.pos, region.selector))
                continue;
            if (sampleMatchesRegion_(sample, region))
            {
                if (matchKindOut != nullptr)
                    *matchKindOut = RegionMatchKind::direct;
                return static_cast<int>(ri);
            }
            if (region.selector.lineIds.empty())
                continue;

            const double d2 = distanceToRegionLines2_(sample.pos, ri, iGeometry);
            if (d2 < bestD2)
            {
                bestD2 = d2;
                bestRegion = static_cast<int>(ri);
            }
        }

        if (bestRegion >= 0 && std::isfinite(bestD2) && bestD2 <= fallbackTol2)
        {
            if (matchKindOut != nullptr)
                *matchKindOut = RegionMatchKind::fallback;
            return bestRegion;
        }
        return -1;
    }

    template<typename Mesh>
    bool sampleTouchesActiveErosionSupport_(const WallImpactSample& sample,
                                            int step,
                                            const ImmersedGeometry<Mesh>& iGeometry) const
    {
        if (dictGeometry_ == nullptr)
            return false;
        const double fallbackTol2 = erosionFallbackLineTol2_();
        for (std::size_t ri = 0; ri < config_.regions.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            if (!regionIsActive_(region, step))
                continue;
            if (regionCaches_[ri].motionMode != "normalerosion")
                continue;
            if (region.boundaryId > 0 && sample.boundaryId != region.boundaryId)
                continue;
            if (region.levelSetIndex >= 0 && sample.levelSetIndex != region.levelSetIndex)
                continue;
            if (!selectorMatchesPoint_(sample.pos, region.selector))
                continue;
            if (region.selector.lineIds.empty() || sampleMatchesRegion_(sample, region))
                return true;
            const double d2 = distanceToRegionLines2_(sample.pos, ri, iGeometry);
            if (std::isfinite(d2) && d2 <= fallbackTol2)
                return true;
        }
        return false;
    }

    bool regionIsActive_(const MovingEBRegion& region, int step) const
    {
        if (step < region.startStep)
            return false;
        if (region.endStep >= 0 && step > region.endStep)
            return false;
        return true;
    }

    int activeStepCountInWindow_(const MovingEBRegion& region,
                                 int windowStartExclusive,
                                 int windowEndInclusive) const
    {
        if (windowEndInclusive <= windowStartExclusive)
            return 0;
        const int firstActiveStep = std::max(region.startStep, std::max(windowStartExclusive, 0) + 1);
        int lastActiveStep = windowEndInclusive;
        if (region.endStep >= 0)
            lastActiveStep = std::min(lastActiveStep, region.endStep);
        if (lastActiveStep < firstActiveStep)
            return 0;
        return lastActiveStep - firstActiveStep + 1;
    }

    bool selectorMatchesPoint_(const Vector2D& p, const MovingEBSelector& selector) const
    {
        const auto lengthToCode = [&](double v) {
            return (units_ != nullptr) ? units_->LengthSIToCode(v) : v;
        };
        if (selector.hasXRange && (p.x < lengthToCode(selector.xMin_m) || p.x > lengthToCode(selector.xMax_m)))
            return false;
        if (selector.hasYRange)
        {
            const double ay = std::abs(p.y);
            const double ymin = std::min(std::abs(lengthToCode(selector.yMin_m)), std::abs(lengthToCode(selector.yMax_m)));
            const double ymax = std::max(std::abs(lengthToCode(selector.yMin_m)), std::abs(lengthToCode(selector.yMax_m)));
            if (ay < ymin || ay > ymax)
                return false;
        }
        return true;
    }

    template<typename Mesh>
    double distanceToRegionLines2_(const Vector2D& p, std::size_t regionIndex, const ImmersedGeometry<Mesh>& iGeometry) const
    {
        auto pointSegDist2 = [](const Vector2D& pt, const Vector2D& a, const Vector2D& b) {
            const Vector2D ab = b - a;
            const double denom = ab.mod2();
            if (!(denom > 0.0))
                return (pt - a).mod2();
            const Vector2D ap = pt - a;
            const double proj = ap.x * ab.x + ap.y * ab.y;
            const double t = std::clamp(proj / denom, 0.0, 1.0);
            return (pt - (a + ab * t)).mod2();
        };

        if (regionIndex >= regionCaches_.size() || regionCaches_[regionIndex].lineIds.empty())
            return 0.0;

        double best = std::numeric_limits<double>::infinity();
        for (int lid : regionCaches_[regionIndex].lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(iGeometry.lines.size()))
                continue;
            const auto& ln = iGeometry.lines[static_cast<std::size_t>(lid)];
            best = std::min(best, pointSegDist2(p, ln.points[0], ln.points[1]));
        }
        return std::isfinite(best) ? best : 0.0;
    }

    void openOutputs_(bool append)
    {
        openOutput_(geometryOut_, config_.metricsPath,
                    "step,time_s,updated,resolved_impacts,unresolved_impacts,gsa_mm,gad_mm,da_mm\n", append);
        openOutput_(consistencyOut_, config_.consistencyPath,
                    "step,time_s,updated,lines_current_geometry_frame,used_advected_path,point_count,inner_point_count,"
                    "max_line_point_delta_um,rms_line_point_delta_um,max_arc_center_delta_um,rms_arc_center_delta_um,"
                    "max_inner_point_delta_um,rms_inner_point_delta_um,gsa_live_mm,gsa_ref_mm,gad_live_mm,gad_ref_mm,"
                    "da_live_mm,da_ref_mm\n", append);
        openOutput_(regionOut_, config_.erosionMetricsPath,
                    "step,time_s,updated,region,boundary_id,level_set_index,sample_count,resolved_ion_count,direct_match_count,fallback_match_count,area_m2,flux_1pm2s,mean_energy_eV,mean_cos_theta,vn_mps,normal_dir_sign,normal_dir_initialized,active_flag\n", append);
        openOutput_(impactOut_, config_.impactStatsPath,
                    "step,time_s,updated,resolved_impacts,unresolved_impacts,ignored_nonion_impacts,region_miss_impacts,total_buffered\n", append);
        std::filesystem::path matchPath(config_.erosionMetricsPath);
        matchPath.replace_filename("erosion_match_stats.csv");
        openOutput_(matchOut_, matchPath.string(),
                    "step,time_s,updated,boundary_id,level_set_index,resolved_ion_count,matched_resolved_ion_count,unmatched_resolved_ion_count,direct_match_count,fallback_match_count\n", append);
    }

    static void openOutput_(std::ofstream& out, const std::string& path, const std::string& header, bool append)
    {
        if (path.empty())
            return;
        std::filesystem::path p(path);
        if (!p.parent_path().empty())
            std::filesystem::create_directories(p.parent_path());
        std::error_code ec;
        const bool haveExistingData =
            append && std::filesystem::exists(p, ec) && !ec && std::filesystem::is_regular_file(p, ec) && !ec &&
            (std::filesystem::file_size(p, ec) > 0u);
        if (haveExistingData)
            out.open(p.string(), std::ios::out | std::ios::app);
        else
            out.open(p.string(), std::ios::out | std::ios::trunc);
        if (out)
        {
            if (!haveExistingData)
                out << header;
            out << std::setprecision(16) << std::fixed;
        }
    }

    void writeGeometryMetrics_(int step, double time_s, bool updated)
    {
        if (!geometryOut_)
            return;
        geometryOut_ << step << ',' << time_s << ',' << (updated ? 1 : 0) << ','
                     << resolvedImpactsStep_ << ',' << unresolvedImpactsStep_ << ','
                     << lastGeometryMetrics_.gsa_mm << ',' << lastGeometryMetrics_.gad_mm << ','
                     << lastGeometryMetrics_.da_mm << "\n";
        geometryOut_.flush();
    }

    void buildExactRigidShiftReference_(std::vector<Line>& refLines,
                                        std::vector<Vector2D>& refInnerPoints) const
    {
        refLines = dictGeometry_->lines;
        refInnerPoints = dictGeometry_->innerPoints;

        for (std::size_t ri = 0; ri < config_.regions.size() && ri < states_.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            if (toLower_(region.motionMode) != "rigidshift" || region.boundaryId <= 0)
                continue;

            auto bit = dictGeometry_->boundaryById.find(region.boundaryId);
            if (bit == dictGeometry_->boundaryById.end())
                continue;

            const Vector2D shift(states_[ri].rigidShiftX_code, states_[ri].rigidShiftY_code);
            for (int lid : bit->second.lineIds)
            {
                if (lid < 0 || lid >= static_cast<int>(refLines.size()))
                    continue;
                auto& ln = refLines[static_cast<std::size_t>(lid)];
                ln.points[0] += shift;
                ln.points[1] += shift;
                if (ln.type() == LineType::arc)
                    ln.setArcCenter(ln.arcCenter() + shift);
            }

            if (region.levelSetIndex >= 0 &&
                region.levelSetIndex < static_cast<int>(refInnerPoints.size()) &&
                region.levelSetIndex < static_cast<int>(dictGeometry_->innerPoints.size()))
            {
                refInnerPoints[static_cast<std::size_t>(region.levelSetIndex)] =
                    dictGeometry_->innerPoints[static_cast<std::size_t>(region.levelSetIndex)] + shift;
            }
        }
    }

    std::vector<Vector2D> collectBoundaryPolylinePointsFromLines_(const std::vector<Line>& lines,
                                                                  int boundaryId) const
    {
        std::vector<Vector2D> pts;
        if (dictGeometry_ == nullptr)
            return pts;
        auto it = dictGeometry_->boundaryById.find(boundaryId);
        if (it == dictGeometry_->boundaryById.end())
            return pts;
        for (int lid : it->second.lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(lines.size()))
                continue;
            const auto& ln = lines[static_cast<std::size_t>(lid)];
            pts.push_back(ln.points[0]);
            pts.push_back(ln.points[1]);
        }
        return pts;
    }

    GeometryMetrics computeGeometryMetricsFromLines_(const std::vector<Line>& lines) const
    {
        const auto codeLengthToSI = [&](double v) {
            return (units_ != nullptr) ? units_->LengthCodeToSI(v) : v;
        };
        GeometryMetrics gm;
        auto minX = [](const std::vector<Vector2D>& pts) {
            double out = std::numeric_limits<double>::infinity();
            for (const auto& p : pts)
                out = std::min(out, p.x);
            return std::isfinite(out) ? out : 0.0;
        };
        auto maxX = [](const std::vector<Vector2D>& pts) {
            double out = -std::numeric_limits<double>::infinity();
            for (const auto& p : pts)
                out = std::max(out, p.x);
            return std::isfinite(out) ? out : 0.0;
        };
        auto minAbsY = [](const std::vector<Vector2D>& pts) {
            double out = std::numeric_limits<double>::infinity();
            for (const auto& p : pts)
                out = std::min(out, std::abs(p.y));
            return std::isfinite(out) ? out : 0.0;
        };

        auto screen = collectBoundaryPolylinePointsFromLines_(lines, 10);
        auto lowerScreen = collectBoundaryPolylinePointsFromLines_(lines, 11);
        auto accel = collectBoundaryPolylinePointsFromLines_(lines, 12);
        auto lowerAccel = collectBoundaryPolylinePointsFromLines_(lines, 13);
        auto decel = collectBoundaryPolylinePointsFromLines_(lines, 14);
        auto lowerDecel = collectBoundaryPolylinePointsFromLines_(lines, 15);
        screen.insert(screen.end(), lowerScreen.begin(), lowerScreen.end());
        accel.insert(accel.end(), lowerAccel.begin(), lowerAccel.end());
        decel.insert(decel.end(), lowerDecel.begin(), lowerDecel.end());
        gm.gsa_mm = 1.0e3 * codeLengthToSI(minX(accel) - maxX(screen));
        gm.gad_mm = 1.0e3 * codeLengthToSI(minX(decel) - maxX(accel));
        gm.da_mm = 2.0e3 * codeLengthToSI(minAbsY(accel));
        return gm;
    }

    template<typename Mesh>
    void writeConsistencyMetrics_(int step, double time_s, bool updated, const ImmersedGeometry<Mesh>& iGeometry)
    {
        if (!consistencyOut_ || dictGeometry_ == nullptr)
            return;

        std::vector<Line> refLines;
        std::vector<Vector2D> refInnerPoints;
        buildExactRigidShiftReference_(refLines, refInnerPoints);

        double maxPointDelta_um = 0.0;
        double sumPointDelta2_um2 = 0.0;
        std::size_t pointCount = 0;
        double maxArcDelta_um = 0.0;
        double sumArcDelta2_um2 = 0.0;
        std::size_t arcCount = 0;

        const auto codeLengthToUm = [&](double v) {
            const double si = (units_ != nullptr) ? units_->LengthCodeToSI(v) : v;
            return 1.0e6 * si;
        };

        const std::size_t lineCount = std::min(iGeometry.lines.size(), refLines.size());
        for (std::size_t li = 0; li < lineCount; ++li)
        {
            const auto& live = iGeometry.lines[li];
            const auto& ref = refLines[li];
            for (int pi = 0; pi < 2; ++pi)
            {
                const double deltaUm = codeLengthToUm((live.points[pi] - ref.points[pi]).mod());
                maxPointDelta_um = std::max(maxPointDelta_um, deltaUm);
                sumPointDelta2_um2 += deltaUm * deltaUm;
                ++pointCount;
            }
            if (live.type() == LineType::arc && ref.type() == LineType::arc)
            {
                const double deltaUm = codeLengthToUm((live.arcCenter() - ref.arcCenter()).mod());
                maxArcDelta_um = std::max(maxArcDelta_um, deltaUm);
                sumArcDelta2_um2 += deltaUm * deltaUm;
                ++arcCount;
            }
        }

        double maxInnerDelta_um = 0.0;
        double sumInnerDelta2_um2 = 0.0;
        std::size_t innerCount = 0;
        const std::size_t innerPointCount = std::min(iGeometry.levelSets.size(), refInnerPoints.size());
        for (std::size_t lsI = 0; lsI < innerPointCount; ++lsI)
        {
            const double deltaUm = codeLengthToUm((iGeometry.levelSets[lsI].innerPoint() - refInnerPoints[lsI]).mod());
            maxInnerDelta_um = std::max(maxInnerDelta_um, deltaUm);
            sumInnerDelta2_um2 += deltaUm * deltaUm;
            ++innerCount;
        }

        const double rmsPointDelta_um = (pointCount > 0)
                                            ? std::sqrt(sumPointDelta2_um2 / static_cast<double>(pointCount))
                                            : 0.0;
        const double rmsArcDelta_um = (arcCount > 0)
                                          ? std::sqrt(sumArcDelta2_um2 / static_cast<double>(arcCount))
                                          : 0.0;
        const double rmsInnerDelta_um = (innerCount > 0)
                                            ? std::sqrt(sumInnerDelta2_um2 / static_cast<double>(innerCount))
                                            : 0.0;

        const GeometryMetrics liveMetrics = computeGeometryMetricsFromLines_(iGeometry.lines);
        const GeometryMetrics refMetrics = computeGeometryMetricsFromLines_(refLines);

        consistencyOut_ << step << ',' << time_s << ',' << (updated ? 1 : 0) << ','
                        << (linesCurrentGeometryFrame_ ? 1 : 0) << ','
                        << (lastUpdateUsedAdvectedPath_ ? 1 : 0) << ','
                        << pointCount << ',' << innerCount << ','
                        << maxPointDelta_um << ',' << rmsPointDelta_um << ','
                        << maxArcDelta_um << ',' << rmsArcDelta_um << ','
                        << maxInnerDelta_um << ',' << rmsInnerDelta_um << ','
                        << liveMetrics.gsa_mm << ',' << refMetrics.gsa_mm << ','
                        << liveMetrics.gad_mm << ',' << refMetrics.gad_mm << ','
                        << liveMetrics.da_mm << ',' << refMetrics.da_mm
                        << "\n";
        consistencyOut_.flush();
    }

    void writeRegionMetrics_(int step, double time_s, bool updated)
    {
        if (!regionOut_)
            return;
        for (std::size_t ri = 0; ri < config_.regions.size(); ++ri)
        {
            const auto& region = config_.regions[ri];
            const auto& state = states_[ri];
            regionOut_ << step << ',' << time_s << ',' << (updated ? 1 : 0) << ',' << region.name << ','
                       << region.boundaryId << ',' << (region.levelSetIndex + 1) << ',' << state.lastSampleCount << ','
                       << state.lastResolvedIonCount << ',' << state.lastDirectMatchCount << ','
                       << state.lastFallbackMatchCount << ','
                       << state.lastArea_m2 << ',' << state.smoothedFlux << ',' << state.smoothedEnergy << ','
                       << state.smoothedCos << ',' << state.lastNormalSpeed_mps << ','
                       << state.normalDirectionSign << ',' << (state.normalDirectionInitialized ? 1 : 0) << ','
                       << (state.activeThisStep ? 1 : 0) << "\n";
        }
        regionOut_.flush();
    }

    void writeImpactStats_(int step, double time_s, bool updated)
    {
        if (!impactOut_)
            return;
        impactOut_ << step << ',' << time_s << ',' << (updated ? 1 : 0) << ','
                   << resolvedImpactsStep_ << ',' << unresolvedImpactsStep_ << ','
                   << ignoredNonIonImpactsStep_ << ',' << regionMissImpactsStep_ << ','
                   << updateWindowSamples_.size() << "\n";
        impactOut_.flush();
    }

    void writeMatchStats_(int step, double time_s, bool updated)
    {
        if (!matchOut_)
            return;
        for (const auto& row : lastMatchStatsRows_)
        {
            matchOut_ << step << ',' << time_s << ',' << (updated ? 1 : 0) << ','
                      << row.boundaryId << ',' << (row.levelSetIndex + 1) << ','
                      << row.resolvedIonCount << ',' << row.matchedResolvedIonCount << ','
                      << row.unmatchedResolvedIonCount << ',' << row.directMatchCount << ','
                      << row.fallbackMatchCount << "\n";
        }
        matchOut_.flush();
    }

    static std::string toLower_(std::string s)
    {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return s;
    }

    static void writeString_(std::ofstream& out, const std::string& value)
    {
        const int n = static_cast<int>(value.size());
        out.write(reinterpret_cast<const char*>(&n), sizeof(int));
        out.write(value.data(), n);
    }

    static std::string readString_(std::ifstream& in)
    {
        int n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(int));
        if (n <= 0)
            return {};
        std::string value(static_cast<std::size_t>(n), '\0');
        in.read(&value[0], n);
        return value;
    }

    static void writeField_(std::ofstream& out, const VertField<double>& field)
    {
        const int n = field.size();
        out.write(reinterpret_cast<const char*>(&n), sizeof(int));
        for (int i = 0; i < n; ++i)
        {
            const int label = field.label(i);
            out.write(reinterpret_cast<const char*>(&label), sizeof(int));
            out.write(reinterpret_cast<const char*>(&field.value[i]), sizeof(double));
        }
    }

    template<typename Mesh>
    static void readField_(std::ifstream& in, const Mesh& mesh, VertField<double>& field)
    {
        int n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(int));
        field.initialize(mesh, 0.0);
        for (int i = 0; i < n; ++i)
        {
            int label = 0;
            double value = 0.0;
            in.read(reinterpret_cast<char*>(&label), sizeof(int));
            in.read(reinterpret_cast<char*>(&value), sizeof(double));
            if (i < field.size())
            {
                field.setLabel(i, label);
                field.value[i] = value;
            }
        }
    }
};
