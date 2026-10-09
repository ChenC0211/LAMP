#include "Modules/Plasma/Physics/PicMeshRefinement.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/ExplicitPIC/Deposit/ParticleStencilQuery.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/Physics/PicChargeDeposition.h"
#include "Modules/Plasma/Physics/PicPoissonFieldSolve.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include "Modules/Plasma/Boundary/PicMovingGeometry.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicPolarization.h"
#include "Modules/Plasma/ImplicitPIC/Orbit/ImplicitPicOrbit.h"
#include "Core/Mesh/QuadtreeMeshRefiner.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <stdexcept>

namespace {
bool isLeafCell(const Cell& c) { return c.son[0] == -1; }
std::vector<int> countLeafLevels(const QuadtreeMesh& mesh)
{
    std::vector<int> counts(static_cast<std::size_t>(std::max(0, mesh.maxLevel()) + 1), 0);
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const Cell& c = mesh.cells(mesh.cellIndexToLabel[ci]);
        if (!isLeafCell(c)) continue;
        const int level = c.level();
        if (level < 0) continue;
        if (level >= static_cast<int>(counts.size())) counts.resize(static_cast<std::size_t>(level + 1), 0);
        counts[static_cast<std::size_t>(level)] += 1;
    }
    return counts;
}
void adapt(QuadtreeMesh& mesh, CellField<double>& indicator, CellField<double>& target,
    std::vector<VertField<double>*>& vertices, std::vector<CellField<double>*> cells)
{
    std::vector<VertField<double>> levelSets;
    std::vector<HEdgeField<double>*> horizontal;
    std::vector<VEdgeField<double>*> vertical;
    mesh.adapt(levelSets, indicator, vertices, horizontal, vertical, cells, &target);
}
double relL2Change(const std::vector<double>& now, const std::vector<double>& prev) {
    if (now.size() != prev.size() || now.empty())
        return std::numeric_limits<double>::infinity();
    long double num = 0.0L;
    long double den = 0.0L;
    for (std::size_t i = 0; i < now.size(); ++i)
    {
        const double a = std::isfinite(now[i]) ? now[i] : 0.0;
        const double b = std::isfinite(prev[i]) ? prev[i] : 0.0;
        const long double d = static_cast<long double>(a - b);
        num += d * d;
        den += static_cast<long double>(b) * static_cast<long double>(b);
    }
    return std::sqrt(static_cast<double>(num / std::max<long double>(1e-300L, den)));
}

}

ExplicitPicRefinement::ExplicitPicRefinement(QuadtreeMesh& mesh)
    : mesh(mesh), amrNeAvg(mesh, 0.0), amrEAvg(mesh, 0.0), amrIndicator(mesh, 0.0), amrTargetLevel(mesh, 0.0) {}

void ExplicitPicRefinement::adaptMesh(std::vector<VertField<double>*>& vertices)
{
    adapt(mesh, amrIndicator, amrTargetLevel, vertices, {&amrNeAvg, &amrEAvg, &amrTargetLevel});
}

bool ExplicitPicRefinement::adaptIfNeeded(const AmrConfig& config, int step, const Adaptation& state)
{
    if (!state.enabled || !pending() || config.adaptEverySteps <= 0
        || step % config.adaptEverySteps != 0) return false;
    if (&state.domain.mesh() != &mesh || &state.fields.mesh_ != &mesh
        || &state.deposition.fields_ != &state.fields || &state.poisson.fields_ != &state.fields
        || &state.coupling.mesh_ != &mesh
        || state.domain.mode_ != PicDomain::Mode::Explicit
        || state.geometry.mode_ != PicAmrGeometry::Mode::Explicit
        || state.moving.mode_ != PicMovingGeometry::Mode::Explicit
        || state.coupling.mode_ != PicCurrentWeights::Policy::Mode::Explicit)
        throw std::invalid_argument("Explicit PIC AMR requires one mesh and consistent explicit owners");
    if (!state.particles.initialized_ || amrTargetLevel.size() != mesh.cellIndexNum())
        throw std::logic_error("Explicit PIC AMR requires initialized particles and current refinement targets");

    AdaptationReport report;
    report.levelsBefore = countLeafLevels(mesh);
    const auto start = std::chrono::steady_clock::now();
    auto* geometry = state.domain.geometry();
    auto& embedded = state.domain.embedded();
    std::vector<VertField<double>*> vertices;
    state.geometry.captureBeforeAdapt(geometry);
    if (geometry)
        for (auto& levelSet : geometry->levelSets) vertices.push_back(&levelSet);
    adaptMesh(vertices);
    if (state.moving.controller()) state.moving.onMeshChanged(mesh);
    state.geometry.refresh(mesh, geometry, embedded, state.moving.controller(),
        {step, 0.0, state.embeddedEnabled, state.geometryRefinementActive, state.domain.refiner().phaseActive(), false});
    state.domain.rebuildFaceGraph();
    if (embedded.ready())
        embedded.rebuildMetrics(mesh, geometry, state.domain.faceGraph(), state.lineToBoundaryId, state.moving.controller());
    state.fields.resetAfterMeshChange();
    state.deposition.resetAfterMeshChange();
    state.particles.relabelAfterMeshChange(mesh, PicSpeciesState::RelabelOptions::fromExplicitConfig(
        config.particleRelabel.mode, config.particleRelabel.fallback));
    state.poisson.resetAfterMeshChange();
    state.coupling.refreshAfterMeshChange();
    adaptationCompleted();
    report.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (state.reportLevels) report.levelsAfter = countLeafLevels(mesh);
    adaptationReport_ = std::move(report);
    return true;
}

ImplicitPicRefinement::ImplicitPicRefinement(QuadtreeMesh& mesh)
    : mesh(mesh), amrTargetLevel(mesh, 0.0), amrIndicator(mesh, 0.0) {}

void ImplicitPicRefinement::adaptMesh(std::vector<VertField<double>*>& vertices, CellField<double>* potential)
{
    std::vector<CellField<double>*> cells{&amrTargetLevel};
    if (potential) cells.push_back(potential);
    adapt(mesh, amrIndicator, amrTargetLevel, vertices, std::move(cells));
}

bool ImplicitPicRefinement::adaptIfNeeded(int step, const Adaptation& state)
{
    if (!state.enabled || !pending() || state.every <= 0 || step % std::max(1, state.every) != 0) return false;
    const auto& input = state.fieldInput;
    if (&state.domain.mesh() != &mesh || &state.fields.mesh_ != &mesh
        || &state.history.fields_ != &state.fields || &state.poisson.fields_ != &state.fields
        || &state.polarization.mesh != &mesh || &state.coupling.mesh_ != &mesh
        || &state.currentWeights.mesh_ != &mesh || &state.orbit.mesh != &mesh
        || &input.graph != &state.domain.faceGraph()
        || (input.cut && input.cut != &state.domain.embedded().metrics())
        || &state.rebuild.species != &state.particles.species()
        || &state.rebuild.macroWeight != &state.particles.weights()
        || &state.rebuild.coupling != &state.coupling || &state.rebuild.embedded != &state.domain.embedded()
        || state.domain.mode_ != PicDomain::Mode::Implicit
        || state.geometry.mode_ != PicAmrGeometry::Mode::Implicit
        || state.moving.mode_ != PicMovingGeometry::Mode::Implicit
        || state.coupling.mode_ != PicCurrentWeights::Policy::Mode::Implicit)
        throw std::invalid_argument("Implicit PIC AMR requires one mesh and consistent implicit owners");
    if (!state.particles.initialized_ || amrTargetLevel.size() != mesh.cellIndexNum()
        || state.polarization.Nstar.size() != state.particles.species().size()
        || state.orbit.cachesBar.size() != state.particles.species().size()
        || !state.refreshOutputCoordinates)
        throw std::logic_error("Implicit PIC AMR requires initialized species, current targets and an output refresh hook");

    using Clock = std::chrono::steady_clock;
    const auto secondsSince = [](Clock::time_point begin) {
        return std::chrono::duration<double>(Clock::now() - begin).count();
    };
    const auto begin = Clock::now();
    AdaptationReport report;
    report.cellsBefore = mesh.cellIndexNum();
    const auto& phi = state.fields.potential();
    report.useFastRemap = state.fastRemap && phi.size() == report.cellsBefore && phi.value != nullptr;
    CellField<double> remapped;
    if (report.useFastRemap) remapped = phi;
    auto* geometry = state.domain.geometry();
    auto& embedded = state.domain.embedded();
    state.geometry.captureBeforeAdapt(geometry);
    std::vector<VertField<double>*> vertices;
    if (geometry) for (auto& levelSet : geometry->levelSets) vertices.push_back(&levelSet);
    const auto meshStart = Clock::now();
    adaptMesh(vertices, report.useFastRemap ? &remapped : nullptr);
    report.meshAdaptS = secondsSince(meshStart);
    const double change = static_cast<double>(std::abs(mesh.cellIndexNum() - report.cellsBefore))
        / static_cast<double>(std::max(1, report.cellsBefore));
    if (state.moving.controller())
    {
        const auto movingStart = Clock::now();
        state.moving.onMeshChanged(mesh);
        report.movingEbMeshChangedS = secondsSince(movingStart);
    }
    state.geometry.refresh(mesh, geometry, embedded, state.moving.controller(),
        {step, change, state.embeddedEnabled, state.geometryRefinementActive,
         state.domain.refiner().phaseActive(), state.moving.updated()});
    if (geometry) geometry->dumpAllLevelSetsScatter(mesh, "s6_post_amr");
    const auto graphStart = Clock::now();
    state.domain.rebuildFaceGraph();
    report.faceGraphS = secondsSince(graphStart);
    if (state.embeddedEnabled && embedded.ready())
    {
        const auto metricsStart = Clock::now();
        embedded.rebuildMetrics(mesh, geometry, state.domain.faceGraph(), state.lineToBoundaryId, state.moving.controller());
        report.embeddedMetricsS = secondsSince(metricsStart);
    }
    const auto resetStart = Clock::now();
    state.fields.resetAfterMeshChange(report.useFastRemap ? &remapped : nullptr);
    state.history.resetAfterMeshChange();
    state.polarization.resetDepositedFields();
    resetFieldsAfterMeshChange();
    state.polarization.resetCoefficients(state.domain.faceGraph().segments().size());
    report.fieldResetS = secondsSince(resetStart);
    state.coupling.refreshAfterMeshChange();
    const auto relabelStart = Clock::now();
    state.particles.relabelAfterMeshChange(mesh, {true, true});
    report.relabelParticlesS = secondsSince(relabelStart);
    const auto rebuildStart = Clock::now();
    if (report.useFastRemap) state.poisson.reseedImplicitFromPotential(input, state.history);
    else
    {
        PicCurrentWeights rebuiltWeights(mesh);
        state.poisson.rebuildImplicit(input, state.rebuild, rebuiltWeights, &state.history);
    }
    report.fieldRebuildS = secondsSince(rebuildStart);
    const auto coordinateStart = Clock::now();
    state.refreshOutputCoordinates(mesh);
    report.coordRefreshS = secondsSince(coordinateStart);
    adaptationCompleted();
    if (state.resetTransients)
    {
        state.currentWeights.clear();
        state.orbit.resetScratch();
        state.polarization.resetScratch();
    }
    report.cellsAfter = mesh.cellIndexNum();
    report.seconds = secondsSince(begin);
    adaptationReport_ = std::move(report);
    return true;
}

void ImplicitPicRefinement::resetFieldsAfterMeshChange()
{
    amrTargetLevel.initialize(mesh, 0.0);
    amrIndicator.initialize(mesh, 0.0);
}

void ImplicitPicRefinement::requestAdaptation()
{
    amrPending = true;
    amrPendingCount = std::max(amrPendingCount, 1);
}

bool ExplicitPicRefinement::update(const AmrConfig& amrCfg, const Sample& sample, double amrTe_eV,
    int stepForPrint, int steps, QuadtreeMeshRefiner& meshRefiner)
{
    const bool useDebye = amrCfg.debye.enabled;
    const bool useEField = amrCfg.indicator.useEField;
    const int nCells = mesh.cellIndexNum();
    if (nCells <= 0)
        return false;

    const auto& norm = sample.normalization;
    const auto& units = sample.units;
    const auto& amrNeInstant = sample.numberAmount;
    const auto& Ecell = sample.electricField;

    const int K_ne = std::max(1, amrCfg.debye.avgWindowSteps);
    const double alphaNe = useDebye ? (2.0 / (static_cast<double>(K_ne) + 1.0)) : 0.0;
    const int K_E = std::max(1, amrCfg.indicator.eFieldAvgWindow);
    const double alphaE = useEField ? (2.0 / (static_cast<double>(K_E) + 1.0)) : 0.0;
    const int K_Eref = std::max(1, amrCfg.indicator.eFieldRefAvgWindow);
    const double alphaEref = useEField ? (2.0 / (static_cast<double>(K_Eref) + 1.0)) : 0.0;

    const double Te_eV = std::max(1e-6, amrTe_eV);
    const double Te_J = Te_eV * PhysicalConstants::eV_J;

    const double cd = std::max(1e-6, amrCfg.debye.CD);
    const double neFloor = std::max(1e-300, amrCfg.debye.neFloor);
    const double eAbsFloor = std::max(0.0, amrCfg.indicator.eFieldAbsFloor);
    double eRefFrac = amrCfg.indicator.eFieldRefFrac;
    if (!std::isfinite(eRefFrac))
        eRefFrac = 0.0;
    if (eRefFrac < 0.0)
        eRefFrac = 0.0;
    if (eRefFrac > 1.0)
        eRefFrac = 1.0;
    const double eps0 = PhysicalConstants::epsilon0_Fpm;
    const double qe = PhysicalConstants::qe_C;

    std::vector<double> ratio(nCells, 1.0);
    std::vector<double> Eavg(nCells, 0.0);
    std::vector<double> debyeLengthCode;
    std::vector<double> eAvgCode;
    std::vector<double> lambdaD_si;
    std::vector<double> h_si_cell;
    if (useDebye)
    {
        debyeLengthCode.assign(nCells, 0.0);
        lambdaD_si.assign(nCells, 0.0);
        h_si_cell.assign(nCells, 0.0);
    }
    if (useEField)
        eAvgCode.assign(nCells, 0.0);
    std::vector<double> eSamples;
    if (useEField)
        eSamples.reserve(nCells / 4);

    const bool hasAlphaAmr = sample.alpha &&
        (sample.alpha->size() == static_cast<size_t>(mesh.cellIndexNum()));
    auto alphaForCellAmr = [&](int ci) -> double {
        if (!hasAlphaAmr)
            return 1.0;
        double alpha = (*sample.alpha)[static_cast<size_t>(ci)];
        if (!std::isfinite(alpha))
            alpha = 1.0;
        if (alpha < 0.0)
            alpha = 0.0;
        else if (alpha > 1.0)
            alpha = 1.0;
        return alpha;
    };

    const int indicatorEvery = std::max(1, amrCfg.indicatorEverySteps);
    const bool doHeavy = (indicatorEvery <= 1) ||
        ((stepForPrint % indicatorEvery) == 0) || (stepForPrint == steps);

    double ratio_sum = 0.0;
    double ratio_max = 0.0;
    double lambdaD_min = std::numeric_limits<double>::max();
    double h_min = std::numeric_limits<double>::max();
    int leafCount = 0;
    double Emax = 0.0;
    int debyeLowNeCount = 0;
    int debyeBadCount = 0;
    [[maybe_unused]] const bool amrUseOmp = (nCells > 2048);

    // First pass: update EMA fields and compute Debye ratios.
#ifdef _OPENMP
#pragma omp parallel for if(amrUseOmp) reduction(+:ratio_sum, leafCount, debyeLowNeCount, debyeBadCount) reduction(max:ratio_max, Emax) reduction(min:lambdaD_min, h_min)
#endif
    for (int ci = 0; ci < nCells; ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;

        const double area_code = std::max(1e-300, c.area());
        const double alpha = alphaForCellAmr(ci);
        if (!(alpha > 0.0))
        {
            ratio[ci] = 0.0;
            if (useDebye)
                amrNeAvg.value[ci] = 0.0;
            if (useEField)
                amrEAvg.value[ci] = 0.0;
            Eavg[ci] = 0.0;
            if (useDebye)
            {
                debyeLengthCode[ci] = 0.0;
                lambdaD_si[ci] = 0.0;
                h_si_cell[ci] = 0.0;
            }
            if (useEField)
                eAvgCode[ci] = 0.0;
            continue;
        }
        const Vector2D sz = c.size();
        const double h_code = std::min(sz.x, sz.y);
        const double h_si = useDebye ? units.LengthCodeToSI(h_code) : 0.0;
        if (useDebye)
            h_si_cell[ci] = h_si;
        double lambdaD = 0.0;
        if (useDebye)
        {
            const double vol_SI = std::max(1e-300, norm.VolumeCodeToSI(area_code * alpha));
            const double ne = amrNeInstant.value[ci] / vol_SI;

            double ne_avg = (amrNeAvgInit ? ((1.0 - alphaNe) * amrNeAvg.value[ci] + alphaNe * ne) : ne);
            if (!std::isfinite(ne_avg) || ne_avg < 0.0)
                ne_avg = 0.0;
            amrNeAvg.value[ci] = ne_avg;

            const bool lowNe = (ne_avg <= neFloor);
            const double ne_eff = lowNe ? neFloor : ne_avg;
            if (lowNe)
                debyeLowNeCount += 1;
            lambdaD = std::sqrt(std::max(0.0, eps0 * Te_J / (ne_eff * qe * qe)));
            if (std::isfinite(lambdaD) && lambdaD > 0.0)
            {
                lambdaD_si[ci] = lambdaD;
            }
            else
            {
                lambdaD_si[ci] = 0.0;
                debyeBadCount += 1;
            }
        }
        else
        {
            amrNeAvg.value[ci] = 0.0;
        }

        double r = 1.0;
        if (useDebye)
        {
            r = (lambdaD > 0.0 && h_si > 0.0) ? (h_si / (cd * lambdaD)) : 1.0;
            if (!std::isfinite(r) || r <= 0.0)
                r = 1.0;
        }
        ratio[ci] = r;

        if (useEField)
        {
            const double Emag_code = Ecell.value[ci].mod();
            const double Emag_si = units.ElectricFieldCodeToSI(Emag_code);
            const double e_avg = (amrEAvgInit ? ((1.0 - alphaE) * amrEAvg.value[ci] + alphaE * Emag_si) : Emag_si);
            amrEAvg.value[ci] = e_avg;
            Eavg[ci] = e_avg;
        }
        else
        {
            Eavg[ci] = 0.0;
            amrEAvg.value[ci] = 0.0;
        }
    }

    if (useDebye)
        amrNeAvgInit = true;
    if (useEField)
        amrEAvgInit = true;

    if (!doHeavy)
        return false;

    // Optional spatial smoothing (leaf-cell neighborhood).
    auto smoothSpatial = [&](const std::vector<double>& in,
                             std::vector<double>& out,
                             int radiusCells,
                             bool ignoreNonPositive) -> void
    {
        if (radiusCells <= 0)
        {
            out = in;
            return;
        }
        out.assign(nCells, 0.0);
        std::vector<int> stencil;
        stencil.reserve(64);
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isLeafCell(c))
            {
                out[ci] = in[ci];
                continue;
            }
            if (!(alphaForCellAmr(ci) > 0.0))
            {
                out[ci] = 0.0;
                continue;
            }
            const Vector2D cc = c.centerCoord();
            const double h = std::min(c.size().x, c.size().y);
            const double hp = static_cast<double>(radiusCells) * h;
            if (!(hp > 0.0))
            {
                out[ci] = in[ci];
                continue;
            }
            ParticleStencilQuery::enumerateLeafCellsIntersectParticleSupport(
                mesh, cc.x, cc.y, hp, &stencil, /*debugPrint=*/false, c.label(), /*deterministic=*/true);
            double sum = 0.0;
            double wsum = 0.0;
            for (int lbl : stencil)
            {
                if (lbl < 0)
                    continue;
                const Cell& nb = mesh.cells(lbl);
                if (!isLeafCell(nb))
                    continue;
                const int ni = nb.index();
                if (ni < 0 || ni >= nCells)
                    continue;
                const double alpha = alphaForCellAmr(ni);
                if (!(alpha > 0.0))
                    continue;
                const double val = in[ni];
                if (!std::isfinite(val))
                    continue;
                if (ignoreNonPositive && val <= 0.0)
                    continue;
                const double w = nb.area() * alpha;
                sum += w * val;
                wsum += w;
            }
            out[ci] = (wsum > 0.0) ? (sum / wsum) : in[ci];
        }
    };

    std::vector<double> lambdaD_use;
    std::vector<double> Eavg_use;
    if (useDebye)
    {
        const int smoothCells = std::max(0, amrCfg.debye.smoothWindowCells);
        smoothSpatial(lambdaD_si, lambdaD_use, smoothCells, /*ignoreNonPositive=*/true);
    }
    if (useEField)
    {
        const int smoothCells = std::max(0, amrCfg.indicator.eFieldSmoothWindowCells);
        smoothSpatial(Eavg, Eavg_use, smoothCells, /*ignoreNonPositive=*/false);
    }
    if (!useDebye)
        lambdaD_use.clear();
    if (!useEField)
        Eavg_use.clear();

    // Recompute Debye ratio stats and build smoothed indicator fields.
    ratio_sum = 0.0;
    ratio_max = 0.0;
    lambdaD_min = std::numeric_limits<double>::max();
    h_min = std::numeric_limits<double>::max();
    leafCount = 0;
    Emax = 0.0;

    for (int ci = 0; ci < nCells; ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        const double alpha = alphaForCellAmr(ci);
        if (!(alpha > 0.0))
        {
            ratio[ci] = 0.0;
            if (useDebye)
                debyeLengthCode[ci] = 0.0;
            if (useEField)
                eAvgCode[ci] = 0.0;
            continue;
        }

        if (useDebye)
        {
            const double lambdaD = (ci < static_cast<int>(lambdaD_use.size())) ? lambdaD_use[ci] : 0.0;
            if (std::isfinite(lambdaD) && lambdaD > 0.0)
            {
                debyeLengthCode[ci] = units.LengthSIToCode(lambdaD);
                const double h_si = h_si_cell[ci];
                double r = (h_si > 0.0) ? (h_si / (cd * lambdaD)) : 1.0;
                if (!std::isfinite(r) || r <= 0.0)
                    r = 1.0;
                ratio[ci] = r;
                ratio_sum += r;
                ratio_max = std::max(ratio_max, r);
                lambdaD_min = std::min(lambdaD_min, lambdaD);
                if (h_si > 0.0)
                    h_min = std::min(h_min, h_si);
            }
            else
            {
                debyeLengthCode[ci] = 0.0;
                ratio[ci] = 1.0;
            }
        }
        else
        {
            ratio[ci] = 1.0;
        }

        if (useEField)
        {
            const double e_avg = (ci < static_cast<int>(Eavg_use.size())) ? Eavg_use[ci] : 0.0;
            eAvgCode[ci] = units.ElectricFieldSIToCode(e_avg);
            if (std::isfinite(e_avg) && e_avg > 0.0)
                Emax = std::max(Emax, e_avg);
        }
        leafCount += 1;
    }

    double Eref_step = Emax;
    if (useEField)
    {
        eSamples.clear();
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (!isLeafCell(c))
                continue;
            if (!(alphaForCellAmr(ci) > 0.0))
                continue;
            const double val = (ci < static_cast<int>(Eavg_use.size())) ? Eavg_use[ci] : 0.0;
            if (std::isfinite(val) && val > 0.0)
                eSamples.push_back(val);
        }
        if (!eSamples.empty())
        {
            double pct = amrCfg.indicator.eFieldPercentile;
            if (pct < 0.0) pct = 0.0;
            if (pct > 100.0) pct = 100.0;
            const size_t n = eSamples.size();
            const size_t k = static_cast<size_t>(std::ceil((pct / 100.0) * static_cast<double>(n - 1)));
            if (k < n)
            {
                std::nth_element(eSamples.begin(), eSamples.begin() + k, eSamples.end());
                Eref_step = eSamples[k];
            }
        }
        if (!(Eref_step > 0.0))
            Eref_step = Emax;
        if (eRefFrac > 0.0 && Emax > 0.0)
            Eref_step = std::max(Eref_step, eRefFrac * Emax);
        if (eAbsFloor > 0.0)
            Eref_step = std::max(Eref_step, eAbsFloor);

        if (amrErefAvgInit)
            amrErefAvg = (1.0 - alphaEref) * amrErefAvg + alphaEref * Eref_step;
        else
        {
            amrErefAvg = Eref_step;
            amrErefAvgInit = true;
        }
    }

    const double Eref_use = useEField ? amrErefAvg : 0.0;
    const double Eref_code = useEField ? units.ElectricFieldSIToCode(Eref_use) : 0.0;
    if (useDebye)
        refinementIndicators.setDebye(debyeLengthCode, cd);
    else
        refinementIndicators.clearDebye();
    if (useEField)
        refinementIndicators.setEField(eAvgCode, Eref_code);
    else
        refinementIndicators.clearEField();

    QuadtreeMeshRefiner::BudgetOptions budgetOpt;
    QuadtreeMeshRefiner::BudgetOptions* budgetPtr = nullptr;
    if (amrCfg.cellBudgetPerAdapt > 0)
    {
        budgetOpt.cellBudget = std::max(0, amrCfg.cellBudgetPerAdapt);
        budgetPtr = &budgetOpt;
    }

    meshRefiner.buildTargetLevels(mesh, amrTargetLevel, &amrIndicator, budgetPtr,
                                  false, nullptr, &refinementIndicators.fields());

    int refineCount = 0;
    int coarsenCount = 0;
#ifdef _OPENMP
#pragma omp parallel for if(amrUseOmp) reduction(+:refineCount, coarsenCount)
#endif
    for (int ci = 0; ci < nCells; ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (!isLeafCell(c))
            continue;
        const int aimLevel = static_cast<int>(std::lround(amrTargetLevel.value[ci]));
        if (aimLevel > c.level())
            refineCount++;
        else if (aimLevel < c.level())
            coarsenCount++;
    }

    const double ratio_mean = (leafCount > 0) ? (ratio_sum / static_cast<double>(leafCount)) : 0.0;
    if (useDebye)
    {
        amrStats.lambdaD_min_m = std::isfinite(lambdaD_min) ? lambdaD_min : 0.0;
        amrStats.h_min_m = std::isfinite(h_min) ? h_min : 0.0;
        amrStats.ratio_max = ratio_max;
        amrStats.ratio_mean = ratio_mean;
    }
    else
    {
        amrStats.lambdaD_min_m = 0.0;
        amrStats.h_min_m = 0.0;
        amrStats.ratio_max = 0.0;
        amrStats.ratio_mean = 0.0;
    }
    amrStats.refine_cells = refineCount;
    amrStats.coarsen_cells = coarsenCount;

    const bool needChange = (refineCount > 0) || (coarsenCount > 0);
    if (amrCfg.persistenceSteps <= 1)
    {
        amrPending = needChange;
        amrPendingCount = needChange ? 1 : 0;
    }
    else
    {
        if (needChange)
            amrPendingCount += 1;
        else
            amrPendingCount = 0;
        amrPending = (amrPendingCount >= amrCfg.persistenceSteps);
    }

    amrStats.debyeLowNeCount = debyeLowNeCount;
    amrStats.debyeBadCount = debyeBadCount;
    amrStats.neFloor = neFloor;
    return true;
}

void ImplicitPicRefinement::update(const Options& options, const Sample& sample, int stepForPrint,
    QuadtreeMeshRefiner& meshRefiner)
{
    auto alphaForCellAmr = [&](int ci) -> double {
        if (!sample.alpha)
            return 1.0;
        if (ci < 0 || ci >= static_cast<int>(sample.alpha->size()))
            return 1.0;
        double a = (*sample.alpha)[static_cast<std::size_t>(ci)];
        if (!std::isfinite(a))
            a = 1.0;
        if (a < 0.0) a = 0.0;
        if (a > 1.0) a = 1.0;
        return a;
    };

    amrSkippedByGateLast = false;
    if (!options.enabled)
        return;
    if (options.indicatorEvery <= 0)
        return;
    if ((stepForPrint % std::max(1, options.indicatorEvery)) != 0)
        return;

    const int nCells = mesh.cellIndexNum();
    if (nCells <= 0)
        return;

    std::vector<double> eMag;
    std::vector<double> lambdaD_code;
    double relEChange = std::numeric_limits<double>::infinity();
    double relDebyeChange = std::numeric_limits<double>::infinity();

    if (options.wantEField)
    {
        eMag.assign(static_cast<std::size_t>(nCells), 0.0);
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (c.son[0] != -1)
                continue;
            if (!(alphaForCellAmr(ci) > 0.0))
                continue;
            const double v = sample.electricField.value[ci].mod();
            eMag[static_cast<std::size_t>(ci)] = (std::isfinite(v) && v > 0.0) ? v : 0.0;
        }
        relEChange = relL2Change(eMag, prevAmrEFieldIndicator);
    }
    else
    {
        refinementIndicators.clearEField();
        prevAmrEFieldIndicator.clear();
    }

    if (options.wantDebye && options.electronIndex >= 0 && options.electronIndex < static_cast<int>(sample.numberAmounts.size()))
    {
        lambdaD_code.assign(static_cast<std::size_t>(nCells), 0.0);
        const double Te_J = std::max(1e-6, options.temperatureEV) * PhysicalConstants::eV_J;
        const double eps0 = PhysicalConstants::epsilon0_Fpm;
        const double qe = PhysicalConstants::qe_C;
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (c.son[0] != -1)
                continue;
            const double alpha = alphaForCellAmr(ci);
            if (!(alpha > 0.0))
                continue;
            const double area_code = std::max(1e-300, c.area());
            const double vol_SI = std::max(1e-300, sample.normalization.VolumeCodeToSI(area_code * alpha));
            const double ne = sample.numberAmounts[static_cast<std::size_t>(options.electronIndex)].value[ci] / vol_SI;
            if (!std::isfinite(ne) || ne <= 0.0)
                continue;
            const double lambdaD_si = std::sqrt(std::max(0.0, eps0 * Te_J / (ne * qe * qe)));
            if (!std::isfinite(lambdaD_si) || !(lambdaD_si > 0.0))
                continue;
            lambdaD_code[static_cast<std::size_t>(ci)] = sample.units->LengthSIToCode(lambdaD_si);
        }
        relDebyeChange = relL2Change(lambdaD_code, prevAmrDebyeIndicator);
    }
    else
    {
        refinementIndicators.clearDebye();
        prevAmrDebyeIndicator.clear();
    }

    const bool gateE = options.wantEField && (options.gateElectricRelative > 0.0) && std::isfinite(relEChange) && (relEChange < options.gateElectricRelative);
    const bool gateDebye = options.wantDebye && (options.gateDebyeRelative > 0.0) && std::isfinite(relDebyeChange) && (relDebyeChange < options.gateDebyeRelative);
    const bool shouldGateAmr = (options.wantEField || options.wantDebye) &&
                               ((!options.wantEField || gateE) && (!options.wantDebye || gateDebye));
    if (shouldGateAmr)
    {
        amrSkippedByGateLast = true;
        amrPending = false;
        amrPendingCount = 0;
        return;
    }

    if (options.wantEField)
    {
        refinementIndicators.setEField(std::move(eMag), /*reference=*/-1.0, /*percentile=*/options.electricPercentile);
        prevAmrEFieldIndicator = refinementIndicators.electricFieldValues();
    }
    if (options.wantDebye)
    {
        refinementIndicators.setDebye(std::move(lambdaD_code), (options.debyeScale > 0.0) ? options.debyeScale : 1.0);
        prevAmrDebyeIndicator = refinementIndicators.debyeValues();
    }

    meshRefiner.buildTargetLevels(mesh, amrTargetLevel, &amrIndicator, nullptr,
                                  false, nullptr, &refinementIndicators.fields());

    struct AmrCandidate
    {
        int ci = -1;
        int cur = 0;
        int aim = 0;
        double score = 0.0;
    };
    std::vector<AmrCandidate> candidates;
    candidates.reserve(static_cast<std::size_t>(nCells));
    int refineCount = 0;
    int coarsenCount = 0;
    for (int ci = 0; ci < nCells; ++ci)
    {
        const int cLbl = mesh.cellIndexToLabel[ci];
        const Cell& c = mesh.cells(cLbl);
        if (c.son[0] != -1)
            continue;
        const int cur = c.level();
        const int aim = static_cast<int>(std::lround(amrTargetLevel.value[ci]));
        if (aim > cur)
            refineCount++;
        else if (aim < cur)
            coarsenCount++;
        if (aim != cur)
        {
            double score = 0.0;
            if (ci >= 0 && ci < mesh.cellIndexNum() && std::isfinite(amrIndicator.value[ci]))
                score = std::abs(amrIndicator.value[ci]);
            candidates.push_back(AmrCandidate{ci, cur, aim, score});
        }
    }

    const int totalChanges = refineCount + coarsenCount;
    if (options.maxChanges > 0 && totalChanges > options.maxChanges)
    {
        std::sort(candidates.begin(), candidates.end(), [](const AmrCandidate& a, const AmrCandidate& b) {
            if (a.score == b.score)
                return a.ci < b.ci;
            return a.score > b.score;
        });

        std::vector<unsigned char> keep(static_cast<std::size_t>(nCells), 0u);
        const int keepN = std::min(options.maxChanges, static_cast<int>(candidates.size()));
        for (int i = 0; i < keepN; ++i)
        {
            const int ci = candidates[static_cast<std::size_t>(i)].ci;
            if (ci >= 0 && ci < nCells)
                keep[static_cast<std::size_t>(ci)] = 1u;
        }
        for (const auto& cand : candidates)
        {
            if (cand.ci < 0 || cand.ci >= nCells)
                continue;
            if (keep[static_cast<std::size_t>(cand.ci)] == 0u)
                amrTargetLevel.value[cand.ci] = static_cast<double>(cand.cur);
        }

        refineCount = 0;
        coarsenCount = 0;
        for (int ci = 0; ci < nCells; ++ci)
        {
            const int cLbl = mesh.cellIndexToLabel[ci];
            const Cell& c = mesh.cells(cLbl);
            if (c.son[0] != -1)
                continue;
            const int aim = static_cast<int>(std::lround(amrTargetLevel.value[ci]));
            if (aim > c.level())
                refineCount++;
            else if (aim < c.level())
                coarsenCount++;
        }
    }

    const bool needChange = (refineCount > 0) || (coarsenCount > 0);
    if (options.persistenceSteps <= 1)
    {
        amrPending = needChange;
        amrPendingCount = needChange ? 1 : 0;
    }
    else
    {
        if (needChange)
            amrPendingCount += 1;
        else
            amrPendingCount = 0;
        amrPending = (amrPendingCount >= options.persistenceSteps);
    }
}
