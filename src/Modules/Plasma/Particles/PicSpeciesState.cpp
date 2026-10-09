#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/Physics/PicDomain.h"
#include "Modules/Plasma/IO/PicCheckpointReader.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/ExplicitPIC/Field/MatchedGather.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/IO/ConsoleControl.h"
#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <utility>
#include <iostream>

void initializePicSpecies(ParticleType& pt, const SpeciesConfig& sp,
    const NormalizationScales& sc, const QuadtreeMesh& mesh, std::mt19937_64& rng,
    bool fromCheckpoint, const PicInitialParticleExclusion& exclude)
{
    pt.particle.reserve(std::max(0, sp.numParticles));
    const double L0 = std::max(1e-300, sc.L0_m);
    const double T0 = std::max(1e-300, sc.T0_s);
    const double Phi0 = sc.Phi0_V;
    double kappa_s = 0.0;
    if (sp.m > 0.0)
        kappa_s = (sp.q / sp.m) * (Phi0 * T0 * T0 / (L0 * L0));
    pt.set_kappa_s(kappa_s);

    const double v0_mps = std::max(1e-300, sc.v0_mps);
    const Vector2D drift_code(sp.drift.x / v0_mps, sp.drift.y / v0_mps);
    const double kT_J = std::max(0.0, sp.Te_eV) * PhysicalConstants::eV_J;
    const double vStd_code = (sp.m > 0.0 && kT_J > 0.0)
                                 ? std::sqrt(kT_J / sp.m) / v0_mps
                                 : 0.0;
    std::normal_distribution<double> gauss(0.0, vStd_code);

    const double x0 = mesh.minCoord.x;
    const double x1 = mesh.maxCoord.x;
    const double y0 = mesh.minCoord.y;
    const double y1 = mesh.maxCoord.y;
    const double mx = 1e-12 + 1e-6 * (x1 - x0);
    const double my = 1e-12 + 1e-6 * (y1 - y0);
    std::uniform_real_distribution<double> unifX(x0 + mx, x1 - mx);
    std::uniform_real_distribution<double> unifY(y0 + my, y1 - my);

    if (!fromCheckpoint)
    {
        for (int p = 0; p < sp.numParticles; ++p)
        {
            const double px = unifX(rng);
            const double py = unifY(rng);
            const double vx = drift_code.x + (vStd_code > 0.0 ? gauss(rng) : 0.0);
            const double vy = drift_code.y + (vStd_code > 0.0 ? gauss(rng) : 0.0);
            const double vz = (vStd_code > 0.0 ? gauss(rng) : 0.0);
            const int inCellL = mesh.findCellLabel(Vector2D(px, py));
            if (inCellL < 0)
                continue;
            if (exclude && exclude(Vector2D(px, py), inCellL))
                continue;

            Particle part;
            part.position = Vector2D(px, py);
            part.velocity = Vector2D(vx, vy);
            part.vz = vz;
            part.w = sp.weight;
            part.inCellL = inCellL;
            part.inBackgroundCellL = mesh.findBackgroundCellLabel(Vector2D(px, py));
            pt.particle.push_back(part);
        }
    }
    pt.number = static_cast<int>(pt.particle.size());
}

void PicSpeciesState::configurePush(PushOptions options)
{
    if (pushOptions_) throw std::logic_error("PIC push policy is already configured");
    pushOptions_ = std::move(options);
}

void PicSpeciesState::configurePush(const PicDomain& domain, bool useEmbeddedGeometry,
    bool periodicX, bool periodicY, bool absorbEmbeddedCrossing)
{
    PushOptions options;
    const auto& embedded = domain.embedded();
    options.absorbAtDomain = true;
    options.periodicX = periodicX;
    options.periodicY = periodicY;
    if (domain.mode_ == PicDomain::Mode::Explicit)
    {
        options.absorbAtEmbedded = embedded.ready() && embedded.absorbing();
        if (options.absorbAtEmbedded)
        {
            options.embeddedLevelSetIds = embedded.levelSetIndices().empty()
                ? std::vector<int>{embedded.levelSetIndex()} : embedded.levelSetIndices();
            options.embeddedInsideIsNegative = embedded.insideIsNegative();
            options.absorbEmbeddedCrossing = absorbEmbeddedCrossing;
        }
    }
    else
    {
        options.absorbAtEmbedded = useEmbeddedGeometry && embedded.ready() && embedded.absorbing();
        options.absorbEmbeddedCrossing = useEmbeddedGeometry && embedded.ready() && absorbEmbeddedCrossing;
        options.embeddedLevelSetIds = embedded.levelSetIndices();
        options.embeddedInsideIsNegative = embedded.insideIsNegative();
    }
    configurePush(std::move(options));
}

const PicSpeciesState::PushOptions& PicSpeciesState::pushOptions() const
{
    if (!pushOptions_) throw std::logic_error("PIC push policy is not configured");
    return *pushOptions_;
}

void PicSpeciesState::preparePushStorage()
{
    if (!initialized_) throw std::logic_error("PIC species are not initialized");
    (void)pushOptions();
    pushGather_.resize(species_.size());
    pushResults_.resize(species_.size());
}

void PicSpeciesState::beginPush()
{
    preparePushStorage();
    for (auto& result : pushResults_)
    {
        result.domain.Clear();
        result.embedded.Clear();
    }
}

void PicSpeciesState::initializeBackwardHalfStep(QuadtreeMesh& mesh,
    const std::vector<CouplingWeightsCache>& weights, const CellField<Vector2D>& field,
    double dt, ImmersedGeometry<QuadtreeMesh>* geometry)
{
    if (!initialized_ || !halfStepInitializationAllowed_)
        throw std::logic_error("PIC backward half-step requires fresh, unadvanced particles");
    if (weights.size() < species_.size() || field.size() != mesh.cellIndexNum())
        throw std::invalid_argument("PIC backward half-step requires current fields and species weights");
    for (std::size_t si = 0; si < species_.size(); ++si)
        if (!species_[si].particle.empty() && weights[si].numParticles() != species_[si].particle.size())
            throw std::invalid_argument("PIC backward half-step particle/cache count mismatch");
    preparePushStorage();
    const auto& options = pushOptions();
    halfStepInitializationAllowed_ = false;
    for (std::size_t si = 0; si < species_.size(); ++si)
    {
        auto& sp = species_[si];
        if (sp.particle.empty()) continue;
        auto& gathered = pushGather_[si];
        ExplicitPIC::Field::MatchedGather::gatherVector2D(weights[si], field, gathered, false);
        ExplicitPIC::Pusher::ParticlePusherLeapfrog::initializeHalfStep<QuadtreeMesh>(
            mesh, sp, gathered, -dt, options, geometry);
    }
}

PicSpeciesState::RelabelOptions PicSpeciesState::RelabelOptions::fromExplicitConfig(
    std::string mode, std::string fallback)
{
    for (auto& c : mode) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (auto& c : fallback) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return {(mode == "localmapthenfindcell" || mode == "localmap"),
            (fallback == "findcelllabel" || fallback == "findcell")};
}

void PicSpeciesState::relabelAfterMeshChange(const QuadtreeMesh& mesh, const RelabelOptions& options)
{
    if (!initialized_) throw std::logic_error("PIC species are not initialized");
    for (auto& sp : species_)
    {
        for (auto& p : sp.particle)
        {
            int cellL = p.inCellL;
            if (options.useHint) cellL = mesh.findCellLabelFromHint(cellL, p.position);
            if (cellL < 0 && options.useFallback) cellL = mesh.findCellLabel(p.position);
            p.inCellL = cellL;
            p.inBackgroundCellL = (cellL >= 0) ? mesh.findBackgroundCellLabel(p.position) : -1;
        }
        sp.number = static_cast<int>(sp.particle.size());
    }
}

const PicSpeciesState::PushResult& PicSpeciesState::pushSpecies(std::size_t index,
    QuadtreeMesh& mesh, const CouplingWeightsCache& weights, const CellField<Vector2D>& field,
    double dt, ImmersedGeometry<QuadtreeMesh>* geometry,
    bool recordEmbedded)
{
    preparePushStorage();
    const auto& options = pushOptions();
    auto& sp = species_.at(index);
    auto& result = pushResults_.at(index);
    halfStepInitializationAllowed_ = false;
    result.domain.Clear();
    result.embedded.Clear();
    if (sp.particle.empty()) return result;
    auto& gathered = pushGather_[index];
    ExplicitPIC::Field::MatchedGather::gatherVector2D(weights, field, gathered, false);
    ExplicitPIC::Pusher::ParticlePusherLeapfrog::advance<QuadtreeMesh>(
        mesh, sp, gathered, dt, options, geometry, &result.domain,
        recordEmbedded ? &result.embedded : nullptr);
    return result;
}

void PicSpeciesState::pushAll(QuadtreeMesh& mesh, const std::vector<CouplingWeightsCache>& weights,
    const CellField<Vector2D>& field, double dt,
    ImmersedGeometry<QuadtreeMesh>* geometry)
{
    if (weights.size() < species_.size())
        throw std::invalid_argument("PIC push requires weights for every species");
    beginPush();
    halfStepInitializationAllowed_ = false;
    for (std::size_t index = 0; index < species_.size(); ++index)
        pushSpecies(index, mesh, weights[index], field, dt, geometry, true);
}

PicSpeciesState::PicSpeciesState(std::size_t capacity, unsigned int randomSeed)
    : rng_(randomSeed == 0u ? 12345u : randomSeed)
{
    species_.reserve(capacity);
    weights_.reserve(capacity);
    ions_.reserve(capacity);
    electrons_.reserve(capacity);
    neutrals_.reserve(capacity);
}

void PicSpeciesState::initialize(const std::vector<SpeciesConfig>& config, const NormalizationScales& scales,
    const QuadtreeMesh& mesh, bool fromCheckpoint, const PicInitialParticleExclusion& exclude,
    const Observer& observer)
{
    if (initialized_) throw std::logic_error("PIC species are already initialized");
    initialized_ = true;
    for (const auto& sp : config)
    {
        if (sp.q > 0.0) ions_.insert(sp.name);
        if (sp.q < 0.0) electrons_.insert(sp.name);
        if (std::abs(sp.q) < 1e-30) neutrals_.insert(sp.name);
        ParticleType pt(sp.name, sp.m, sp.q);
        initializePicSpecies(pt, sp, scales, mesh, rng_, fromCheckpoint, exclude);
        if (observer) observer(sp, pt);
        species_.push_back(std::move(pt));
        weights_.push_back(sp.weight);
    }
    halfStepInitializationAllowed_ = !fromCheckpoint;
}

std::size_t PicSpeciesState::cullEmbedded(const PicDomain& domain, const char* tag,
    const std::vector<unsigned char>* speciesMask, std::vector<long long>* macroLoss,
    std::vector<double>* physicalLoss)
{
    if (!initialized_) throw std::logic_error("PIC embedded cull requires initialized species");
    const bool implicit = domain.mode_ == PicDomain::Mode::Implicit;
    const auto result = ParticleManager::cullEmbedded(domain.mesh(), domain.geometry(), domain.embedded(),
        species_, weights_, implicit, speciesMask, macroLoss, physicalLoss);
    if (result.removed > 0 && ConsoleControl::enabled(implicit ? ConsoleControl::Level::kInfo : ConsoleControl::Level::kDebug))
        std::cout << (implicit ? "[ImplicitPIC]" : "[ExplicitPIC]") << "[EB] Culled " << result.removed
                  << " particles inside embedded conductor (" << tag
                  << "); remaining=" << result.kept << std::endl;
    return result.removed;
}

std::optional<int> PicSpeciesState::initializeForRun(const std::vector<SpeciesConfig>& config,
    const NormalizationScales& scales, const PicDomain& domain,
    const PicCheckpointReader* restart, bool restoreRng)
{
    const bool implicit = domain.mode_ == PicDomain::Mode::Implicit;
    const bool consoleDebug = !implicit && ConsoleControl::enabled(ConsoleControl::Level::kDebug);
    const auto& embedded = domain.embedded();
    initialize(config, scales, domain.mesh(), restart != nullptr,
        [&](const Vector2D& position, int cellLabel) {
            return (implicit || (embedded.ready() && embedded.absorbing())) &&
                embedded.containsPoint(domain.mesh(), domain.geometry(), position, cellLabel);
        },
        [&](const SpeciesConfig& sp, const ParticleType& pt) {
            if (consoleDebug)
                std::cout << "  init species: " << pt.name() << " N=" << pt.number
                          << " q[C]=" << pt.charge() << " m[kg]=" << pt.mass()
                          << " Te[eV]=" << sp.Te_eV
                          << " kappa_s=" << pt.kappa_s() << std::endl;
        });
    if (consoleDebug) std::cerr << "[ExplicitPIC] init: species ready.\n";
    if (!restart) cullEmbedded(domain, "init");
    if (consoleDebug) std::cerr << "[ExplicitPIC] init: cull done.\n";
    int startStep = 0;
    if (restart)
    {
        if (!restart->readParticles(domain.mesh(), species_, rng_, restoreRng)) return std::nullopt;
        cullEmbedded(domain, "restart");
        startStep = restart->nextStep();
        restart->report();
    }
    if (consoleDebug) std::cerr << "[ExplicitPIC] init: restart handling done.\n";
    return startStep;
}

void PicSpeciesState::rebuildIndex()
{
    if (!initialized_) throw std::logic_error("PIC species are not initialized");
    index_.clear();
    index_.reserve(species_.size());
    for (std::size_t si = 0; si < species_.size(); ++si)
        index_.emplace(species_[si].name(), si);
}

int PicSpeciesState::electronIndex() const
{
    for (std::size_t si = 0; si < species_.size(); ++si)
    {
        std::string nm = species_[si].name();
        for (auto& c : nm) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (nm == "e" || nm == "electron" || nm.find("electron") != std::string::npos)
            return static_cast<int>(si);
    }
    for (std::size_t si = 0; si < species_.size(); ++si)
    {
        if (species_[si].charge() < 0.0)
            return static_cast<int>(si);
    }
    return -1;
}
