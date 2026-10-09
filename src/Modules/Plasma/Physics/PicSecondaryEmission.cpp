#include "Modules/Plasma/Physics/PicSecondaryEmission.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include "Modules/Plasma/ExplicitPIC/Pusher/ParticlePusherLeapfrog.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Modules/Plasma/ExplicitPIC/Geometry/LevelSetQuery.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <utility>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
std::string lower(std::string s)
{
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
int sideIndex(ExplicitPIC::Pusher::DomainSide side)
{
    using Side = ExplicitPIC::Pusher::DomainSide;
    switch (side) {
    case Side::kLeft: return 0;
    case Side::kRight: return 1;
    case Side::kBottom: return 2;
    case Side::kTop: return 3;
    default: return -1;
    }
}
}

PicSecondaryEmission::PicSecondaryEmission(SEEConfig config, Mode mode)
    : config_(std::move(config)), mode_(mode), yieldModel_(lower(config_.gammaModel)),
      energyModel_(lower(config_.emittedEnergyModel)), angular_(lower(config_.angular)),
      monoJ_(config_.emittedEnergy_eV * PhysicalConstants::eV_J),
      thermalJ_(config_.emittedTe_eV * PhysicalConstants::eV_J),
      meanJ_(config_.emittedEmean_eV * PhysicalConstants::eV_J)
{}

void PicSecondaryEmission::beginStep(int step)
{
    if (mode_ == Mode::Explicit) normal_.reset();
    stats_ = {};
    pendingSides_.clear();
    step_ = step;
}

double PicSecondaryEmission::yield(double ionEnergyEV) const
{
    if (yieldModel_ == "powerlaw")
    {
        const double Eeff = std::max(0.0, ionEnergyEV - std::max(0.0, config_.Emin_eV));
        if (Eeff <= 0.0)
            return 0.0;
        const double denom = std::max(1e-12, config_.Eref_eV);
        double g = config_.gamma0 * std::pow(Eeff / denom, config_.alpha);
        if (config_.gammaMax > 0.0)
            g = std::min(g, config_.gammaMax);
        return std::max(0.0, g);
    }
    if (yieldModel_ == "piecewise")
    {
        const auto& E = config_.piecewiseE_eV;
        const auto& G = config_.piecewiseGamma;
        const size_t n = std::min(E.size(), G.size());
        if (n < 2)
            return std::max(0.0, config_.gamma); // fallback

        // Retain the two original endpoint rules for unbalanced legacy tables.
        const double endE = mode_ == Mode::Explicit ? E.back() : E[n - 1];
        const double endG = mode_ == Mode::Explicit ? G.back() : G[n - 1];
        if (ionEnergyEV <= E.front())
            return std::max(0.0, G.front());
        if (ionEnergyEV >= endE)
            return std::max(0.0, endG);

        // Find interval [i, i+1] containing ionEnergyEV (linear scan; n is tiny).
        for (size_t i = 0; i + 1 < n; ++i)
        {
            if (ionEnergyEV < E[i + 1])
            {
                const double t = (E[i + 1] > E[i]) ? ((ionEnergyEV - E[i]) / (E[i + 1] - E[i])) : 0.0;
                const double g = (1.0 - t) * G[i] + t * G[i + 1];
                return std::max(0.0, g);
            }
        }
        return std::max(0.0, endG);
    }
    // constant (default)
    return std::max(0.0, config_.gamma);
}

int PicSecondaryEmission::sampleCount(double expectedMacro, Engine& rng, Uniform& uniform) const
{
    if (expectedMacro <= 0.0) return 0;
    int count = static_cast<int>(std::floor(expectedMacro));
    const double fraction = expectedMacro - static_cast<double>(count);
    // Keep the original draw even when the fractional part is exactly zero.
    if (uniform(rng) < fraction) count += 1;
    if (config_.maxEmitMacro > 0) count = std::min(count, config_.maxEmitMacro);
    return count;
}

double PicSecondaryEmission::sampleEnergyJ(Engine& rng, Uniform& uniform)
{
    // NOTE: sampled energies are in SI (J).
    if (energyModel_ == "maxwell")
    {
        const double kT = std::max(0.0, thermalJ_);
        if (kT <= 0.0)
            return 0.0;
        const double g1 = normal_(rng);
        const double g2 = normal_(rng);
        const double g3 = normal_(rng);
        const double chi2 = g1 * g1 + g2 * g2 + g3 * g3;
        return 0.5 * kT * chi2; // mean = 1.5 kT
    }
    if (energyModel_ == "exponential")
    {
        const double mean = std::max(0.0, meanJ_);
        if (mean <= 0.0)
            return 0.0;
        const double u = std::max(1e-12, uniform(rng));
        return -mean * std::log(u);
    }
    // mono (default)
    return std::max(0.0, monoJ_);
}

std::array<double,3> PicSecondaryEmission::sampleDirection(
    const std::array<double,3>& normal, Engine& rng, Uniform& uniform) const
{
    std::array<double, 3> n = normal;
    const double nmag = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (nmag <= 0.0)
        n = {1.0, 0.0, 0.0};
    else
    {
        n[0] /= nmag;
        n[1] /= nmag;
        n[2] /= nmag;
    }

    if (angular_ == "fixed")
        return n;

    const std::array<double, 3> t1{-n[1], n[0], 0.0};
    const std::array<double, 3> t2{0.0, 0.0, 1.0};

    constexpr double kPi = 3.141592653589793238462643383279502884;
    const double u1 = uniform(rng);
    const double u2 = uniform(rng);
    const double mu = (angular_ == "isotropic") ? u1 : std::sqrt(u1);
    const double phi = 2.0 * kPi * u2;
    const double sinT = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    const double c = std::cos(phi);
    const double s = std::sin(phi);

    std::array<double, 3> d{
        mu * n[0] + sinT * c * t1[0] + sinT * s * t2[0],
        mu * n[1] + sinT * c * t1[1] + sinT * s * t2[1],
        mu * n[2] + sinT * c * t1[2] + sinT * s * t2[2],
    };
    const double dn = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (dn > 0.0)
    {
        d[0] /= dn;
        d[1] /= dn;
        d[2] /= dn;
    }
    return d;
}

std::array<double,3> PicSecondaryEmission::sampleDirection(
    ExplicitPIC::Pusher::DomainSide side, Engine& rng, Uniform& uniform) const
{
    std::array<double, 3> n{0.0, 0.0, 0.0};
    switch (side)
    {
        case ExplicitPIC::Pusher::DomainSide::kLeft:
            n = {1.0, 0.0, 0.0};
            break;
        case ExplicitPIC::Pusher::DomainSide::kRight:
            n = {-1.0, 0.0, 0.0};
            break;
        case ExplicitPIC::Pusher::DomainSide::kBottom:
            n = {0.0, 1.0, 0.0};
            break;
        case ExplicitPIC::Pusher::DomainSide::kTop:
            n = {0.0, -1.0, 0.0};
            break;
        default:
            n = {1.0, 0.0, 0.0};
            break;
    }
    return sampleDirection(n, rng, uniform);
}

bool PicSecondaryEmission::embeddedNormal(const QuadtreeMesh& mesh, ImmersedGeometry<QuadtreeMesh>* iGeometry,
    const DictionaryGeometry& dictGeometry, const PicEmbeddedGeometry& embeddedGeometry,
    int bid, const ExplicitPIC::Pusher::EmbeddedImpact& imp, Vector2D& inwardNormal, double& gmag)
{
    const auto& cutMetrics = embeddedGeometry.metrics();
    int levelSetI = embeddedGeometry.levelSetIndex();
    auto bit = dictGeometry.boundaryById.find(bid);
    if (bit != dictGeometry.boundaryById.end() && bit->second.levelSetIndex >= 0)
        levelSetI = bit->second.levelSetIndex;
    if (!iGeometry || levelSetI < 0 || levelSetI >= static_cast<int>(iGeometry->levelSets.size()))
        return false;

    Vector2D grad = ExplicitPIC::Geometry::LevelSetQuery::gradAt(mesh, *iGeometry, levelSetI, imp.pos);
    gmag = std::sqrt(grad.x * grad.x + grad.y * grad.y);
    if (!(gmag > 0.0) && embeddedGeometry.metricsReady() && imp.cellLabel >= 0)
    {
        const int cellIndex = mesh.cells(imp.cellLabel).index();
        const auto& cm = cutMetrics.cellMetricsByIndex();
        if (cellIndex >= 0 && cellIndex < static_cast<int>(cm.size()))
        {
            grad = cm[static_cast<std::size_t>(cellIndex)].n_gamma;
            gmag = std::sqrt(grad.x * grad.x + grad.y * grad.y);
        }
    }
    if (!(gmag > 0.0))
        return false;

    inwardNormal = embeddedGeometry.insideIsNegative() ? (grad * (1.0 / gmag))
                                                            : (grad * (-1.0 / gmag));
    return true;
}

std::optional<PicSecondaryEmission::Spawn> PicSecondaryEmission::spawnDomain(const QuadtreeMesh& mesh,
    const ExplicitPIC::Pusher::DomainImpact& imp, const DomainBox& box,
    const Kinematics& settings, Engine& rng, Uniform& unif01)
{
    const double seeEmittedMass_kg = settings.massKg, emit_w_clamped = settings.weight;
    const double v0_mps = settings.velocityScale, seeInjectEps_code = settings.epsilon;
    const double x0 = box.min.x, y0 = box.min.y, x1 = box.max.x, y1 = box.max.y;
    const double mx = box.margin.x, my = box.margin.y;
    Spawn result;
    const double E_emit_J = result.energyJ = sampleEnergyJ(rng, unif01);
    const double v_emit_si = (E_emit_J > 0.0)
                                ? std::sqrt(2.0 * E_emit_J / seeEmittedMass_kg)
                                : 0.0;
    const double v_emit_code = (v0_mps > 0.0) ? (v_emit_si / v0_mps) : 0.0;
    Particle& np = result.particle;
    np.position = imp.particle.position;
    np.w = emit_w_clamped;
    // Inject slightly inside the domain.
    const double eps = std::max(1e-15, seeInjectEps_code);
    const auto emitSide = imp.side;
    if (emitSide == ExplicitPIC::Pusher::DomainSide::kLeft)
        np.position.x = x0 + std::max(eps, mx);
    else if (emitSide == ExplicitPIC::Pusher::DomainSide::kRight)
        np.position.x = x1 - std::max(eps, mx);
    else if (emitSide == ExplicitPIC::Pusher::DomainSide::kBottom)
        np.position.y = y0 + std::max(eps, my);
    else if (emitSide == ExplicitPIC::Pusher::DomainSide::kTop)
        np.position.y = y1 - std::max(eps, my);

    // Clamp the other coordinate inside the domain.
    np.position.x = std::min(std::max(np.position.x, x0 + mx), x1 - mx);
    np.position.y = std::min(std::max(np.position.y, y0 + my), y1 - my);

    const auto dir = sampleDirection(emitSide, rng, unif01);
    np.velocity = Vector2D(v_emit_code * dir[0], v_emit_code * dir[1]);
    np.vz = v_emit_code * dir[2];

    int cellL = -1;
    if (imp.particle.inCellL >= 0)
        cellL = mesh.findCellLabelFromHint(imp.particle.inCellL, np.position);
    if (cellL < 0)
        cellL = mesh.findCellLabel(np.position);
    if (cellL < 0)
        return std::nullopt;
    np.inCellL = cellL;
    np.inBackgroundCellL = mesh.findBackgroundCellLabel(np.position);
    return result;
}

std::optional<PicSecondaryEmission::Spawn> PicSecondaryEmission::spawnEmbedded(const QuadtreeMesh& mesh,
    const ExplicitPIC::Pusher::EmbeddedImpact& imp, const Vector2D& inwardNormal,
    const Kinematics& settings, Engine& rng, Uniform& unif01)
{
    const double seeEmittedMass_kg = settings.massKg, emit_w_clamped = settings.weight;
    const double v0_mps = settings.velocityScale, seeInjectEps_code = settings.epsilon;
    Spawn result;
    const double E_emit_J = result.energyJ = sampleEnergyJ(rng, unif01);
    const double v_emit_si = (E_emit_J > 0.0)
                                ? std::sqrt(2.0 * E_emit_J / seeEmittedMass_kg)
                                : 0.0;
    const double v_emit_code = (v0_mps > 0.0) ? (v_emit_si / v0_mps) : 0.0;
    Particle& np = result.particle;
    np.position = imp.pos + inwardNormal * std::max(1e-15, seeInjectEps_code);
    np.w = emit_w_clamped;

    const auto dir = sampleDirection({inwardNormal.x, inwardNormal.y, 0.0}, rng, unif01);
    np.velocity = Vector2D(v_emit_code * dir[0], v_emit_code * dir[1]);
    np.vz = v_emit_code * dir[2];

    int cellL = -1;
    if (imp.particle.inCellL >= 0)
        cellL = mesh.findCellLabelFromHint(imp.particle.inCellL, np.position);
    if (cellL < 0 && imp.cellLabel >= 0)
        cellL = mesh.findCellLabelFromHint(imp.cellLabel, np.position);
    if (cellL < 0)
        cellL = mesh.findCellLabel(np.position);
    if (cellL < 0)
        return std::nullopt;
    np.inCellL = cellL;
    np.inBackgroundCellL = mesh.findBackgroundCellLabel(np.position);
    return result;
}

void PicSecondaryEmission::initialize(const std::unordered_set<int>& boundaryIds,
    const std::vector<ParticleType>& species, const std::vector<double>& weights,
    const std::unordered_map<std::string,std::size_t>& index,
    double lengthScale, double velocityScale, bool debug)
{
    if (initialized_) throw std::logic_error("PIC emission initialization may only be attempted once");
    initialized_ = true;
    boundaryIds_ = boundaryIds;
    velocityScale_ = velocityScale;
    debug_ = mode_ == Mode::Implicit && debug;
    const char* prefix = mode_ == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]";
    if (mode_ == Mode::Explicit) {
        auto names = config_.cathodeSides;
        if (names.empty() && !config_.cathodeSide.empty()) names.push_back(config_.cathodeSide);
        if (names.empty()) names.push_back("left");
        for (const auto& name : names) {
            const auto token = lower(name);
            const int side = token == "left" ? 0 : token == "right" ? 1 : token == "bottom" ? 2 : token == "top" ? 3 : -1;
            if (side < 0) {
                std::cerr << "[ExplicitPIC][SEE] WARNING: unknown cathode side '" << name << "'; ignoring." << std::endl;
                continue;
            }
            sideMask_[static_cast<std::size_t>(side)] = true;
        }
        const bool haveSides = sideMask_[0] || sideMask_[1] || sideMask_[2] || sideMask_[3];
        enabled_ = config_.enabled && (haveSides || !boundaryIds_.empty());
        if (config_.enabled && !haveSides && boundaryIds_.empty())
            std::cerr << "[ExplicitPIC][SEE] WARNING: no valid cathode sides or boundaryIds configured; disabling SEE." << std::endl;
    } else {
        enabled_ = config_.enabled && !boundaryIds_.empty();
        if (config_.enabled && boundaryIds_.empty())
            std::cerr << "[ImplicitPIC][SEE] WARNING: no embedded boundaryId has SEE enabled; disabling SEE." << std::endl;
    }
    if (yieldModel_ == "constant" && config_.gamma <= 0.0) enabled_ = false;
    if (enabled_) {
        const auto found = index.find(config_.emittedSpecies);
        if (found == index.end()) {
            std::cerr << prefix << "[SEE] WARNING: emittedSpecies '" << config_.emittedSpecies
                      << "' is not defined; disabling SEE." << std::endl;
            enabled_ = false;
        } else {
            emittedIndex_ = found->second;
            macroWeight_ = weights[emittedIndex_];
            mass_ = species[emittedIndex_].mass();
            if (mode_ == Mode::Explicit) {
                epsilon_ = config_.injectEps_m / ((lengthScale > 0.0) ? lengthScale : 1.0);
                if (epsilon_ < 0.0) epsilon_ = 0.0;
            } else epsilon_ = std::max(0.0, config_.injectEps_m / lengthScale);
        }
    }
    if (debug_) {
        std::ostringstream oss;
        oss << "[ImplicitPIC][SEE][init] enabled=" << (enabled_ ? "true" : "false")
            << " cfgEnabled=" << (config_.enabled ? "true" : "false")
            << " boundaryCount=" << boundaryIds.size()
            << " gammaModel=" << config_.gammaModel << " gamma=" << config_.gamma
            << " emittedSpecies=" << config_.emittedSpecies
            << " emittedIndex=" << (emittedIndex_ == static_cast<std::size_t>(-1) ? -1 : static_cast<int>(emittedIndex_))
            << " emittedTe_eV=" << config_.emittedTe_eV << " injectEps_m=" << config_.injectEps_m
            << " maxEmitMacro=" << config_.maxEmitMacro << " boundaryIds={";
        bool first = true;
        for (const int id : boundaryIds) { if (!first) oss << ','; first = false; oss << id; }
        oss << '}';
        std::cout << oss.str() << std::endl;
    }
}

void PicSecondaryEmission::requireMode(Mode mode) const
{
    if (!initialized_ || mode_ != mode) throw std::logic_error("PIC emission stage requires matching initialized mode");
}

double PicSecondaryEmission::impactEnergyJ(const Particle& particle, double mass) const
{
    const double vx = particle.velocity.x * velocityScale_;
    const double vy = particle.velocity.y * velocityScale_;
    const double vz = particle.vz * velocityScale_;
    return 0.5 * mass * (vx * vx + vy * vy + vz * vz);
}

void PicSecondaryEmission::recordImpact(double weight, double energyJ)
{
    stats_.ionImpactMacroCathode += 1;
    stats_.ionImpactPhysCathode += weight;
    stats_.ionImpactEnergyJCathode += energyJ * weight;
}

PicSecondaryEmission::Burst PicSecondaryEmission::sampleBurst(double weight, double energyJ,
    int boundaryId, double gradientMagnitude, Engine& rng, Uniform& uniform)
{
    const double energyEV = PhysicalConstants::eV_J > 0.0 ? energyJ / PhysicalConstants::eV_J : 0.0;
    const double gamma = yield(energyEV);
    const double denominator = macroWeight_ > 0.0 ? macroWeight_ : 1.0;
    const double expected = gamma * weight / denominator;
    if (debug_ && debugBudget_ > 0) {
        --debugBudget_;
        std::cout << "[ImplicitPIC][SEE][yield] step=" << step_ << " boundaryId=" << boundaryId
                  << " gammaEff=" << gamma << " expectedMacro=" << expected
                  << " seeMacroWeight=" << macroWeight_ << " gradMag=" << gradientMagnitude << std::endl;
    }
    if (expected <= 0.0) return {};
    const int count = sampleCount(expected, rng, uniform);
    // Explicit checks mass AFTER count sampling; implicit checks it before
    // recording impacts/resolving normals. Do not combine these two gates.
    if (count <= 0 || (mode_ == Mode::Explicit && mass_ <= 0.0)) return {};
    const double emitWeight = count > 0 ? gamma * weight / count : denominator;
    const double clamped = emitWeight > 0.0 ? emitWeight : denominator;
    if (debug_ && debugBudget_ > 0) {
        --debugBudget_;
        std::cout << "[ImplicitPIC][SEE][emit] step=" << step_ << " boundaryId=" << boundaryId
                  << " emitCount=" << count << " emitWeight=" << clamped << std::endl;
    }
    return {count, clamped};
}

void PicSecondaryEmission::appendExplicit(const Spawn& spawned, int side, std::vector<Particle>& pending)
{
    pending.push_back(spawned.particle);
    pendingSides_.push_back(side);
    const double weight = spawned.particle.w, energyJ = spawned.energyJ;
    stats_.emittedMacro += 1;
    stats_.emittedPhys += weight;
    stats_.emittedEnergyJ += energyJ * weight;
    stats_.emittedEnergyJ2 += (energyJ * energyJ) * weight;
}

void PicSecondaryEmission::processExplicit(const BoundaryGeometry& g,
    const std::vector<ParticleType>& species, std::size_t incidentIndex, double defaultWeight,
    const ExplicitPIC::Pusher::DomainImpactBuffer& domainImpacts,
    const ExplicitPIC::Pusher::EmbeddedImpactBuffer* embeddedImpacts,
    const DomainBox& box, const std::vector<int>& lineToBoundaryId,
    std::vector<Particle>& pending, Engine& rng, Uniform& uniform)
{
    requireMode(Mode::Explicit);
    const auto& incident = species[incidentIndex];
    if (!(incident.charge() > 0.0)) return;
    for (const auto& impact : domainImpacts.impacts()) {
        const int side = sideIndex(impact.side);
        if (side < 0 || !sideMask_[static_cast<std::size_t>(side)]) continue;
        const double weight = impact.particle.w > 0.0 ? impact.particle.w : defaultWeight;
        const double energyJ = impactEnergyJ(impact.particle, incident.mass());
        recordImpact(weight, energyJ);
        if (!enabled_ || emittedIndex_ == static_cast<std::size_t>(-1) || macroWeight_ <= 0.0) continue;
        const auto burst = sampleBurst(weight, energyJ, 0, 0., rng, uniform);
        for (int i = 0; i < burst.count; ++i) {
            const auto spawned = spawnDomain(g.mesh, impact, box, {mass_,burst.weight,velocityScale_,epsilon_}, rng, uniform);
            if (!spawned) continue;
            appendExplicit(*spawned, side, pending);
        }
    }
    if (!embeddedImpacts) return;
    for (const auto& impact : embeddedImpacts->impacts()) {
        const int boundaryId = g.embedded.boundaryIdAt(g.mesh,g.geometry,lineToBoundaryId,impact.pos,impact.cellLabel);
        if (boundaryIds_.find(boundaryId) == boundaryIds_.end()) continue;
        const double weight = impact.particle.w > 0.0 ? impact.particle.w : defaultWeight;
        const double energyJ = impactEnergyJ(impact.particle, incident.mass());
        recordImpact(weight, energyJ);
        if (!enabled_ || emittedIndex_ == static_cast<std::size_t>(-1) || macroWeight_ <= 0.0) continue;
        Vector2D normal; double magnitude = 0.;
        if (!embeddedNormal(g.mesh,g.geometry,g.dictionary,g.embedded,boundaryId,impact,normal,magnitude)) continue;
        const auto burst = sampleBurst(weight, energyJ, boundaryId, magnitude, rng, uniform);
        for (int i = 0; i < burst.count; ++i) {
            const auto spawned = spawnEmbedded(g.mesh,impact,normal,{mass_,burst.weight,velocityScale_,epsilon_},rng,uniform);
            if (spawned) appendExplicit(*spawned,-1,pending);
        }
    }
}

void PicSecondaryEmission::processImplicitEmbedded(const BoundaryGeometry& g,
    const ExplicitPIC::Pusher::EmbeddedImpact& impact, int boundaryId,
    const std::vector<ParticleType>& species, std::size_t incidentIndex,
    double weight, double energyJ, ParticleSpawnBuffer& pending, Engine& rng, Uniform& uniform)
{
    requireMode(Mode::Implicit);
    const auto& incident = species[incidentIndex];
    const bool boundaryEligible = boundaryId > 0 && boundaryIds_.find(boundaryId) != boundaryIds_.end();
    if (debug_ && debugBudget_ > 0 && incident.charge() > 0.0 && boundaryEligible) {
        --debugBudget_;
        std::cout << "[ImplicitPIC][SEE][impact] step=" << step_ << " species=" << incident.name()
                  << " boundaryId=" << boundaryId << " seeEnabled=" << (enabled_ ? "true" : "false")
                  << " w=" << weight << " ionEnergy_eV="
                  << (PhysicalConstants::eV_J > 0.0 ? energyJ / PhysicalConstants::eV_J : 0.0) << std::endl;
    }
    if (incident.charge() <= 0.0 || !enabled_ || emittedIndex_ == static_cast<std::size_t>(-1) ||
        macroWeight_ <= 0.0 || mass_ <= 0.0 || !boundaryEligible) return;
    recordImpact(weight, energyJ);
    Vector2D normal; double magnitude = 0.;
    if (!embeddedNormal(g.mesh,g.geometry,g.dictionary,g.embedded,boundaryId,impact,normal,magnitude)) return;
    const auto burst = sampleBurst(weight, energyJ, boundaryId, magnitude, rng, uniform);
    for (int i = 0; i < burst.count; ++i) {
        const auto spawned = spawnEmbedded(g.mesh,impact,normal,{mass_,burst.weight,velocityScale_,epsilon_},rng,uniform);
        if (spawned) pending.Push(species[emittedIndex_].name(),spawned->particle);
    }
}

PicSpawnManagement::Stats PicSecondaryEmission::thinExplicit(const PicSpawnManagement& manager,
    std::vector<Particle>& pending, bool enabled, double probability, Engine& rng, Uniform& uniform)
{
    requireMode(Mode::Explicit);
    return manager.thin(pending,emittedIndex_,enabled,probability,rng,uniform,&pendingSides_);
}

void PicSecondaryEmission::accountRetained(const std::vector<Particle>& pending,
    const std::vector<ParticleType>& species, const std::vector<double>& weights)
{
    requireMode(Mode::Explicit);
    if (pending.size() != pendingSides_.size())
        throw std::logic_error("SEE accounting requires retained boundary provenance");
    stats_.emittedMacro = 0;
    stats_.emittedPhys = stats_.emittedEnergyJ = stats_.emittedEnergyJ2 = 0.;
    stats_.injectedCharge = {};
    if (!enabled_ || emittedIndex_ == static_cast<std::size_t>(-1) || pending.empty()) return;
    const double mass = emittedIndex_ < species.size() && species[emittedIndex_].mass() > 0.0
        ? species[emittedIndex_].mass() : mass_;
    const double defaultWeight = emittedIndex_ < weights.size() ? weights[emittedIndex_] : 1.0;
    double physical = 0., energy = 0., energy2 = 0.;
    for (std::size_t i = 0; i < pending.size(); ++i) {
        const auto& particle = pending[i];
        const double weight = particle.w > 0.0 ? particle.w : defaultWeight;
        const double ke = mass > 0.0 ? impactEnergyJ(particle,mass) : 0.0;
        physical += weight;
        energy += ke * weight;
        energy2 += (ke * ke) * weight;
        const int side = pendingSides_[i];
        if (side >= 0) stats_.injectedCharge[static_cast<std::size_t>(side)] += species[emittedIndex_].charge() * weight;
    }
    stats_.emittedMacro = static_cast<long long>(pending.size());
    stats_.emittedPhys = physical;
    stats_.emittedEnergyJ = energy;
    stats_.emittedEnergyJ2 = energy2;
}
