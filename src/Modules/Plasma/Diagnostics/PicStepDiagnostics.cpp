#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/Boundary/PicFieldBoundaryState.h"
#include "Modules/Plasma/ExplicitPIC/Field/ElectricFieldEvaluator.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Base/Interpolation.h"
#include "Core/Geometry/LevelSet.h"
#include "Modules/Plasma/Diagnostics/PhysicsDiagnostics.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/Physics/ParticleSpawnBuffer.h"
#include "Modules/Plasma/Boundary/MovingEBController.h"
#include <algorithm>
#include <cmath>
#include <utility>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include "Core/Field/CellField.h"
#include "Core/Mesh/PeriodicFaceUtils.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Poisson/PoissonBC.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Particles/PicSpeciesState.h"
#include "Modules/Plasma/Particles/ParticleManager.h"
#include "Modules/Plasma/Physics/PhysicalConstants.h"
#include <cctype>

PicParticleWeightStatistics summarizePicParticleWeights(const QuadtreeMesh& mesh,
    const ParticleType& species, double defaultWeight)
{
    const std::size_t n_macro = species.particle.size();
    if (n_macro == 0) return {};
    auto percentile = [](std::vector<double>& v, double p) -> double {
        if (v.empty())
            return 0.0;
        const double clamped = std::min(1.0, std::max(0.0, p));
        const std::size_t n = v.size();
        const std::size_t k = static_cast<std::size_t>(std::llround(clamped * static_cast<double>(n - 1)));
        const auto it = v.begin() + static_cast<std::vector<double>::difference_type>(k);
        std::nth_element(v.begin(), it, v.end());
        return v[k];
    };

    const int cellCount = mesh.cellIndexNum();
    std::vector<double> weights;
    weights.reserve(n_macro);
    double w_min = std::numeric_limits<double>::max();
    double w_max = 0.0;
    double w_sum = 0.0;

    std::vector<int> cellMacro;
    std::vector<double> cellPhys;
    cellMacro.assign(cellCount, 0);
    cellPhys.assign(cellCount, 0.0);

    for (const auto& p : species.particle)
    {
        const double w = (p.w > 0.0) ? p.w : defaultWeight;
        weights.push_back(w);
        w_min = std::min(w_min, w);
        w_max = std::max(w_max, w);
        w_sum += w;

        if (p.inCellL >= 0)
        {
            const int ci = mesh.cells(p.inCellL).index();
            if (ci >= 0 && ci < cellCount)
            {
                cellMacro[ci] += 1;
                cellPhys[ci] += w;
            }
        }
    }

    const double w_mean = (n_macro > 0) ? (w_sum / static_cast<double>(n_macro)) : 0.0;
    const double w_p50 = percentile(weights, 0.50);
    const double w_p90 = percentile(weights, 0.90);
    const double w_p99 = percentile(weights, 0.99);

    int cellMacroMax = 0;
    int cellMacroIdx = -1;
    double cellPhysMax = 0.0;
    int cellPhysIdx = -1;
    for (int ci = 0; ci < cellCount; ++ci)
    {
        if (cellMacro[ci] > cellMacroMax)
        {
            cellMacroMax = cellMacro[ci];
            cellMacroIdx = ci;
        }
        if (cellPhys[ci] > cellPhysMax)
        {
            cellPhysMax = cellPhys[ci];
            cellPhysIdx = ci;
        }
    }

    return {n_macro, w_sum, w_min, w_mean, {w_p50, w_p90, w_p99}, w_max,
        cellMacroMax, cellPhysMax, cellMacroIdx, cellPhysIdx};
}

PicEmissionEnergy summarizePicEmissionEnergy(double physical, double energyJ, double energyJ2)
{
    double emittedEnergyMean_eV = 0.0;
    double emittedEnergyStd_eV = 0.0;
    if (physical > 0.0 && PhysicalConstants::eV_J > 0.0)
    {
        const double meanJ = energyJ / physical;
        const double meanE2 = energyJ2 / physical;
        const double varJ2 = std::max(0.0, meanE2 - meanJ * meanJ);
        emittedEnergyMean_eV = meanJ / PhysicalConstants::eV_J;
        emittedEnergyStd_eV = std::sqrt(varJ2) / PhysicalConstants::eV_J;
    }

    return {emittedEnergyMean_eV, emittedEnergyStd_eV};
}

std::vector<PicRfPhase> summarizePicRfPhases(const PicDomainBC& domain,
    const PicEmbeddedBC& embedded, double timeSI, const UnitsContext* outputUnits)
{
    std::vector<PicRfPhase> rows;
    auto sample = [&](const char* target, const PicBoundarySide::TimeFunction& tf) {
        if (!tf.enabled || tf.freq_hz <= 0.0) return;
        double phaseRad = 0.0, phaseDeg = 0.0;
        const double potential = PicFieldBoundaryState::evaluate(tf, timeSI, &phaseRad, &phaseDeg);
        rows.push_back({target, tf.type, tf.freq_hz, phaseRad, phaseDeg,
            outputUnits ? outputUnits->PotentialCodeToSI(potential) : potential});
    };
    sample("left", domain.left.timeDirichlet);
    sample("right", domain.right.timeDirichlet);
    sample("bottom", domain.bottom.timeDirichlet);
    sample("top", domain.top.timeDirichlet);
    if (embedded.enabled) sample("embedded", embedded.timeDirichlet);
    return rows;
}

PicBoundaryLosses::PicBoundaryLosses(std::size_t count, double velocityScale,
    bool periodicX, bool periodicY, Mode mode)
    : velocityScale_(velocityScale), periodicX_(periodicX), periodicY_(periodicY), mode_(mode)
{
    totals_.domainMacro.resize(count);
    totals_.domainPhysical.resize(count);
    totals_.domainEnergy.resize(count);
    totals_.embeddedMacro.resize(count);
    totals_.embeddedPhysical.resize(count);
    totals_.embeddedEnergy.resize(count);
}

PicBoundaryLosses::ImplicitWall PicBoundaryLosses::summarizeImplicitWall(
    const std::vector<ParticleType>& species, const QuadtreeMesh& mesh,
    const UnitsContext* units, const NormalizationContext* normalization, double dtSI) const
{
    if (mode_ != Mode::Implicit || species.size() != totals_.domainPhysical.size())
        throw std::invalid_argument("Implicit wall summary requires matching implicit boundary losses");
    const double slabCfg_m = normalization ? normalization->config().slabThickness_m : 1.0;
    const double slab_m = slabCfg_m > 0.0 ? slabCfg_m : 1.0;
    const double Lx_m = units ? units->LengthCodeToSI(mesh.maxCoord.x - mesh.minCoord.x)
                             : mesh.maxCoord.x - mesh.minCoord.x;
    const double Ly_m = units ? units->LengthCodeToSI(mesh.maxCoord.y - mesh.minCoord.y)
                             : mesh.maxCoord.y - mesh.minCoord.y;
    ImplicitWall wall;
    wall.areaM2 = {std::max(0.0, Ly_m * slab_m), std::max(0.0, Ly_m * slab_m),
                   std::max(0.0, Lx_m * slab_m), std::max(0.0, Lx_m * slab_m)};
    for (int side = 0; side < 4; ++side)
    {
        double qOutE_C = 0.0, qOutI_C = 0.0;
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const double n_phys = totals_.domainPhysical[si][side];
            const double q_C = species[si].charge() * n_phys;
            if (species[si].charge() < 0.0) qOutE_C += q_C;
            else if (species[si].charge() > 0.0) qOutI_C += q_C;
        }
        wall.electronCharge[side] = qOutE_C;
        wall.ionCharge[side] = qOutI_C;
        const double A = wall.areaM2[side];
        const double dt_out = std::max(1e-300, dtSI);
        if (A > 0.0)
        {
            wall.electronCurrent[side] = qOutE_C / (A * dt_out);
            wall.ionCurrent[side] = qOutI_C / (A * dt_out);
        }
    }
    return wall;
}

PicBoundaryEnergy PicBoundaryEnergy::sample(PicSpeciesCount out,
    double energyJ, double areaM2, double dtSI)
{
    const double powerW = energyJ / dtSI;
    const double meanEV = (out.physical > 0.0 && PhysicalConstants::eV_J > 0.0)
        ? energyJ / (out.physical * PhysicalConstants::eV_J) : 0.0;
    const double flux = areaM2 > 0.0 ? powerW / areaM2 : std::numeric_limits<double>::quiet_NaN();
    return {out, energyJ, powerW, meanEV, flux, areaM2};
}

PicBoundaryLosses::Sample PicBoundaryLosses::sample(const ParticleType& species,
    const Particle& particle, double defaultWeight) const
{
    double w = particle.w > 0.0 ? particle.w : defaultWeight;
    if (mode_ == Mode::Implicit && (!std::isfinite(w) || !(w > 0.0))) w = defaultWeight;
    double energy = 0.0;
    const double mass = species.mass();
    if (mass > 0.0 && (mode_ == Mode::Explicit || velocityScale_ > 0.0))
    {
        const double vx = particle.velocity.x * velocityScale_;
        const double vy = particle.velocity.y * velocityScale_;
        const double vz = particle.vz * velocityScale_;
        const double v2 = vx * vx + vy * vy + vz * vz;
        energy = 0.5 * mass * v2;
    }
    return {true, w, species.charge() * w, energy};
}

bool PicBoundaryLosses::countEnergy(const ParticleType& species, double weightedEnergy) const
{
    return mode_ == Mode::Explicit ||
        (species.mass() > 0.0 && velocityScale_ > 0.0 && std::isfinite(weightedEnergy));
}

PicBoundaryLosses::Sample PicBoundaryLosses::recordDomain(std::size_t si, const ParticleType& species,
    const Particle& particle, double defaultWeight, int side)
{
    const int sign = species.charge() < 0.0 ? -1 : (species.charge() > 0.0 ? 1 : 0);
    if (side < 0 || side > 3)
    {
        if (sign < 0) ++totals_.unknownE;
        else if (sign > 0) ++totals_.unknownI;
        return {false, 0.0, 0.0, 0.0};
    }
    if ((periodicX_ && side < 2) || (periodicY_ && side >= 2))
    {
        if (mode_ == Mode::Explicit)
        {
            if (sign < 0) ++totals_.periodicLeakE[side];
            else if (sign > 0) ++totals_.periodicLeakI[side];
        }
        return {false, 0.0, 0.0, 0.0};
    }
    const auto impact = sample(species, particle, defaultWeight);
    if (sign < 0) ++totals_.domainMacroE[side];
    else if (sign > 0) ++totals_.domainMacroI[side];
    ++totals_.domainMacro[si][side];
    totals_.domainPhysical[si][side] += impact.weight;
    // Implicit side charge/current uses the original later species reduction.
    if (mode_ == Mode::Explicit)
    {
        if (sign < 0) totals_.domainChargeE[side] += impact.chargeC;
        else if (sign > 0) totals_.domainChargeI[side] += impact.chargeC;
    }
    const double energy = impact.energyPerParticleJ * impact.weight;
    if (countEnergy(species, energy))
    {
        totals_.domainEnergyJ += energy;
        totals_.domainEnergy[si][side] += energy;
    }
    return impact;
}

PicBoundaryLosses::Sample PicBoundaryLosses::recordEmbedded(std::size_t si, const ParticleType& species,
    const Particle& particle, double defaultWeight)
{
    const auto impact = sample(species, particle, defaultWeight);
    ++totals_.embeddedMacro[si];
    totals_.embeddedPhysical[si] += impact.weight;
    const double energy = impact.energyPerParticleJ * impact.weight;
    if (countEnergy(species, energy))
    {
        totals_.embeddedEnergyJ += energy;
        totals_.embeddedEnergy[si] += energy;
    }
    return impact;
}

void PicBoundaryLosses::recordEmbeddedBoundary(int boundaryId, const Sample& impact,
    bool crossing, bool explicitClassification)
{
    if (boundaryId > 0)
    {
        totals_.embeddedChargeById[boundaryId] += impact.chargeC;
        if (mode_ == Mode::Implicit || explicitClassification) ++totals_.embeddedCounts.byId[boundaryId];
    }
    if (mode_ == Mode::Explicit && explicitClassification)
    {
        if (boundaryId <= 0) ++totals_.embeddedCounts.unknown;
        if (crossing) ++totals_.embeddedCounts.crossing;
        else ++totals_.embeddedCounts.inside;
    }
}

namespace
{
Boundary selectElectrode(std::string_view side, const Solvers::Poisson::DomainBC& boundaries)
{
    std::string name(side);
    for (char& c : name) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (name.empty() || name == "auto")
    {
        const std::array<const PicBoundarySide*, 4> sides{
            &boundaries.left, &boundaries.right, &boundaries.bottom, &boundaries.top};
        const std::array<Boundary, 4> labels{Boundary::west, Boundary::east, Boundary::south, Boundary::north};
        for (std::size_t i = 0; i < sides.size(); ++i)
            if (sides[i]->type == PicBCType::Dirichlet && sides[i]->timeDirichlet.enabled) return labels[i];
        for (std::size_t i = 0; i < sides.size(); ++i)
            if (sides[i]->type == PicBCType::Dirichlet) return labels[i];
        return Boundary::non;
    }
    if (name == "left" || name == "west") return Boundary::west;
    if (name == "right" || name == "east") return Boundary::east;
    if (name == "bottom" || name == "south") return Boundary::south;
    if (name == "top" || name == "north") return Boundary::north;
    return Boundary::non;
}
int electrodeIndex(Boundary boundary)
{
    switch (boundary)
    {
    case Boundary::west: return 0;
    case Boundary::east: return 1;
    case Boundary::south: return 2;
    case Boundary::north: return 3;
    default: return -1;
    }
}
}

PicElectrodeCurrent::PicElectrodeCurrent(std::string_view side, const Solvers::Poisson::DomainBC& boundaries)
    : boundary_(selectElectrode(side, boundaries)), sideIndex_(electrodeIndex(boundary_))
{}

void PicElectrodeCurrent::sampleNormalFieldCode(const FaceSegmentGraph& graph,
    const std::vector<double>& normalFields)
{
    normalCode_ = averagePicBoundaryNormalField(graph.segments(), normalFields, boundary_);
}
void PicElectrodeCurrent::sampleFieldSI(double normalFieldVpm, double voltageV)
{
    normalSI_ = normalFieldVpm;
    voltageSI_ = voltageV;
}
void PicElectrodeCurrent::sampleFieldSI(const PicDomainBC& boundaries, double timeSI,
    const UnitsContext* units, std::optional<double> fieldScaleWithoutUnits)
{
    const PicBoundarySide* side = nullptr;
    switch (boundary_)
    {
    case Boundary::west: side = &boundaries.left; break;
    case Boundary::east: side = &boundaries.right; break;
    case Boundary::south: side = &boundaries.bottom; break;
    case Boundary::north: side = &boundaries.top; break;
    default: break;
    }
    double voltageCode = 0.0;
    if (side) voltageCode = PicFieldBoundaryState::voltage(*side, timeSI);
    sampleFieldSI(units ? units->ElectricFieldCodeToSI(normalCode_)
                       : (fieldScaleWithoutUnits ? *fieldScaleWithoutUnits * normalCode_ : normalCode_),
        units ? units->PotentialCodeToSI(voltageCode) : voltageCode);
}
void PicElectrodeCurrent::sampleConductionSI(const std::vector<ParticleType>& species,
    const std::vector<std::array<double, 4>>& physicalOutflow, double areaM2, double dtS)
{
    if (sideIndex_ < 0) return;
    conductionSI_ = 0.0;
    if (areaM2 > 0.0 && dtS > 0.0)
    {
        double qOut_C = 0.0;
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const double n_phys = physicalOutflow[si][static_cast<std::size_t>(sideIndex_)];
            qOut_C += species[si].charge() * n_phys;
        }
        conductionSI_ = qOut_C / (areaM2 * dtS);
    }
}
void PicElectrodeCurrent::finishStep(double dtS)
{
    if (sideIndex_ < 0)
    {
        current_ = Sample{};
        return;
    }
    double displacement = 0.0;
    if (havePreviousNormal_ && dtS > 0.0)
        displacement = PhysicalConstants::epsilon0_Fpm * (normalSI_ - previousNormalSI_) / dtS;
    current_ = {true, voltageSI_, conductionSI_, displacement, conductionSI_ + displacement};
    previousNormalSI_ = normalSI_;
    havePreviousNormal_ = true;
}

void PicBoundaryWork::addCurrent(double voltageV, double currentA)
{
    if (trackPotentialCurrentSum_) result_.potentialCurrentW += voltageV * currentA;
    result_.powerW += (-voltageV) * currentA;
}
void PicBoundaryWork::sampleNormalFieldsSI(const FaceSegmentGraph& graph,
    const std::vector<double>& normalFields, const UnitsContext& units)
{
    const auto averages = averagePicBoundaryNormalFields(graph.segments(), normalFields);
    for (int s = 0; s < 4; ++s)
        sampleNormalFieldSI(static_cast<std::size_t>(s), units.ElectricFieldCodeToSI(averages[static_cast<std::size_t>(s)]));
}
void PicBoundaryWork::sampleVoltagesSI(const PicDomainBC& boundaries, const UnitsContext& units,
    std::optional<double> evaluateAtTimeSI)
{
    const auto voltage = [&](const PicBoundarySide& side) {
        if (side.type != PicBCType::Dirichlet) return 0.0;
        return units.PotentialCodeToSI(evaluateAtTimeSI
            ? PicFieldBoundaryState::voltage(side, *evaluateAtTimeSI) : side.value);
    };
    sampleDirichlet({boundaries.left.type == PicBCType::Dirichlet, boundaries.right.type == PicBCType::Dirichlet,
                    boundaries.bottom.type == PicBCType::Dirichlet, boundaries.top.type == PicBCType::Dirichlet});
    sampleVoltagesSI({voltage(boundaries.left), voltage(boundaries.right), voltage(boundaries.bottom), voltage(boundaries.top)});
}
void PicBoundaryWork::addEmbedded(const std::unordered_map<int, double>& chargeOutC,
    const PicEmbeddedBC& fallback, const PicEmbeddedBCMap& byId, const UnitsContext& units, double timeSI)
{
    addEmbedded(chargeOutC, [&](int boundaryId, double& potentialV) {
        double valueCode = 0.0;
        if (!PicFieldBoundaryState::embeddedVoltage(fallback, byId, boundaryId, timeSI, valueCode)) return false;
        potentialV = units.PotentialCodeToSI(valueCode);
        return true;
    });
}
void PicBoundaryWork::addEmbedded(const std::unordered_map<int, double>& chargeOutC,
    const PicFieldBoundaryState& boundaries, const UnitsContext& units)
{
    addEmbedded(chargeOutC, [&](int boundaryId, double& potentialV) {
        double valueCode = 0.0;
        if (!boundaries.embeddedVoltage(boundaryId, valueCode)) return false;
        potentialV = units.PotentialCodeToSI(valueCode);
        return true;
    });
}

void PicBoundaryWork::addDomain(const std::array<double, 4>& areaM2,
    const std::vector<ParticleType>& species,
    const std::vector<std::array<double, 4>>& physicalOutflow, const Injection* injection)
{
    if (!(dtS_ > 0.0)) return;
    for (int s = 0; s < 4; ++s)
    {
        const auto side = static_cast<std::size_t>(s);
        if (!dirichlet_[side]) continue;
        const double area = areaM2[side];
        if (!(area > 0.0)) continue;
        double qOut_C = 0.0;
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            const double n_phys = physicalOutflow[si][side];
            qOut_C += species[si].charge() * n_phys;
        }
        double Jcond_Am2;
        if (injection)
        {
            const double qInj_C = injection->sourceChargeC[side] + injection->secondaryChargeC[side];
            const double qNetOut_C = qOut_C - qInj_C;
            Jcond_Am2 = qNetOut_C / (area * dtS_);
        }
        else
            Jcond_Am2 = qOut_C / (area * dtS_);
        double Jdisp_Am2 = 0.0;
        if (havePreviousNormal_[side])
        {
            const double dEn = normalSI_[side] - previousNormalSI_[side];
            Jdisp_Am2 = PhysicalConstants::epsilon0_Fpm * dEn / dtS_;
        }
        const double Jtot_Am2 = Jcond_Am2 + Jdisp_Am2;
        const double V_V = voltageSI_[side];
        addCurrent(V_V, Jtot_Am2 * area);
    }
}
void PicBoundaryWork::finishStep()
{
    if (!(dtS_ > 0.0)) return;
    for (int s = 0; s < 4; ++s)
    {
        previousNormalSI_[static_cast<std::size_t>(s)] = normalSI_[static_cast<std::size_t>(s)];
        havePreviousNormal_[static_cast<std::size_t>(s)] = true;
    }
    result_.energyJ = result_.powerW * dtS_;
}

PicProcessSummary summarizePicProcesses(const PicProcessSummaryInput& input)
{
    const double currentDenom = (std::isfinite(input.conductionCurrent) ? std::abs(input.conductionCurrent) : 0.0)
                              + (std::isfinite(input.displacementCurrent) ? std::abs(input.displacementCurrent) : 0.0);
    const double displacementFraction = (currentDenom > 0.0 && std::isfinite(input.displacementCurrent))
                                      ? (std::abs(input.displacementCurrent) / currentDenom)
                                      : 0.0;

    std::string dominantSource = "none";
    double dominantSourceEnergyJ = 0.0;
    const double sourceDriveJ = std::max(0.0, input.driveEnergyJ);
    const double sourceInjectionJ = std::max(0.0, input.injectionEnergyJ);
    if (sourceDriveJ >= sourceInjectionJ && sourceDriveJ > 0.0)
    {
        dominantSource = "drive";
        dominantSourceEnergyJ = sourceDriveJ;
    }
    else if (sourceInjectionJ > 0.0)
    {
        dominantSource = "injection";
        dominantSourceEnergyJ = sourceInjectionJ;
    }
    const double totalSourceEnergyJ = sourceDriveJ + sourceInjectionJ;
    const double dominantSourceFraction = (totalSourceEnergyJ > 0.0)
                                        ? (dominantSourceEnergyJ / totalSourceEnergyJ)
                                        : 0.0;

    static const char* kBoundaryNames[5] = {"L", "R", "B", "T", "embedded"};
    std::array<double, 5> boundaryEnergyJ = {0.0, 0.0, 0.0, 0.0, 0.0};
    for (std::size_t si = 0; si < input.speciesCount; ++si)
    {
        for (int s = 0; s < 4; ++s)
            boundaryEnergyJ[static_cast<std::size_t>(s)] += input.domainAbsorptionBySpecies[si][static_cast<std::size_t>(s)];
        boundaryEnergyJ[4] += input.embeddedAbsorptionBySpecies[si];
    }
    std::size_t dominantBoundaryIndex = 0;
    for (std::size_t bi = 1; bi < boundaryEnergyJ.size(); ++bi)
    {
        if (boundaryEnergyJ[bi] > boundaryEnergyJ[dominantBoundaryIndex])
            dominantBoundaryIndex = bi;
    }
    const double dE_abs_total_J = input.domainAbsorptionJ + input.embeddedAbsorptionJ;
    const double dominantBoundaryEnergyJ = boundaryEnergyJ[dominantBoundaryIndex];
    const double dominantBoundaryFraction = (dE_abs_total_J > 0.0)
                                         ? (dominantBoundaryEnergyJ / dE_abs_total_J)
                                         : 0.0;
    const double dominantBoundaryPowerW = (input.timeWindowS > 0.0)
                                        ? (dominantBoundaryEnergyJ / input.timeWindowS)
                                        : 0.0;

    std::string dominantReactionChannel = "none";
    double dominantReactionEnergyJ = 0.0;
    for (const auto& kv : input.reactionEnergyByChannel)
    {
        if (kv.second > dominantReactionEnergyJ)
        {
            dominantReactionChannel = kv.first;
            dominantReactionEnergyJ = kv.second;
        }
    }
    const double dominantReactionFraction = (input.reactionEnergyLossJ > 0.0)
                                           ? (dominantReactionEnergyJ / input.reactionEnergyLossJ)
                                           : 0.0;
    const double dominantReactionPowerW = (input.timeWindowS > 0.0)
                                         ? (dominantReactionEnergyJ / input.timeWindowS)
                                         : 0.0;

    return {
        displacementFraction,
        std::move(dominantSource),
        dominantSourceFraction,
        kBoundaryNames[dominantBoundaryIndex],
        dominantBoundaryPowerW,
        dominantBoundaryFraction,
        std::move(dominantReactionChannel),
        dominantReactionPowerW,
        dominantReactionFraction
    };
}

PicGeometrySummary summarizePicGeometry(bool cutMetricsReady,
    const std::vector<double>& cutAlphaByCellIndex, int cellCount,
    const MovingEBController* movingEb, int step)
{
    std::size_t activeCutCellCount = 0;
    const bool hasAlphaGeom = cutMetricsReady &&
        (cutAlphaByCellIndex.size() == static_cast<std::size_t>(cellCount));
    if (hasAlphaGeom)
    {
        for (int ci = 0; ci < cellCount; ++ci)
        {
            const double alpha = cutAlphaByCellIndex[static_cast<std::size_t>(ci)];
            if (alpha > 1e-12 && alpha < (1.0 - 1e-12))
                ++activeCutCellCount;
        }
    }

    std::size_t movingRegionCount = 0;
    std::size_t movingLevelsetCount = 0;
    double movingRegionTotalAreaM2 = 0.0;
    double movingRegionTotalFlux = 0.0;
    double movingRegionMeanEnergyEV = 0.0;
    double movingRegionMeanCos = 1.0;
    double movingRegionMeanAbsSpeedMps = 0.0;
    double movingRegionMaxAbsSpeedMps = 0.0;
    double geometryGsaMm = 0.0;
    double geometryGadMm = 0.0;
    double geometryDaMm = 0.0;
    if (movingEb)
    {
        const auto geometryDiag = movingEb->geometryDiagnostics();
        geometryGsaMm = geometryDiag.gsa_mm;
        geometryGadMm = geometryDiag.gad_mm;
        geometryDaMm = geometryDiag.da_mm;

        const auto regionDiags = movingEb->activeRegionDiagnostics(step);
        const auto levelSetIndices = movingEb->activeMovingLevelSetIndices(step);
        movingRegionCount = regionDiags.size();
        movingLevelsetCount = levelSetIndices.size();

        double areaWeightedEnergyEV = 0.0;
        double areaWeightedCos = 0.0;
        double sumAbsSpeedMps = 0.0;
        for (const auto& regionDiag : regionDiags)
        {
            const double areaM2 = std::max(0.0, regionDiag.area_m2);
            const double absSpeedMps = std::abs(regionDiag.vn_mps);
            movingRegionTotalAreaM2 += areaM2;
            movingRegionTotalFlux += regionDiag.flux_1pm2s;
            areaWeightedEnergyEV += areaM2 * regionDiag.meanEnergy_eV;
            areaWeightedCos += areaM2 * regionDiag.meanCos;
            sumAbsSpeedMps += absSpeedMps;
            movingRegionMaxAbsSpeedMps = std::max(movingRegionMaxAbsSpeedMps, absSpeedMps);
        }
        if (movingRegionTotalAreaM2 > 0.0)
        {
            movingRegionMeanEnergyEV = areaWeightedEnergyEV / movingRegionTotalAreaM2;
            movingRegionMeanCos = areaWeightedCos / movingRegionTotalAreaM2;
        }
        else if (!regionDiags.empty())
        {
            double sumEnergyEV = 0.0;
            double sumCos = 0.0;
            for (const auto& regionDiag : regionDiags)
            {
                sumEnergyEV += regionDiag.meanEnergy_eV;
                sumCos += regionDiag.meanCos;
            }
            movingRegionMeanEnergyEV = sumEnergyEV / static_cast<double>(regionDiags.size());
            movingRegionMeanCos = sumCos / static_cast<double>(regionDiags.size());
        }
        if (!regionDiags.empty())
            movingRegionMeanAbsSpeedMps = sumAbsSpeedMps / static_cast<double>(regionDiags.size());
    }

    return {
        activeCutCellCount,
        movingRegionCount,
        movingLevelsetCount,
        movingRegionTotalAreaM2,
        movingRegionTotalFlux,
        movingRegionMeanEnergyEV,
        movingRegionMeanCos,
        movingRegionMeanAbsSpeedMps,
        movingRegionMaxAbsSpeedMps,
        geometryGsaMm,
        geometryGadMm,
        geometryDaMm
    };
}

double averagePicBoundaryNormalField(const std::vector<FaceSegment>& segments,
    const std::vector<double>& normalFields, Boundary boundary)
{
    const size_t n = std::min(segments.size(), normalFields.size());
    double sumEnA = 0.0;
    double sumA = 0.0;
    for (size_t si = 0; si < n; ++si)
    {
        const FaceSegment& s = segments[si];
        if (s.boundary != boundary)
            continue;
        const double Af = s.Af;
        if (!(Af > 0.0))
            continue;
        sumEnA += normalFields[si] * Af;
        sumA += Af;
    }
    return (sumA > 0.0) ? (sumEnA / sumA) : 0.0;
}

std::array<double, 4> averagePicBoundaryNormalFields(const std::vector<FaceSegment>& segments,
    const std::vector<double>& normalFields)
{
    const size_t n = std::min(segments.size(), normalFields.size());
    std::array<double, 4> sumEnA = {0.0, 0.0, 0.0, 0.0};
    std::array<double, 4> sumA = {0.0, 0.0, 0.0, 0.0};
    auto bToIndex = [](Boundary b) -> int {
        switch (b)
        {
        case Boundary::west:  return 0;
        case Boundary::east:  return 1;
        case Boundary::south: return 2;
        case Boundary::north: return 3;
        default:              return -1;
        }
    };
    for (size_t si = 0; si < n; ++si)
    {
        const FaceSegment& s = segments[si];
        const int idx = bToIndex(s.boundary);
        if (idx < 0)
            continue;
        const double Af = s.Af;
        if (!(Af > 0.0))
            continue;
        sumEnA[static_cast<std::size_t>(idx)] += normalFields[si] * Af;
        sumA[static_cast<std::size_t>(idx)] += Af;
    }

    std::array<double, 4> result;
    for (int s = 0; s < 4; ++s)
    {
        const double EnAvg_code = (sumA[static_cast<std::size_t>(s)] > 0.0)
                                      ? (sumEnA[static_cast<std::size_t>(s)] / sumA[static_cast<std::size_t>(s)])
                                      : 0.0;
        result[static_cast<std::size_t>(s)] = EnAvg_code;
    }
    return result;
}

PicPeriodicJumps measurePicPeriodicJumps(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph, const CellField<double>& phiIn,
    bool periodicX, bool periodicY, const Solvers::Poisson::CutCellMetrics* cutMetrics,
    PicPeriodicMatePolicy matePolicy, bool debug)
{
    auto boundaryName = [](Boundary b) -> const char* {
        switch (b)
        {
        case Boundary::west:  return "west";
        case Boundary::east:  return "east";
        case Boundary::south: return "south";
        case Boundary::north: return "north";
        default:              return "non";
        }
    };
    PicPeriodicJumps out;
    if (!(periodicX || periodicY))
        return out;
    std::unordered_set<std::uint64_t> seenPairsX;
    std::unordered_set<std::uint64_t> seenPairsY;
    auto packPair = [](const int a, const int b) -> std::uint64_t {
        const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
        const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
        return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
    };

    const bool hasCutPeriodicMask =
        (cutMetrics != nullptr) &&
        (cutMetrics->betaBySegment().size() == faceGraph.segments().size()) &&
        (cutMetrics->cellMetricsByIndex().size() == static_cast<std::size_t>(mesh.cellIndexNum()));
    constexpr double periodicOpenTol = 1e-6;

    for (std::size_t sid = 0; sid < faceGraph.segments().size(); ++sid)
    {
        const auto& seg = faceGraph.segments()[sid];
        if (seg.boundary == Boundary::non)
            continue;

        bool isPX = false;
        bool isPY = false;
        if (periodicX && (seg.boundary == Boundary::west || seg.boundary == Boundary::east))
            isPX = true;
        if (periodicY && (seg.boundary == Boundary::south || seg.boundary == Boundary::north))
            isPY = true;

        if (!(isPX || isPY))
            continue;

        const int kLab = seg.K;
        const int kIdx = mesh.cells(kLab).index();
        if (kIdx < 0)
            continue;
        const Vector2D shift = PeriodicFaceUtils::shiftToMate(mesh, seg.boundary);
        const int lLab = PeriodicFaceUtils::findMateCellLabel(mesh, seg, shift);
        if (lLab < 0)
            continue;
        const int lIdx = mesh.cells(lLab).index();
        if (lIdx < 0)
            continue;
        if (matePolicy == PicPeriodicMatePolicy::CleanMate
            && !PeriodicFaceUtils::isCleanCellMate(mesh, seg, lLab, shift))
            continue;

        if (hasCutPeriodicMask)
        {
            double beta = cutMetrics->betaBySegment()[sid];
            if (!std::isfinite(beta))
                beta = 1.0;
            beta = std::clamp(beta, 0.0, 1.0);
            if (!(beta > 0.0))
                continue;

            const auto& cm = cutMetrics->cellMetricsByIndex();
            const auto& cmK = cm[static_cast<std::size_t>(kIdx)];
            const auto& cmL = cm[static_cast<std::size_t>(lIdx)];
            double alphaK = cmK.alpha;
            double alphaL = cmL.alpha;
            if (!std::isfinite(alphaK))
                alphaK = 1.0;
            if (!std::isfinite(alphaL))
                alphaL = 1.0;
            alphaK = std::clamp(alphaK, 0.0, 1.0);
            alphaL = std::clamp(alphaL, 0.0, 1.0);
            if (!(alphaK > 0.0) || !(alphaL > 0.0))
                continue;
            if ((beta < (1.0 - periodicOpenTol)) || cmK.isCut || cmL.isCut)
                continue;
        }

        const std::uint64_t pairKey = packPair(kIdx, lIdx);
        if (isPX && !seenPairsX.insert(pairKey).second)
            isPX = false;
        if (isPY && !seenPairsY.insert(pairKey).second)
            isPY = false;
        if (!(isPX || isPY))
            continue;

        const double phiK = phiIn.value[kIdx];
        const double phiL = phiIn.value[lIdx];
        const double jump = std::abs(phiK - phiL);

        auto maybeUpdateDebug = [&](bool isPeriodicAxis, double currentMax, std::string& slot) {
            if (!debug || !isPeriodicAxis || !(jump >= currentMax))
                return;

            std::ostringstream oss;
            oss << "boundary=" << boundaryName(seg.boundary)
                << " sid=" << sid
                << " K(label=" << kLab << ", index=" << kIdx << ", center=("
                << mesh.cells(kLab).centerCoord().x << "," << mesh.cells(kLab).centerCoord().y << "))"
                << " L(label=" << lLab << ", index=" << lIdx << ", center=("
                << mesh.cells(lLab).centerCoord().x << "," << mesh.cells(lLab).centerCoord().y << "))"
                << " faceCenter=(" << seg.faceCenter.x << "," << seg.faceCenter.y << ")"
                << " phiK=" << phiK
                << " phiL=" << phiL
                << " jump=" << jump;

            if (hasCutPeriodicMask)
            {
                double beta = cutMetrics->betaBySegment()[sid];
                if (!std::isfinite(beta))
                    beta = 1.0;
                const auto& cm = cutMetrics->cellMetricsByIndex();
                const auto& cmK = cm[static_cast<std::size_t>(kIdx)];
                const auto& cmL = cm[static_cast<std::size_t>(lIdx)];
                oss << " beta=" << beta
                    << " alphaK=" << cmK.alpha
                    << " alphaL=" << cmL.alpha
                    << " cutK=" << (cmK.isCut ? 1 : 0)
                    << " cutL=" << (cmL.isCut ? 1 : 0)
                    << " gammaBoundaryK=" << cmK.gammaBoundaryId
                    << " gammaBoundaryL=" << cmL.gammaBoundaryId;
            }

            slot = oss.str();
        };

        if (isPX)
        {
            maybeUpdateDebug(true, out.maxJumpX, out.debugX);
            out.maxJumpX = std::max(out.maxJumpX, jump);
        }
        if (isPY)
        {
            maybeUpdateDebug(true, out.maxJumpY, out.debugY);
            out.maxJumpY = std::max(out.maxJumpY, jump);
        }
    }
    return out;
}

PicSpeciesBudget summarizePicSpeciesBudget(const PicSpeciesBudgetInput& input)
{
    long long macroOutDomain = 0;
    double physOutDomain = 0.0;
    for (int s = 0; s < 4; ++s)
    {
        macroOutDomain += static_cast<long long>(input.domainMacro[static_cast<std::size_t>(s)]);
        physOutDomain += input.domainPhysical[static_cast<std::size_t>(s)];
    }
    const long long macroResidual = input.begin.macro + input.source.macro + input.mccSpawn.macro + input.seeSource.macro
        - macroOutDomain - input.outEmbedded.macro - input.end.macro;
    const double physResidual = input.begin.physical + input.source.physical + input.mccSpawn.physical + input.seeSource.physical
        - physOutDomain - input.outEmbedded.physical - input.end.physical;
    const long long macroReactionNet = macroResidual - input.particleManagement.macro;
    const double physReactionNet = physResidual - input.particleManagement.physical;
    return {input.begin, input.source, input.mccSpawn, input.seeSource, {macroOutDomain, physOutDomain},
        input.outEmbedded, input.end, input.particleManagement, {macroReactionNet, physReactionNet}, {macroResidual, physResidual}};
}
PicDiagnosticWindow::Totals::Totals(std::size_t count)
    : speciesMacroBegin(count, 0ll),
      speciesPhysBegin(count, 0.0),
      sourceMacroBySpecies(count, 0ll),
      sourcePhysBySpecies(count, 0.0),
      mccMacroBySpecies(count, 0ll),
      mccPhysBySpecies(count, 0.0),
      seeMacroBySpecies(count, 0ll),
      seePhysBySpecies(count, 0.0),
      particleMgmtMacroLoss(count, 0ll),
      particleMgmtPhysLoss(count, 0.0),
      wallOutMacroBySpecies(count, {0, 0, 0, 0}),
      wallOutPhysBySpecies(count, {0.0, 0.0, 0.0, 0.0}),
      wallAbsEnergyJBySpecies(count, {0.0, 0.0, 0.0, 0.0}),
      ebOutPhysBySpecies(count, 0.0),
      ebOutMacroBySpecies(count, 0ll),
      ebAbsEnergyJBySpecies(count, 0.0)
{}

std::uint64_t countPicMacroParticles(const std::vector<ParticleType>& species)
{
    std::uint64_t total = 0;
    for (const auto& sp : species)
        total += static_cast<std::uint64_t>(sp.particle.size());
    return total;
}

void PicDiagnosticWindow::beginFromParticles(double timeCode, const std::vector<ParticleType>& species,
    const std::vector<double>& weights)
{
    if (active_) return;
    std::vector<long long> macro(species.size(), 0ll);
    std::vector<double> physical(species.size(), 0.0);
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto count = PicSpeciesCount::sample(species[si], si < weights.size() ? weights[si] : 1.0);
        macro[si] = count.macro;
        physical[si] = count.physical;
    }
    begin(timeCode, macro, physical);
}

void PicDiagnosticWindow::begin(double timeCode, const std::vector<long long>& macroBegin,
    const std::vector<double>& physicalBegin)
{
    if (active_) return;
    data_.energy.driveJ = 0.0;
    data_.energy.injectionJ = 0.0;
    data_.energy.domainAbsorptionJ = 0.0;
    data_.energy.embeddedAbsorptionJ = 0.0;
    data_.energy.reactionLossJ = 0.0;
    data_.energy.seeEmissionJ = 0.0;
    data_.energy.ionImpactJ = 0.0;
    std::fill(data_.sourceMacroBySpecies.begin(), data_.sourceMacroBySpecies.end(), 0ll);
    std::fill(data_.sourcePhysBySpecies.begin(), data_.sourcePhysBySpecies.end(), 0.0);
    std::fill(data_.mccMacroBySpecies.begin(), data_.mccMacroBySpecies.end(), 0ll);
    std::fill(data_.mccPhysBySpecies.begin(), data_.mccPhysBySpecies.end(), 0.0);
    std::fill(data_.seeMacroBySpecies.begin(), data_.seeMacroBySpecies.end(), 0ll);
    std::fill(data_.seePhysBySpecies.begin(), data_.seePhysBySpecies.end(), 0.0);
    std::fill(data_.particleMgmtMacroLoss.begin(), data_.particleMgmtMacroLoss.end(), 0ll);
    std::fill(data_.particleMgmtPhysLoss.begin(), data_.particleMgmtPhysLoss.end(), 0.0);
    for (auto& values : data_.wallOutMacroBySpecies)
        values = {0, 0, 0, 0};
    for (auto& values : data_.wallOutPhysBySpecies)
        values = {0.0, 0.0, 0.0, 0.0};
    for (auto& values : data_.wallAbsEnergyJBySpecies)
        values = {0.0, 0.0, 0.0, 0.0};
    std::fill(data_.ebOutPhysBySpecies.begin(), data_.ebOutPhysBySpecies.end(), 0.0);
    std::fill(data_.ebOutMacroBySpecies.begin(), data_.ebOutMacroBySpecies.end(), 0ll);
    std::fill(data_.ebAbsEnergyJBySpecies.begin(), data_.ebAbsEnergyJBySpecies.end(), 0.0);
    data_.reactionEnergyJByChannel.clear();
    active_ = true;
    startTimeCode_ = timeCode;
    data_.speciesMacroBegin = macroBegin;
    data_.speciesPhysBegin = physicalBegin;
}
void PicDiagnosticWindow::addEnergy(const Energy& energy)
{
    data_.energy.driveJ += energy.driveJ;
    data_.energy.injectionJ += energy.injectionJ;
    data_.energy.domainAbsorptionJ += energy.domainAbsorptionJ;
    data_.energy.embeddedAbsorptionJ += energy.embeddedAbsorptionJ;
    data_.energy.reactionLossJ += energy.reactionLossJ;
    data_.energy.seeEmissionJ += energy.seeEmissionJ;
    data_.energy.ionImpactJ += energy.ionImpactJ;
}
void PicDiagnosticWindow::addBoundary(const BoundarySample& boundary)
{
    for (std::size_t si = 0; si < speciesCount_; ++si)
    {
        for (int s = 0; s < 4; ++s)
        {
            data_.wallOutMacroBySpecies[si][static_cast<std::size_t>(s)] += boundary.domainMacro[si][static_cast<std::size_t>(s)];
            data_.wallOutPhysBySpecies[si][static_cast<std::size_t>(s)] += boundary.domainPhysical[si][static_cast<std::size_t>(s)];
            data_.wallAbsEnergyJBySpecies[si][static_cast<std::size_t>(s)] += boundary.domainEnergyJ[si][static_cast<std::size_t>(s)];
        }
        data_.ebOutPhysBySpecies[si] += boundary.embeddedPhysical[si];
        data_.ebOutMacroBySpecies[si] += boundary.embeddedMacro[si];
        data_.ebAbsEnergyJBySpecies[si] += boundary.embeddedEnergyJ[si];
        if (si < boundary.managementMacro.size())
            data_.particleMgmtMacroLoss[si] += boundary.managementMacro[si];
        if (si < boundary.managementPhysical.size())
            data_.particleMgmtPhysLoss[si] += boundary.managementPhysical[si];
    }
}
void PicDiagnosticWindow::addSpawns(Spawn kind, const ParticleSpawnBuffer& pending,
    const std::unordered_map<std::string, std::size_t>& speciesIndexByName,
    const std::vector<double>& macroWeight)
{
    auto& macro = kind == Spawn::Source ? data_.sourceMacroBySpecies : data_.mccMacroBySpecies;
    auto& physical = kind == Spawn::Source ? data_.sourcePhysBySpecies : data_.mccPhysBySpecies;
    for (const auto& it : pending.items())
    {
        auto mit = speciesIndexByName.find(it.species_name);
        if (mit == speciesIndexByName.end())
            continue;
        const std::size_t si = mit->second;
        const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        const double w = (it.particle.w > 0.0) ? it.particle.w : wDefault;
        macro[si] += 1;
        physical[si] += w;
    }
}
void PicDiagnosticWindow::addSecondary(std::size_t seeEmittedIndex, const std::vector<Particle>& seeSpawn,
    const std::vector<double>& macroWeight)
{
    const double wDefault = (seeEmittedIndex < macroWeight.size()) ? macroWeight[seeEmittedIndex] : 1.0;
    for (const auto& p : seeSpawn)
    {
        const double w = (p.w > 0.0) ? p.w : wDefault;
        data_.seeMacroBySpecies[seeEmittedIndex] += 1;
        data_.seePhysBySpecies[seeEmittedIndex] += w;
    }
}
void PicDiagnosticWindow::addReactionEnergy(const std::unordered_map<std::string, double>& energy)
{
    for (const auto& kv : energy)
        data_.reactionEnergyJByChannel[kv.first] += kv.second;
}
PicSpeciesBudget PicDiagnosticWindow::budget(std::size_t si, PicSpeciesCount end) const
{
    return summarizePicSpeciesBudget({
        {data_.speciesMacroBegin[si], data_.speciesPhysBegin[si]},
        {data_.sourceMacroBySpecies[si], data_.sourcePhysBySpecies[si]},
        {data_.mccMacroBySpecies[si], data_.mccPhysBySpecies[si]},
        {data_.seeMacroBySpecies[si], data_.seePhysBySpecies[si]},
        data_.wallOutMacroBySpecies[si], data_.wallOutPhysBySpecies[si],
        {data_.ebOutMacroBySpecies[si], data_.ebOutPhysBySpecies[si]}, end,
        {data_.particleMgmtMacroLoss[si], data_.particleMgmtPhysLoss[si]}});
}

void observeExplicitPicBulkEnergy(const std::vector<ParticleType>& species,
    const std::unordered_map<std::string, std::size_t>& index,
    const std::vector<double>& weights, double velocityScale, const PicBulkEnergyObserver& observer)
{
    if (!observer) return;
    for (const auto& sp : species)
    {
        const double wDefault = [&]() {
            auto it = index.find(sp.name());
            if (it == index.end() || it->second >= weights.size())
                return 1.0;
            return weights[it->second];
        }();
        for (const auto& p : sp.particle)
        {
            const double w = (p.w > 0.0) ? p.w : wDefault;
            observer(sp.name(), sp.charge(), ParticleManager::kineticEnergyEV(p, sp.mass(), velocityScale), w);
        }
    }
}

PicImplicitParticleSummary summarizeImplicitPicParticles(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, double v0_mps, const PicBulkEnergyObserver& observeEnergy)
{
    PicImplicitParticleSummary result;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const double q = sp.charge();
        const double m = sp.mass();
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        for (const auto& p : sp.particle)
        {
            const double w = (p.w > 0.0) ? p.w : w_default;
            result.macroCharge += q * w;
            if (observeEnergy)
                observeEnergy(sp.name(), q, ParticleManager::kineticEnergyEV(p, m, v0_mps), w);
            result.kinetic.addImplicit(p, m, q, w, v0_mps);
        }
    }
    return result;
}

PicKineticEnergy summarizeExplicitPicKineticEnergy(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight, double velocityScale)
{
    PicKineticEnergy result;
    const double v0_mps = std::max(1e-300, velocityScale);
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const double m_kg = sp.mass();
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        double sum = 0.0;
        for (const auto& p : sp.particle)
        {
            const double vx = p.velocity.x * v0_mps;
            const double vy = p.velocity.y * v0_mps;
            const double vz = p.vz * v0_mps;
            const double v2 = vx * vx + vy * vy + vz * vz;
            const double ke_J = 0.5 * m_kg * v2;
            const double w = (p.w > 0.0) ? p.w : w_default;
            sum += ke_J * w;
        }
        if (sp.charge() < 0.0)
            result.electronJ += sum;
        else if (sp.charge() > 0.0)
            result.ionJ += sum;
        else
            result.neutralJ += sum;
    }
    return result;
}

void PicKineticEnergy::addImplicit(const Particle& p, double m, double q, double w, double v0_mps)
{
    if (!(m > 0.0) || !(w > 0.0)) return;
    const double vx = p.velocity.x * v0_mps;
    const double vy = p.velocity.y * v0_mps;
    const double vz = p.vz * v0_mps;
    const double v2 = vx * vx + vy * vy + vz * vz;
    const double Ek = 0.5 * m * v2 * w;
    if (!std::isfinite(Ek)) return;
    if (std::abs(q) < 1e-30)
        neutralJ += Ek;
    else if (q < 0.0)
        electronJ += Ek;
    else
        ionJ += Ek;
}

void PicPowerWindow::add(const PicEnergyTransfers& sample)
{
    transfers_.driveJ += sample.driveJ;
    transfers_.injectionJ += sample.injectionJ;
    transfers_.domainAbsorptionJ += sample.domainAbsorptionJ;
    transfers_.embeddedAbsorptionJ += sample.embeddedAbsorptionJ;
    transfers_.reactionLossJ += sample.reactionLossJ;
    transfers_.seeEmissionJ += sample.seeEmissionJ;
    transfers_.ionImpactJ += sample.ionImpactJ;
}

PicPowerWindow::Report PicPowerWindow::sample(double time, double fieldEnergy, double kineticEnergy)
{
    Report report{};
    auto& balance = report.balance;
    const double interval = havePrevious_ ? (time - previousTime_) : 0.0;
    if (havePrevious_ && interval > 0.0)
        balance.fieldRate = (fieldEnergy - previousField_) / interval;
    previousField_ = fieldEnergy;
    previousTime_ = time;
    balance.absorbedJ = transfers_.domainAbsorptionJ + transfers_.embeddedAbsorptionJ;
    balance.totalJ = fieldEnergy + kineticEnergy;
    if (havePrevious_)
        balance.residualJ = (balance.totalJ - previousTotal_) - (transfers_.driveJ + transfers_.injectionJ - balance.absorbedJ);
    previousTotal_ = balance.totalJ;
    havePrevious_ = true;
    report.potentialEnergyJ = -transfers_.driveJ;
    report.potentialCurrentW = (interval > 0.0) ? (report.potentialEnergyJ / interval) : 0.0;
    return report;
}

PicEnergyBalance PicStepEnergyHistory::sample(double dt, double fieldEnergy, double kineticEnergy,
    double driveEnergy, double injectionEnergy, double domainAbsorption, double embeddedAbsorption)
{
    PicEnergyBalance balance;
    balance.totalJ = fieldEnergy + kineticEnergy;
    balance.absorbedJ = domainAbsorption + embeddedAbsorption;
    const bool havePrevious = std::isfinite(previousField_) && std::isfinite(previousKinetic_);
    const double previousTotal = havePrevious ? (previousField_ + previousKinetic_) : 0.0;
    balance.residualJ = havePrevious
        ? ((balance.totalJ - previousTotal) - (driveEnergy + injectionEnergy - balance.absorbedJ)) : 0.0;
    balance.fieldRate = (havePrevious && dt > 0.0) ? ((fieldEnergy - previousField_) / dt) : 0.0;
    balance.kineticRate = (havePrevious && dt > 0.0) ? ((kineticEnergy - previousKinetic_) / dt) : 0.0;
    balance.totalRate = balance.fieldRate + balance.kineticRate;
    previousField_ = fieldEnergy;
    previousKinetic_ = kineticEnergy;
    return balance;
}

PicExplicitCurrentBalance summarizeExplicitPicCurrentBalance(
    const ParticleSpawnBuffer& sourceSpawnPending,
    const std::unordered_map<std::string, std::size_t>& speciesIndexByName,
    const std::vector<ParticleType>& species, const std::vector<double>& macroWeight,
    const std::vector<std::array<double, 4>>& wallOutPhysBySpecies,
    const std::vector<double>& ebOutPhysBySpecies, double dtSI)
{
    double qInjSigned_C = 0.0;
    double qInjIonAbs_C = 0.0;
    double qInjElectronAbs_C = 0.0;
    for (const auto& it : sourceSpawnPending.items())
    {
        auto sit = speciesIndexByName.find(it.species_name);
        if (sit == speciesIndexByName.end())
            continue;
        const std::size_t si = sit->second;
        const double q = species[si].charge();
        const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        const double w = (it.particle.w > 0.0) ? it.particle.w : wDefault;
        const double qC = q * w;
        qInjSigned_C += qC;
        if (q > 0.0)
            qInjIonAbs_C += std::abs(qC);
        else if (q < 0.0)
            qInjElectronAbs_C += std::abs(qC);
    }

    double qOutDomainSigned_C = 0.0;
    double qOutDomainIonAbs_C = 0.0;
    double qOutDomainElectronAbs_C = 0.0;
    double qOutEbSigned_C = 0.0;
    double qOutEbIonAbs_C = 0.0;
    double qOutEbElectronAbs_C = 0.0;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const double q = species[si].charge();

        double nPhysDomain = 0.0;
        for (int s = 0; s < 4; ++s)
            nPhysDomain += wallOutPhysBySpecies[si][static_cast<std::size_t>(s)];
        const double qDomain = q * nPhysDomain;
        qOutDomainSigned_C += qDomain;
        if (q > 0.0)
            qOutDomainIonAbs_C += std::abs(qDomain);
        else if (q < 0.0)
            qOutDomainElectronAbs_C += std::abs(qDomain);

        const double qEb = q * ebOutPhysBySpecies[si];
        qOutEbSigned_C += qEb;
        if (q > 0.0)
            qOutEbIonAbs_C += std::abs(qEb);
        else if (q < 0.0)
            qOutEbElectronAbs_C += std::abs(qEb);
    }

    const double dt_si_diag = dtSI;
    const double invDt = (dt_si_diag > 0.0) ? (1.0 / dt_si_diag) : 0.0;
    const double Ii_inj_A = qInjIonAbs_C * invDt;
    const double Ie_inj_A = qInjElectronAbs_C * invDt;
    const double IinjMag_A = Ii_inj_A + Ie_inj_A;
    const double IoutDomain_A = (qOutDomainIonAbs_C + qOutDomainElectronAbs_C) * invDt;
    const double IwallEb_A = (qOutEbIonAbs_C + qOutEbElectronAbs_C) * invDt;
    const double deltaMag_A = IinjMag_A - IoutDomain_A - IwallEb_A;
    const double deltaMagRel = deltaMag_A / std::max(1e-30, Ii_inj_A);

    const double IinjSigned_A = qInjSigned_C * invDt;
    const double IoutDomainSigned_A = qOutDomainSigned_C * invDt;
    const double IwallEbSigned_A = qOutEbSigned_C * invDt;
    const double deltaSigned_A = IinjSigned_A - IoutDomainSigned_A - IwallEbSigned_A;
    const double deltaSignedRel = deltaSigned_A / std::max(1e-30, Ii_inj_A);

    return {qInjSigned_C, qOutDomainSigned_C, qOutEbSigned_C, Ii_inj_A, Ie_inj_A, IinjMag_A, IoutDomain_A, IwallEb_A, deltaMag_A, deltaMagRel, IinjSigned_A, IoutDomainSigned_A, IwallEbSigned_A, deltaSigned_A, deltaSignedRel};
}

double sumPicMacroCharge(const std::vector<ParticleType>& species,
    const std::vector<double>& macroWeight)
{
    double qSum = 0.0;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const double q = species[si].charge();
        const double wDefault = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        for (const auto& p : species[si].particle)
        {
            const double w = (p.w > 0.0) ? p.w : wDefault;
            qSum += q * w;
        }
    }
    return qSum;
}

PicStepParticleBalance::PicStepParticleBalance(double beginC, std::size_t speciesCount)
    : begin_(speciesCount)
{
    result_.macroBeginC = beginC;
    for (auto& totals : spawns_) totals.species.resize(speciesCount);
}

PicStepParticleBalance::PicStepParticleBalance(const std::vector<ParticleType>& species,
    const std::vector<double>& weights)
    : PicStepParticleBalance(sumPicMacroCharge(species,weights),species.size())
{
    for (std::size_t i = 0; i < species.size(); ++i)
        begin_[i] = PicSpeciesCount::sample(species[i],i < weights.size() ? weights[i] : 1.0);
}

PicSpeciesCount PicSpeciesCount::sample(const ParticleType& species, double defaultWeight)
{
    PicSpeciesCount result{static_cast<long long>(species.particle.size()),0.};
    for (const auto& p : species.particle)
        result.physical += p.w > 0.0 ? p.w : defaultWeight;
    return result;
}

void PicStepParticleBalance::addSpawns(Spawn kind, const ParticleSpawnBuffer& pending,
    const std::vector<ParticleType>& species, const std::vector<double>& weights,
    const std::unordered_map<std::string,std::size_t>& index, double velocityScale)
{
    auto& totals = spawns_.at(static_cast<std::size_t>(kind));
    for (const auto& item : pending.items()) {
        const auto found = index.find(item.species_name);
        if (found == index.end()) continue;
        const auto si = found->second;
        const double defaultWeight = si < weights.size() ? weights[si] : 1.0;
        const double weight = item.particle.w > 0.0 ? item.particle.w : defaultWeight;
        auto& count = totals.species.at(si);
        count.macro += 1;
        count.physical += weight;
        if (kind == Spawn::Collision) continue;
        const double charge = species[si].charge() * weight;
        if (kind == Spawn::Secondary) addSecondary(charge);
        else addSource(charge);
        const double energy = ParticleManager::kineticEnergyEV(item.particle,species[si].mass(),velocityScale) * PhysicalConstants::eV_J;
        totals.energyJ += energy * weight;
    }
}

void PicStepParticleBalance::finishStage(double qGrid)
{
    result_.stageGridC = qGrid;
    result_.stageResidualC = result_.stageGridC - result_.stageMacroC;
    result_.stageRelative = (std::abs(result_.stageMacroC) > 0.0)
                         ? (result_.stageResidualC / result_.stageMacroC)
                         : 0.0;
}

void PicStepParticleBalance::finish(double Qmacro)
{
    result_.macroEndC = Qmacro;
    result_.expectedEndC = result_.macroBeginC + result_.sourceC - result_.domainC - result_.embeddedC;
    result_.residualC = result_.macroEndC - result_.expectedEndC;
    const double dQphysNorm =
        std::max(1e-30,
                 std::max(std::abs(result_.macroBeginC),
                          std::max(std::abs(result_.macroEndC),
                                   std::abs(result_.sourceC) +
                                       std::abs(result_.domainC) +
                                       std::abs(result_.embeddedC))));
    result_.relative = result_.residualC / dQphysNorm;

}

PicImplicitTimingSummary summarizePicImplicitTiming(const PicImplicitStepTiming& stepTiming)
{
    const double buildFallbackFrac =
        (stepTiming.buildWeightsN_particleCount > 0u)
            ? (static_cast<double>(stepTiming.buildWeightsN_cwFallbackCount) /
               static_cast<double>(stepTiming.buildWeightsN_particleCount))
            : 0.0;
    const double orbitAvgSampleFrac =
        (stepTiming.orbitAvg_sampleCount > 0u)
            ? (static_cast<double>(stepTiming.orbitAvg_sampleCwFallbackCount) /
               static_cast<double>(stepTiming.orbitAvg_sampleCount))
            : 0.0;
    const std::uint64_t totalFallbackCount =
        stepTiming.buildWeightsN_cwFallbackCount + stepTiming.orbitAvg_sampleCwFallbackCount;
    const std::uint64_t totalFallbackSamples =
        stepTiming.buildWeightsN_particleCount + stepTiming.orbitAvg_sampleCount;
    const double totalFallbackFrac =
        (totalFallbackSamples > 0u)
            ? (static_cast<double>(totalFallbackCount) / static_cast<double>(totalFallbackSamples))
            : 0.0;
    using NamedTime = PicImplicitTimingSummary::Stage;
    std::array<NamedTime, 8> topStageBuf = {{
        {"mcc_pre", stepTiming.mcc_pre_s},
        {"buildWeightsN", stepTiming.buildWeightsN_s},
        {"predictorSegment", stepTiming.predictorSegment_s},
        {"orbitAvgWeights", stepTiming.orbitAvgWeights_s},
        {"polPoisson_total", stepTiming.polPoisson_total_s},
        {"pm", stepTiming.pm_s},
        {"amr", stepTiming.amr_s},
        {"output", stepTiming.packedFieldOutput_s},
    }};
    std::sort(topStageBuf.begin(), topStageBuf.end(), [](const NamedTime& a, const NamedTime& b) {
        return a.t > b.t;
    });
    return {buildFallbackFrac, orbitAvgSampleFrac, totalFallbackFrac, totalFallbackCount,
        {topStageBuf[0], topStageBuf[1], topStageBuf[2]}};
}

PicCollisionElectronEnergy summarizePicCollisionElectronEnergy(
    const std::vector<ParticleType>& species, const std::vector<double>& weights, double velocityScale)
{
    constexpr double kTailThreshold_eV = 11.59; // Ar EXCITATION_1 threshold

    const double v0_mps = std::max(1e-300, velocityScale);
    double sum_w = 0.0;
    double sum_wE = 0.0;
    double sum_w_tail = 0.0;
    long long n_macro = 0;

    for (std::size_t si = 0; si < species.size(); ++si)
    {
        if (species[si].charge() >= 0.0)
            continue;
        const double m = species[si].mass();
        const double w_default = (si < weights.size()) ? weights[si] : 1.0;
        if (m <= 0.0)
            continue;

        for (const auto& p : species[si].particle)
        {
            const double vx = p.velocity.x * v0_mps;
            const double vy = p.velocity.y * v0_mps;
            if (!(std::isfinite(vx) && std::isfinite(vy)))
                continue;
            const double ke_J = 0.5 * m * (vx * vx + vy * vy);
            const double ke_eV = ke_J / PhysicalConstants::eV_J;
            if (!std::isfinite(ke_eV))
                continue;

            ++n_macro;
            double w = (p.w > 0.0) ? p.w : w_default;
            if (!std::isfinite(w) || !(w > 0.0))
                w = w_default;
            if (!std::isfinite(w) || !(w > 0.0))
                continue;
            sum_w += w;
            sum_wE += w * ke_eV;
            if (ke_eV > kTailThreshold_eV)
                sum_w_tail += w;
        }
    }

    const double mean_eV = (sum_w > 0.0) ? (sum_wE / sum_w) : 0.0;
    const double tail_frac = (sum_w > 0.0) ? (sum_w_tail / sum_w) : 0.0;
    const double n_phys = sum_w;
    return {mean_eV, tail_frac, n_macro, n_phys};
}

PicNeutralEnergy summarizePicNeutralEnergy(const ParticleType& pt, double defaultWeight, double velocityScale)
{
    const double m_kg = pt.mass();
    const double w_default = defaultWeight;

    const std::size_t n_macro = pt.particle.size();
    double n_phys = 0.0;
    double Etot_J = 0.0;
    for (const auto& p : pt.particle)
    {
        const double vx = p.velocity.x * velocityScale;
        const double vy = p.velocity.y * velocityScale;
        const double vz = p.vz * velocityScale;
        const double v2 = vx * vx + vy * vy + vz * vz;
        const double w = (p.w > 0.0) ? p.w : w_default;
        n_phys += w;
        Etot_J += 0.5 * m_kg * v2 * w;
    }
    const double meanE_eV = (n_phys > 0.0)
                               ? (Etot_J / n_phys) / PhysicalConstants::eV_J
                               : 0.0;

    return {n_macro, n_phys, Etot_J, meanE_eV};
}

PicConsoleValues summarizePicConsole(const QuadtreeMesh& mesh,
    const CellField<double>& phi, const UnitsContext* units,
    const std::vector<ParticleType>& species, const std::vector<double>& macroWeight,
    const MCCCollisionOperator::StepStats& mccStepStats,
    const PicConsoleValues::ProjectileGroups* configuredGroups)
{
    double phiMinCode = std::numeric_limits<double>::infinity();
    double phiMaxCode = -std::numeric_limits<double>::infinity();
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        phiMinCode = std::min(phiMinCode, phi.value[ci]);
        phiMaxCode = std::max(phiMaxCode, phi.value[ci]);
    }
    const double phiMinV = (units != nullptr) ? units->PotentialCodeToSI(phiMinCode) : phiMinCode;
    const double phiMaxV = (units != nullptr) ? units->PotentialCodeToSI(phiMaxCode) : phiMaxCode;

    const double tailThr_eV = PicConsoleValues::tailThresholdEV;
    const double v0_mps = units ? units->scales().v0_mps : 1.0;
    double electronNphys = 0.0;
    double electronEtot_J = 0.0;
    double electronTailNphys = 0.0;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        if (sp.charge() >= 0.0)
            continue;
        const double m_kg = sp.mass();
        const double w_default = (si < macroWeight.size()) ? macroWeight[si] : 1.0;
        for (const auto& p : sp.particle)
        {
            const double vx = p.velocity.x * v0_mps;
            const double vy = p.velocity.y * v0_mps;
            const double vz = p.vz * v0_mps;
            const double ke_J = 0.5 * m_kg * (vx * vx + vy * vy + vz * vz);
            const double ke_eV = ke_J / PhysicalConstants::eV_J;
            const double w = (p.w > 0.0) ? p.w : w_default;
            electronEtot_J += ke_J * w;
            electronNphys += w;
            if (ke_eV > tailThr_eV)
                electronTailNphys += w;
        }
    }
    const double eMeanE_eV = (electronNphys > 0.0)
                                 ? ((electronEtot_J / electronNphys) / PhysicalConstants::eV_J)
                                 : 0.0;
    const double Te_eV = (2.0 / 3.0) * eMeanE_eV;
    const double eTailFrac = (electronNphys > 0.0) ? (electronTailNphys / electronNphys) : 0.0;

    auto classifyChannelLower = [](const std::string& chanLower, std::uint64_t cnt,
                                   std::uint64_t* nElastic, std::uint64_t* nExc, std::uint64_t* nIoniz,
                                   std::uint64_t* nAttach, std::uint64_t* nCx, std::uint64_t* nOther) {
        if (!nElastic || !nExc || !nIoniz || !nAttach || !nCx || !nOther)
            return;
        if (chanLower == "elastic")
            *nElastic += cnt;
        else if (chanLower == "excitation" || chanLower.rfind("excitation_", 0) == 0)
            *nExc += cnt;
        else if (chanLower == "ionization")
            *nIoniz += cnt;
        else if (chanLower == "attachment")
            *nAttach += cnt;
        else if (chanLower == "cx" || chanLower == "charge_exchange" || chanLower == "chargeexchange" ||
                 chanLower == "charge-exchange")
            *nCx += cnt;
        else
            *nOther += cnt;
    };

    std::uint64_t nElastic = 0, nExc = 0, nIoniz = 0, nAttach = 0, nCx = 0, nOther = 0;
    for (const auto& kv : mccStepStats.per_channel)
    {
        const std::size_t pos = kv.first.rfind("::");
        const std::string chanRaw = (pos == std::string::npos) ? kv.first : kv.first.substr(pos + 2);
        classifyChannelLower(PhysicsDiagnostics::lowerCopy(chanRaw), kv.second,
                             &nElastic, &nExc, &nIoniz, &nAttach, &nCx, &nOther);
    }

    using MccSplit = PicConsoleValues::MccSplit;
    auto accProj = [&](MccSplit& a, const MCCCollisionOperator::StepStats::ProjectileStats& ps) {
        a.att += ps.n_attempts;
        a.coll += ps.n_collisions;
        a.nul += ps.n_null;
        a.spawn += ps.n_spawned;
        a.maxNuDt = std::max(a.maxNuDt, ps.max_nu_dt_majorant);
        for (const auto& ck : ps.per_channel)
            classifyChannelLower(PhysicsDiagnostics::lowerCopy(ck.first), ck.second,
                                 &a.el, &a.ex, &a.ion, &a.at, &a.cx, &a.other);
    };
    auto projectileGroup = [&](const std::string& name) -> char {
        if (configuredGroups)
            return configuredGroups->electrons.count(name) ? 'e' : configuredGroups->ions.count(name) ? 'i' : 'n';
        for (const auto& sp : species)
        {
            if (sp.name() != name)
                continue;
            if (sp.charge() < 0.0)
                return 'e';
            if (sp.charge() > 0.0)
                return 'i';
            return 'n';
        }
        return 'n';
    };

    MccSplit mccE, mccI, mccN;
    for (const auto& kv : mccStepStats.per_projectile)
    {
        const char group = projectileGroup(kv.first);
        if (group == 'e')
            accProj(mccE, kv.second);
        else if (group == 'i')
            accProj(mccI, kv.second);
        else
            accProj(mccN, kv.second);
    }

    return {phiMinV, phiMaxV, eMeanE_eV, Te_eV, eTailFrac, nElastic, nExc, nIoniz, nAttach, nCx, nOther, mccE, mccI, mccN};
}

PicExplicitConsoleValues summarizeExplicitPicConsole(const PicExplicitConsoleInput& input, double dtSI)
{
    const auto& species = input.particles.species();
    const auto& weights = input.particles.weights();
    const PicConsoleValues::ProjectileGroups groups{input.particles.electrons(), input.particles.ions()};
    PicExplicitConsoleValues values;
    values.state = summarizePicConsole(input.mesh, input.potential, input.units, species, weights, input.collisions, &groups);
    values.species.reserve(species.size());
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const auto& sp = species[si];
        const auto macro = static_cast<std::uint64_t>(sp.particle.size());
        const double fallback = si < weights.size() ? weights[si] : 1.0;
        double physical = 0.0;
        for (const auto& particle : sp.particle)
            physical += particle.w > 0.0 ? particle.w : fallback;
        values.macroCount += macro;
        values.species.push_back({macro, physical});
    }
    values.chargeDifference = input.gridCharge - input.particleCharge;
    values.chargeRelative = std::abs(input.particleCharge) > 0.0
        ? values.chargeDifference / std::abs(input.particleCharge) : 0.0;
    const double widthSI = input.units ? input.units->LengthCodeToSI(input.widthCode) : input.widthCode;
    const double heightSI = input.units ? input.units->LengthCodeToSI(input.heightCode) : input.heightCode;
    const double areaLR = input.units ? heightSI * input.slabMeters : heightSI;
    const double areaBT = input.units ? widthSI * input.slabMeters : widthSI;
    for (int side = 0; side < 4; ++side)
    {
        const double area = side < 2 ? areaLR : areaBT;
        if (area > 0.0 && dtSI > 0.0)
        {
            values.electronCurrent[side] = input.losses.domainChargeE[side] / (area * dtSI);
            values.ionCurrent[side] = input.losses.domainChargeI[side] / (area * dtSI);
        }
    }
    for (int side = 0; side < 4; ++side)
        values.totalCurrent += values.electronCurrent[side] + values.ionCurrent[side];
    for (int side = 0; side < 4; ++side)
        values.wallLossMacro += input.losses.domainMacroE[side] + input.losses.domainMacroI[side];
    if (input.embeddedStatistics)
        values.embeddedAbsorbed = input.losses.embeddedCounts.inside + input.losses.embeddedCounts.crossing;
    return values;
}

PicSpeciesAmount summarizePicSpeciesAmount(const CellField<double>* amount,
    int cellCount, double physicalCount)
{
    double stage = 0.0;
    if (amount && amount->size() == cellCount)
        for (int ci = 0; ci < cellCount; ++ci) stage += amount->value[ci];
    double relative = 0.0;
    if (std::abs(physicalCount) > 0.0)
        relative = (stage - physicalCount) / physicalCount;
    else if (std::abs(stage) > 0.0)
        relative = std::numeric_limits<double>::infinity();
    return {stage, physicalCount, relative};
}

PicEmbeddedBandField summarizeExplicitPicEmbeddedBandField(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faces, const Solvers::Poisson::CutCellMetrics* metrics,
    int bandLayers, const CellField<Vector2D>& field, double fieldScale)
{
    PicEmbeddedBandField result{std::max(0, bandLayers), 0, 0.0, 0.0};
    if (!metrics)
        return result;
    if (field.size() != mesh.cellIndexNum() || field.value == nullptr)
        return result;

    const auto& cm = metrics->cellMetricsByIndex();
    const int nCell = mesh.cellIndexNum();
    if (cm.size() != static_cast<std::size_t>(nCell))
        return result;

    std::vector<unsigned char> band(static_cast<std::size_t>(nCell), 0u);
    std::vector<int> frontier;
    frontier.reserve(256);

    for (int ci = 0; ci < nCell; ++ci)
    {
        if (!cm[static_cast<std::size_t>(ci)].isCut)
            continue;
        band[static_cast<std::size_t>(ci)] = 1u;
        frontier.push_back(ci);
    }
    result.cutCells = static_cast<int>(frontier.size());
    if (frontier.empty())
    {
        result.rmsVpm = 0.0 * fieldScale;
        result.meanVpm = 0.0 * fieldScale;
        return result;
    }

    const int layers = result.layers;
    for (int layer = 0; layer < layers; ++layer)
    {
        std::vector<unsigned char> inFront(static_cast<std::size_t>(nCell), 0u);
        for (const int ci : frontier)
        {
            if (ci >= 0 && ci < nCell)
                inFront[static_cast<std::size_t>(ci)] = 1u;
        }

        std::vector<int> next;
        next.reserve(frontier.size() * 2 + 16);
        for (const auto& seg : faces.segments())
        {
            if (seg.boundary != Boundary::non || seg.L < 0)
                continue;
            const int kIndex = mesh.cells(seg.K).index();
            const int lIndex = mesh.cells(seg.L).index();
            if (kIndex < 0 || lIndex < 0 || kIndex >= nCell || lIndex >= nCell)
                continue;

            const bool kFront = inFront[static_cast<std::size_t>(kIndex)] != 0u;
            const bool lFront = inFront[static_cast<std::size_t>(lIndex)] != 0u;
            if (kFront && !band[static_cast<std::size_t>(lIndex)])
            {
                band[static_cast<std::size_t>(lIndex)] = 1u;
                next.push_back(lIndex);
            }
            if (lFront && !band[static_cast<std::size_t>(kIndex)])
            {
                band[static_cast<std::size_t>(kIndex)] = 1u;
                next.push_back(kIndex);
            }
        }
        frontier.swap(next);
        if (frontier.empty())
            break;
    }

    double meanCode = 0.0, rmsCode = 0.0;
    double sumE = 0.0;
    double sumE2 = 0.0;
    int nBand = 0;
    for (int ci = 0; ci < nCell; ++ci)
    {
        if (!band[static_cast<std::size_t>(ci)])
            continue;
        const Vector2D e = field.value[ci];
        const double em = std::sqrt(e.x * e.x + e.y * e.y);
        sumE += em;
        sumE2 += em * em;
        ++nBand;
    }
    result.cells = nBand;
    if (nBand > 0)
    {
        meanCode = sumE / static_cast<double>(nBand);
        rmsCode = std::sqrt(sumE2 / static_cast<double>(nBand));
    }
    result.rmsVpm = rmsCode * fieldScale;
    result.meanVpm = meanCode * fieldScale;
    return result;
}

PicEmbeddedBandField summarizePicEmbeddedBandField(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faces, bool useCutMetrics, const std::vector<double>& alpha,
    int bandLayers, const CellField<Vector2D>& field, double fieldScale)
{
    std::vector<int> bandCells;
    if (useCutMetrics && (alpha.size() == static_cast<std::size_t>(mesh.cellIndexNum())))
    {
        std::vector<unsigned char> inBand(static_cast<std::size_t>(mesh.cellIndexNum()), 0u);
        std::vector<int> frontier;
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const double a = alpha[static_cast<std::size_t>(ci)];
            if (std::isfinite(a) && (a < 1.0 - 1e-12))
            {
                inBand[static_cast<std::size_t>(ci)] = 1u;
                frontier.push_back(ci);
            }
        }

        const int layers = std::max(0, bandLayers);
        if (!frontier.empty() && layers > 0)
        {
            std::vector<std::vector<int>> adj(static_cast<std::size_t>(mesh.cellIndexNum()));
            for (const auto& seg : faces.segments())
            {
                if (seg.K < 0 || seg.L < 0)
                    continue;
                const int kIdx = mesh.cells(seg.K).index();
                const int lIdx = mesh.cells(seg.L).index();
                if (kIdx < 0 || lIdx < 0 || kIdx == lIdx)
                    continue;
                adj[static_cast<std::size_t>(kIdx)].push_back(lIdx);
                adj[static_cast<std::size_t>(lIdx)].push_back(kIdx);
            }

            std::vector<int> nextFrontier;
            for (int layer = 0; layer < layers; ++layer)
            {
                nextFrontier.clear();
                for (const int ci : frontier)
                {
                    for (const int nb : adj[static_cast<std::size_t>(ci)])
                    {
                        if (!inBand[static_cast<std::size_t>(nb)])
                        {
                            inBand[static_cast<std::size_t>(nb)] = 1u;
                            nextFrontier.push_back(nb);
                        }
                    }
                }
                frontier.swap(nextFrontier);
                if (frontier.empty())
                    break;
            }
        }

        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            if (inBand[static_cast<std::size_t>(ci)])
                bandCells.push_back(ci);
        }
    }

    const double E0 = fieldScale;
    double sumE = 0.0;
    double sumE2 = 0.0;
    for (const int ci : bandCells)
    {
        if (ci < 0 || ci >= mesh.cellIndexNum() || ci >= static_cast<int>(field.size()))
            continue;
        const Vector2D e = field.value[ci];
        const double Emag = std::sqrt(e.x * e.x + e.y * e.y) * E0;
        if (!std::isfinite(Emag))
            continue;
        sumE += Emag;
        sumE2 += Emag * Emag;
    }
    const double nCells = static_cast<double>(bandCells.size());
    const double Emean = (nCells > 0.0) ? (sumE / nCells) : 0.0;
    const double Erms = (nCells > 0.0) ? std::sqrt(sumE2 / nCells) : 0.0;
    return {std::max(0, bandLayers), static_cast<long long>(bandCells.size()), Erms, Emean};
}

namespace {
std::vector<double> cutcellNormalField(const PicCutcellStateView::Input& input)
{
    ExplicitPIC::Field::ElectricFieldEvaluator::Options options;
    if (input.skipMissingEmbeddedBC)
        options.embeddedBCMissingPolicy = ExplicitPIC::Field::ElectricFieldEvaluator::Options::EmbeddedBCMissingPolicy::Skip;
    return ExplicitPIC::Field::ElectricFieldEvaluator::computeEmbeddedNormalField(
        input.mesh, input.metrics, input.embedded, input.embeddedById, input.potential, options);
}
bool hasEnabledCutcellBoundary(const PicEmbeddedBCMap* map)
{
    if (!map) return false;
    for (const auto& entry : *map) if (entry.second.enabled) return true;
    return false;
}
}

PicCutcellStateView::PicCutcellStateView(const Input& input)
    : input_(input), normalField_(cutcellNormalField(input)),
      useEmbeddedMap_(hasEnabledCutcellBoundary(input.embeddedById)),
      useEmbeddedFallback_(input.embedded && input.embedded->enabled)
{}
int PicCutcellStateView::size() const { return input_.mesh.cellIndexNum(); }
int PicCutcellStateView::levelSetCount() const
{ return input_.levelSets ? static_cast<int>(input_.levelSets->size()) : 0; }
bool PicCutcellStateView::active(int ci) const
{ return input_.metrics ? input_.metrics->isActiveCellIndex(ci) : true; }
PicCutcellStateView::Row PicCutcellStateView::sample(int ci) const
{
    const auto& mesh = input_.mesh;
    const auto* cutMetricsPtr = input_.metrics;
    const auto* cmVec = cutMetricsPtr ? &cutMetricsPtr->cellMetricsByIndex() : nullptr;
    const auto& rhoStar = input_.density;
    const auto& kappaByCell = input_.kappa;
    const auto& phiNow = input_.potential;
    const auto& eNow = input_.electric;
    const auto& enGammaByCellIndex = normalField_;
    const auto* embeddedBCByIdPtr = input_.embeddedById;
    const auto* embeddedBCPtr = input_.embedded;
    const bool useEmbeddedMap = useEmbeddedMap_;
    const bool useEmbeddedFallback = useEmbeddedFallback_;
    const bool active = this->active(ci);
    const int cLabel = mesh.cellIndexToLabel[ci];
    const Cell& c = mesh.cells(cLabel);
    const Vector2D xc = c.centerCoord();
    const double areaCode = c.area();
    const double x_m = (input_.units != nullptr) ? input_.units->LengthCodeToSI(xc.x) : xc.x;
    const double y_m = (input_.units != nullptr) ? input_.units->LengthCodeToSI(xc.y) : xc.y;
    const double lenScale_m = (input_.units != nullptr) ? input_.units->LengthCodeToSI(1.0) : 1.0;
    const double area_m2 = (input_.units != nullptr)
        ? areaCode * lenScale_m * lenScale_m
        : areaCode;

    Solvers::Poisson::CutCellMetrics::CellMetric defaultMetric;
    const auto& cm = (cmVec && static_cast<std::size_t>(ci) < cmVec->size())
        ? (*cmVec)[static_cast<std::size_t>(ci)]
        : defaultMetric;

    const double rhoCode = (ci < rhoStar.size()) ? rhoStar.value[ci] : 0.0;
    const double rhoSI = (input_.normalization != nullptr)
        ? input_.normalization->RhoTildeToSI(rhoCode)
        : rhoCode;
    const double kappa = (static_cast<std::size_t>(ci) < kappaByCell.size())
        ? kappaByCell[static_cast<std::size_t>(ci)]
        : 1.0;
    const double phiCode = (ci < phiNow.size()) ? phiNow.value[ci] : 0.0;
    const double phiV = (input_.units != nullptr) ? input_.units->PotentialCodeToSI(phiCode) : phiCode;
    const Vector2D eCode = (ci < eNow.size()) ? eNow.value[ci] : Vector2D(0.0, 0.0);
    const double exVpm = (input_.units != nullptr) ? input_.units->ElectricFieldCodeToSI(eCode.x) : eCode.x;
    const double eyVpm = (input_.units != nullptr) ? input_.units->ElectricFieldCodeToSI(eCode.y) : eCode.y;
    const double enGammaCode =
        (static_cast<std::size_t>(ci) < enGammaByCellIndex.size())
            ? enGammaByCellIndex[static_cast<std::size_t>(ci)]
            : 0.0;
    const double enGammaVpm = (input_.units != nullptr)
        ? input_.units->ElectricFieldCodeToSI(enGammaCode)
        : enGammaCode;

    const PicEmbeddedBC* eb = nullptr;
    const char* ebSource = "none";
    if (useEmbeddedMap && cm.gammaBoundaryId > 0)
    {
        auto it = embeddedBCByIdPtr->find(cm.gammaBoundaryId);
        if (it != embeddedBCByIdPtr->end() && it->second.enabled)
        {
            eb = &it->second;
            ebSource = "map";
        }
    }
    if (!eb && useEmbeddedFallback)
    {
        eb = embeddedBCPtr;
        ebSource = "fallback";
    }

    const double aGamma_m = (input_.units != nullptr)
        ? input_.units->LengthCodeToSI(cm.A_gamma)
        : cm.A_gamma;
    const double dGamma_m = (input_.units != nullptr)
        ? input_.units->LengthCodeToSI(cm.d_gamma)
        : cm.d_gamma;

    std::vector<double> psi;
    if (input_.levelSets)
        for (const auto& levelSet : *input_.levelSets)
            psi.push_back(Interpolation::biLinearInterpolation(mesh, levelSet, c, xc));
    return {ci, cLabel, c.level(), static_cast<int>(c.boundary()), xc, x_m, y_m, areaCode, area_m2,
        active, cm, rhoCode, rhoSI, kappa, phiCode, phiV, eCode, exVpm, eyVpm, enGammaCode, enGammaVpm,
        eb, ebSource, aGamma_m, dGamma_m, std::move(psi)};
}
