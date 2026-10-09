#pragma once

#include <cstddef>
#include <vector>
#include "Modules/Plasma/ExplicitPIC/Deposit/CouplingWeightsCache.h"

class QuadtreeMesh;
class ParticleType;
class NormalizationContext;
template<class T> class CellField;

// Prepared display fields. Inputs are borrowed for the call; no RNG or
// physical state is changed, and field preparation performs no IO.
// Output fields must be distinct, zero-initialized fields on the current mesh.
// A supplied CSR cache must describe these particles in their current order.
namespace PicParticleFields
{
void depositNumberAmount(const ParticleType& species, double defaultWeight,
    const CouplingWeightsCache& cache, CellField<double>& amount);

// Null/wrong-sized alpha means full volume. Non-finite alpha retains the old
// full-volume fallback; finite values are clamped to [0,1].
void numberDensitySI(const QuadtreeMesh& mesh, const NormalizationContext& norm,
    const CellField<double>& amount, const std::vector<double>* alpha,
    CellField<double>& density);

// Preserve the different legacy output policies. Positive is the implicit
// path: after fallback it tests positivity only. FinitePositive is explicit.
enum class ElectronWeightPolicy { Positive, FinitePositive };
std::size_t electronEnergyCSR(const QuadtreeMesh& mesh, const ParticleType& species,
    double defaultWeight, const NormalizationContext& norm, const CouplingWeightsCache& cache,
    ElectronWeightPolicy policy, CellField<double>& meanEnergyEV, CellField<double>& temperatureEV);
void electronEnergyHostCells(const QuadtreeMesh& mesh, const ParticleType& species,
    double defaultWeight, const NormalizationContext& norm,
    CellField<double>& meanEnergyEV, CellField<double>& temperatureEV);

void particleWeights(const QuadtreeMesh& mesh, const std::vector<ParticleType>& species,
    const std::vector<double>& defaultWeights, CellField<double>& sum,
    CellField<double>& mean, CellField<double>& maximum, CellField<double>& count);
}
