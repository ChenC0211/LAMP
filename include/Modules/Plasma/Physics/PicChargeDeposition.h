#pragma once

#include "Core/Field/CellField.h"
#include <vector>

class PicFieldState;
class PicEmbeddedGeometry;
class NormalizationContext;
class ParticleType;
class CouplingWeightsCache;
class QuadtreeMesh;
template<class Mesh> class ImmersedGeometry;

// Explicit charge/number deposition and active-volume density conversion.
// Owns the electron amount, reusable particle buffers and last charge totals;
// writes the bound field state through the original DepositQK numerical kernel.
class ExplicitPicChargeDeposition
{
public:
    struct Totals { double grid = 0.0, particles = 0.0; };
    explicit ExplicitPicChargeDeposition(PicFieldState& fields);
    ExplicitPicChargeDeposition(const ExplicitPicChargeDeposition&) = delete;
    ExplicitPicChargeDeposition& operator=(const ExplicitPicChargeDeposition&) = delete;

    void deposit(const std::vector<ParticleType>& species, const std::vector<double>& macroWeight,
        const std::vector<CouplingWeightsCache>& caches, const NormalizationContext& normalization,
        const PicEmbeddedGeometry& embedded, const ImmersedGeometry<QuadtreeMesh>* geometry,
        bool collectElectronAmount, int electronIndex);
    const CellField<double>& electronAmount() const { return amrNeInstant; }
    const Totals& totals() const { return totals_; }
    // The field owner is reset separately at its original AMR stage. Charge
    // statistics retain their old lifetime; only the number field is reset.
    void resetAfterMeshChange();

private:
    friend class ExplicitPicRefinement;
    friend class PicPoissonFieldSolve;
    PicFieldState& fields_;
    CellField<double> amrNeInstant;
    std::vector<std::vector<double>> depositQpScratch, depositWpScratch;
    Totals totals_;
};
