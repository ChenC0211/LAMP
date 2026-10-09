#pragma once

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include <filesystem>
#include <vector>

class QuadtreeMesh;
class ImplicitPicFieldSolve;
class ImplicitPicFieldHistory;
class PicPoissonFieldSolve;
class ExplicitPicChargeDeposition;
class ParticleType;
class CouplingWeightsCache;
class NormalizationContext;

// Live cell fields, bound to the session's mesh. Charge deposition writes
// through module operations; drivers borrow const charge/RHS fields. Field
// solvers and history still run at the original engine stages.
class PicFieldState
{
public:
    explicit PicFieldState(QuadtreeMesh& mesh);
    PicFieldState(const PicFieldState&) = delete;
    PicFieldState& operator=(const PicFieldState&) = delete;

    const CellField<double>& potential() const { return phi_; }
    const CellField<double>& poissonRhs() const { return rhs_; }
    const CellField<double>& chargeDensity() const { return rho_; }
    const CellField<double>& chargeAmount() const { return q_; }
    const CellField<Vector2D>& electricField() const { return electric_; }

    // Original implicit initialization/rebuild closure: full cell volume.
    // Explicit cut-volume deposition is owned by ExplicitPicChargeDeposition.
    void depositFullCellCharge(const std::vector<ParticleType>& species,
        const std::vector<double>& macroWeight, const std::vector<CouplingWeightsCache>& caches,
        const NormalizationContext& normalization);
    void copyChargeDensityToRhs();

    // AMR discards deposited fields. The implicit fast path may retain the
    // potential already transferred by mesh.adapt, including its cell labels.
    void resetAfterMeshChange(const CellField<double>* remappedPotential = nullptr);

private:
    friend class ImplicitPicRefinement;
    friend class ExplicitPicRefinement;
    friend class ImplicitPicFieldSolve;
    friend class ImplicitPicFieldHistory;
    friend class PicPoissonFieldSolve;
    friend class ExplicitPicChargeDeposition;
    const QuadtreeMesh& mesh_;
    CellField<double> phi_, rhs_, rho_, q_;
    CellField<Vector2D> electric_;
};

// Created at the original implicit-history initialization stage, after the
// initial Poisson solve. Initially zero; rebuild/restore are explicit actions.
// The bound state and mesh must outlive this object. Acceptance copies values
// into the current layer; it does not swap storage or overwrite the half layer.
class ImplicitPicFieldHistory
{
public:
    explicit ImplicitPicFieldHistory(PicFieldState& fields);
    ImplicitPicFieldHistory(const ImplicitPicFieldHistory&) = delete;
    ImplicitPicFieldHistory& operator=(const ImplicitPicFieldHistory&) = delete;

    const CellField<Vector2D>& nextElectricField() const { return next_; }
    const CellField<Vector2D>& halfElectricField() const { return half_; }
    void resetAfterMeshChange();
    void seedFromCurrent();
    void formHalfStep();
    void acceptNext();
    bool readCheckpoint(const std::filesystem::path& directory);
    void writeCheckpoint(const std::filesystem::path& directory) const;

private:
    friend class ImplicitPicRefinement;
    friend class PicPoissonFieldSolve;
    friend class ImplicitPicFieldSolve;
    void requireCurrentSizes() const;
    PicFieldState& fields_;
    CellField<Vector2D> next_, half_;
};
