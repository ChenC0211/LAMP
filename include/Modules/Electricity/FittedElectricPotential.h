#pragma once

#include "Core/Field/CellField.h"
#include <memory>
#include <stdexcept>
#include <utility>

namespace LevelSetHermiteGeometry { class Snapshot; }
namespace Quantum::Electricity { struct FittedCellHarmonicAccess; }
class FittedSharpDielectricSolver;
class FittedPotentialInput;

// Owning read-only cell potential. CellField's const unique_ptr still exposes
// writable elements, so no Field reference or raw data pointer is published.
// Publishing transfers storage; copying a result shares this immutable field.
class FittedElectricPotential
{
public:
    FittedElectricPotential() = default;
    int size() const { return field_ ? field_->size() : 0; }
    double operator[](int index) const { return field_->value[index]; }
    int label(int index) const { return field_->label(index); }
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source() const { return source_; }
    // AMR needs independent mutable interpolation storage, only after an
    // actual mesh-change request. Mutating this copy cannot alter the result.
    CellField<double> copyForMeshTransfer() const {
        if (!field_) throw std::logic_error("electric potential has not been published");
        return CellField<double>(*field_);
    }
private:
    friend class FittedSharpDielectricSolver;
    friend struct Quantum::Electricity::FittedCellHarmonicAccess;
    friend class FittedPotentialInput;
    FittedElectricPotential(std::shared_ptr<const CellField<double>> field,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source)
        : field_(std::move(field)), source_(std::move(source)) {}
    std::shared_ptr<const CellField<double>> field_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
};

// Synchronous read-only input: either an initial/AMR Field or a prior immutable
// result. Never retained by a solver. Only the solver can access Field storage;
// it makes the same single working copy as the original solve.
class FittedPotentialInput
{
public:
    FittedPotentialInput(const CellField<double>& field) : field_(&field) {}
    FittedPotentialInput(const FittedElectricPotential& field) : field_(field.field_.get()) {}
    int size() const { return field_ && field_->value ? field_->size() : 0; }
    double operator[](int index) const { return field_->value[index]; }
    int label(int index) const { return field_->label(index); }
private:
    friend class FittedSharpDielectricSolver;
    friend struct Quantum::Electricity::FittedCellHarmonicAccess;
    const CellField<double>* field_;
};
