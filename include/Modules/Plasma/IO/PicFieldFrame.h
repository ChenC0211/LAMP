#pragma once

#include "Core/Field/CellField.h"
#include "Modules/Plasma/IO/PicFieldSelection.h"
#include "Modules/Plasma/IO/PicOutputReconstruction.h"
#include <functional>

class ParticleType;
class NormalizationContext;
class CouplingWeightsCache;
class UnitsContext;
class LevelSet;
class Vector2D;

// One prepared display frame. Owns names, field values and masks; borrows the
// current mesh and immutable selection until assembly finishes. No file IO or
// simulation-state mutation. The caller still owns scheduling and publication.
class PicFieldFrame
{
public:
    using BuildCache = std::function<void(std::size_t, CouplingWeightsCache&)>;

    PicFieldFrame(const QuadtreeMesh& mesh, const PicFieldSelection& selection,
        std::vector<unsigned char> validCells, const CellField<double>& potential,
        bool cellConstant = false);

    // Explicit: reconstruct, convert/clamp, publish vertices. Implicit:
    // convert/clamp, fill marked cells, publish cells or reconstruct vertices.
    template<class Convert>
    void addScalar(const std::string& name, const CellField<double>& source,
        const Convert& codeToSI, bool clampNonNegative = false)
    {
        if (selection_.mode() == PicFieldMode::Explicit)
        {
            auto vf = PicOutputReconstruction::buildMaskedVertField(mesh_, source, &validCells);
            VertField<double> prepared = *vf;
            const auto negatives = PicOutputReconstruction::convertScalarValues(
                prepared, prepared, codeToSI, clampNonNegative);
            if (clampNonNegative && negatives > 0)
                reportClamped(negatives, name);
            names.push_back(name);
            vertices.push_back(prepared);
            return;
        }
        CellField<double> prepared(mesh_, 0.0);
        PicOutputReconstruction::convertScalarValues(source, prepared, codeToSI, clampNonNegative);
        PicOutputReconstruction::fillMarkedCells(mesh_, validCells, markedCells_, prepared);
        if (cellConstant_)
        {
            names.push_back(name);
            cells.push_back(std::move(prepared));
            return;
        }
        auto vf = PicOutputReconstruction::buildMaskedVertField(mesh_, prepared, &validCells);
        if (!vf)
            return;
        names.push_back(name);
        vertices.push_back(std::move(*vf));
    }

    // The builder is called synchronously and is never retained. Explicit
    // density/energy requests rebuild separately; implicit requests reuse one
    // cache per species for this call. Fast implicit energy uses host cells
    // even when a missing/wrong-sized deposited amount falls back to CSR.
    void addParticles(const std::vector<ParticleType>& species,
        const std::vector<double>& defaultWeights, const NormalizationContext& norm,
        const std::vector<double>* alpha, const BuildCache& buildCache,
        bool fastDensity = false, const std::vector<CellField<double>>* depositedAmounts = nullptr);

    // Read current physical fields, append only selected display fields. Each
    // call finishes preparation synchronously and retains no input references.
    void addElectrostatic(const CellField<double>& potential, const CellField<double>& charge,
        const UnitsContext& units, const NormalizationContext& norm);
    void addElectric(const CellField<Vector2D>& electric, const UnitsContext& units);
    // Null composite means the legacy zero field. Null levelSets means no
    // individual level-set fields; invalid indices are skipped in input order.
    void addLevelSets(const UnitsContext& units, const VertField<double>* composite,
        const std::vector<LevelSet>* levelSets = nullptr,
        const std::vector<int>& indices = {}, const std::vector<std::string>& levelSetNames = {});
    // Implicit diagnostic coefficients: preserve per-segment (not area-weighted)
    // averaging and the original defaults for missing cell/segment values.
    void addSusceptibility(const FaceSegmentGraph& faceGraph, const std::vector<double>& muCell,
        const std::vector<double>& kappaByCell, const std::vector<double>& kappaBySegment);

    // Call once after assembly. Explicit order is fixed; implicit selection
    // ordering precedes periodic vertex reconciliation (Y then X).
    void finish(bool periodicX, bool periodicY);
    bool cellConstant() const { return cellConstant_; }

    std::vector<std::string> names;
    std::vector<VertField<double>> vertices;
    std::vector<CellField<double>> cells;
    const std::vector<unsigned char> validCells;

private:
    static void reportClamped(std::uint64_t count, const std::string& name);
    const QuadtreeMesh& mesh_;
    const PicFieldSelection& selection_;
    const bool cellConstant_;
    std::vector<bool> markedCells_;
};
