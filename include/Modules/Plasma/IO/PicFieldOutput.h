#pragma once

#include <string>
#include <vector>
#include <filesystem>
#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"

class QuadtreeMesh;
class FaceSegmentGraph;
class LevelSet;
class PicProfiles;
namespace Solvers { namespace Poisson { class CutCellMetrics; } }
template<class T> class VertField;
template<class T> class CellField;

// PIC output preparation and legacy packed-field formats. Raw serializers take
// already converted coordinates/values and do not modify physical state.
// Masks retain their original behavior: only a correctly sized mask filters.
namespace PicFieldOutput
{
struct ExplicitInput {
    const QuadtreeMesh& mesh;
    const FaceSegmentGraph& graph;
    const PicEmbeddedGeometry& embedded;
    const CellField<double>& potential;
    const CellField<double>& charge;
    const CellField<Vector2D>& electric;
    const std::vector<ParticleType>& species;
    const std::vector<double>& weights;
    const UnitsContext* units;
    const NormalizationContext* normalization;
    const PicCurrentWeights::Policy& coupling;
    const std::vector<std::string>& fields;
};
// Immediate explicit frame at the caller's original stage (including step 0).
// Profile1D is sampled before packed fields, using separate diagnostic caches;
// averaged profiles remain at their solve-stage caller. Only the profile uses
// writeSIProfile; packed fields always use SI. No input reference is retained.
void writeExplicit(int step, double dt, const ExplicitInput& input,
    PicProfiles& profiles, bool writeSIProfile);

// One implicit cut-cell diagnostic at the caller's original field time layer.
// Resolve the legacy path and create its parent before sampling normal fields;
// file-open failure warns and returns, directory/diagnostic errors propagate.
// The sample is borrowed synchronously, never queued across AMR or a new step.
void writeCutcellState(const std::filesystem::path& outputRoot,
    const std::filesystem::path& configuredPath, int step, double timeSI,
    bool onlyActive, const PicCutcellStateView::Input& input);

// Optional implicit mask diagnostic, ./output/meta/packed_mask_diag_<step>.csv.
// This legacy diagnostic writes raw mesh coordinates/potential despite the
// x_m/y_m column names. Keep default stream precision, no unit conversion and
// the original silent file-open failure; directory errors still propagate.
// validCells/potential must match the current mesh. Wrong-sized alpha is ignored.
void writeMaskDiagnostics(const QuadtreeMesh& mesh, int step,
    const std::vector<unsigned char>& validCells, const CellField<double>& potential,
    const FaceSegmentGraph& faceGraph, const Solvers::Poisson::CutCellMetrics* cutMetrics,
    const std::vector<double>* alpha, const std::vector<LevelSet>* levelSets,
    const std::vector<int>& levelSetIndices, bool insideIsNegative);

// ./output/fields/<name>.dat, vertices in native mesh-index order.
void writeExplicitVertices(const QuadtreeMesh& mesh, const std::string& nameNoExt,
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<std::string>& names, const std::vector<VertField<double>>& fields,
    const std::vector<unsigned char>* cellMask = nullptr);

// ./output/<name>.dat, vertices sorted by the original coordinate tolerance.
void writeImplicitVertices(const QuadtreeMesh& mesh, const std::string& nameNoExt,
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<std::string>& names, const std::vector<VertField<double>>& fields,
    const std::vector<unsigned char>* cellMask = nullptr);

// ./output/<name>.dat, four separate vertices per active cell, no averaging.
void writeImplicitCells(const QuadtreeMesh& mesh, const std::string& nameNoExt,
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<std::string>& names, const std::vector<CellField<double>>& fields,
    const std::vector<unsigned char>* cellMask = nullptr);
}
