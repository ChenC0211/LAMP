#pragma once

#include "Core/Field/VertField.h"
#include <cstdint>
#include <memory>
#include <vector>

class QuadtreeMesh;
class FaceSegmentGraph;
namespace Solvers { namespace Poisson { class CutCellMetrics; } }
template<class T> class CellField;

// Visualization-only reconstruction. Inputs must use the current mesh indices;
// the returned field owns its values and does not change the physical cell field.
// Conversion rules, field selection and output scheduling belong to the caller.
namespace PicOutputReconstruction
{
// Valid quadtree topology and current cell indices are required. Cut metrics,
// when present, decide activity even if the face graph has no segment there.
// Alpha is an optional additional filter; a wrong-sized vector is ignored.
std::vector<unsigned char> buildActiveCellMask(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& faceGraph, const Solvers::Poisson::CutCellMetrics* cutMetrics,
    const std::vector<double>* alpha = nullptr);

// Implicit display policy: identify connected zero-potential cells, then fill
// marked cells in each prepared scalar. Masks must match the current mesh.
// Retains vertex-neighbour multiplicity, in-place traversal and the 20-pass cap.
// This is display extrapolation, not a correction of the physical potential.
std::vector<bool> markZeroPotentialCells(const QuadtreeMesh& mesh,
    const CellField<double>& potential, const std::vector<unsigned char>& validCells);
void fillMarkedCells(const QuadtreeMesh& mesh, const std::vector<unsigned char>& validCells,
    const std::vector<bool>& markedCells, CellField<double>& prepared);

// Source and destination must have the same size on the current mesh; in-place
// conversion is supported. The caller chooses cell-before-interpolation or
// vertex-after-interpolation, preserving its original order and unit rule.
template<class ScalarField, class Convert>
std::uint64_t convertScalarValues(const ScalarField& source, ScalarField& destination,
    const Convert& codeToSI, bool clampNonNegative = false)
{
    std::uint64_t negCount = 0;
    for (int i = 0; i < source.size(); ++i)
    {
        double val = codeToSI(source.value[i]);
        if (clampNonNegative && val < 0.0)
        {
            val = 0.0;
            negCount += 1;
        }
        destination.value[i] = val;
    }
    return negCount;
}

// Only a correctly sized mask filters the stencil. Retains the original bounded
// least-squares fit, sparse-stencil fallbacks and zero for unsupported vertices.
std::unique_ptr<VertField<double>> buildMaskedVertField(const QuadtreeMesh& mesh,
    const CellField<double>& cells, const std::vector<unsigned char>* cellMask = nullptr);

// Reconciles matching boundary vertices in prepared output fields in place,
// Y before X (including corners). Unmatched vertices retain their values.
// This does not impose periodicity on simulation state or rebuild the mesh.
void canonicalizePeriodicBoundaryVertices(const QuadtreeMesh& mesh,
    bool periodicX, bool periodicY, std::vector<VertField<double>>& fields);
}
