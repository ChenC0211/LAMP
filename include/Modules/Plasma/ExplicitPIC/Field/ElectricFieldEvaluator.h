#pragma once

#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include "Solvers/Poisson/PoissonBC.h"

namespace ExplicitPIC
{
namespace Field
{

// ============================================================================
// ElectricFieldEvaluator
// ----------------------------------------------------------------------------
// 1) Face-normal evaluation:
//      E_{n,f} = -(phi_L - phi_K) / d_f
//    with boundary variants (Dirichlet/Neumann/Robin-as-Dirichlet).
//
// 2) Embedded-boundary normal evaluation (Robin BC):
//      a*phi + b*dphi/dn = g  on EB,
//      dphi/dn = (g - a*phi_K)/(a*d_gamma + b),  E_n = -dphi/dn.
//
// 3) Weighted least-squares reconstruction of cell-centered E_K from
//    (E_{n,f}, n_f) pairs. In the common axis-aligned case, the normal
//    equations are diagonal; we still support a full 2x2 solve when EB
//    normals introduce off-diagonal coupling.
//
// Weight model: omega_f = beta_f * A_f.
//   - beta_f from CutCellMetrics::betaBySegment() when provided; else 1.
//   - A_f and d_f are from FaceSegmentGraph::segments().
// ============================================================================

class ElectricFieldEvaluator {
public:
    struct Options
    {
        bool includeEmbedded = true;          // include EB equations in LS
        bool axisAlignedSimplification = true; // use diagonal solve when possible
        double eps = 1e-14;

        enum class EmbeddedBCMissingPolicy
        {
            Error,
            Skip
        };
        EmbeddedBCMissingPolicy embeddedBCMissingPolicy = EmbeddedBCMissingPolicy::Error;
    };

    // Workaround for GCC: avoid Options{} as a default argument inside the enclosing class.
    static Options defaultOptions() { return Options{}; }


    // Normal-field outputs.
    struct NormalFields
    {
        std::vector<double> EnBySegment;       // size = faceGraph.segments().size()
        std::vector<double> EnGammaByCellIndex; // size = mesh.cellIndexNum()
    };

    // Compute E_n per face segment (aligned with faceGraph.segments()).
    static std::vector<double> computeFaceNormalField(
        const QuadtreeMesh& mesh,
        const FaceSegmentGraph& faceGraph,
        const PicDomainBC& domainBC,
        const CellField<double>& phi,
        const Options& opt = defaultOptions());

    // Compute E_n on EB per cell index (size = mesh.cellIndexNum()).
    // Returns a vector filled with 0 for non-cut/non-EB cells.
    static std::vector<double> computeEmbeddedNormalField(
        const QuadtreeMesh& mesh,
        const Solvers::Poisson::CutCellMetrics* metrics,
        const PicEmbeddedBC* embeddedBC,
        const PicEmbeddedBCMap* embeddedBCById,
        const CellField<double>& phi,
        const Options& opt = defaultOptions());

    // Convenience: call computeFaceNormalField + computeEmbeddedNormalField.
    static NormalFields evaluate(
        const QuadtreeMesh& mesh,
        const FaceSegmentGraph& faceGraph,
        const Solvers::Poisson::CutCellMetrics* metrics,
        const PicDomainBC& domainBC,
        const PicEmbeddedBC* embeddedBC,
        const PicEmbeddedBCMap* embeddedBCById,
        const CellField<double>& phi,
        const Options& opt = defaultOptions());

    // Reconstruct cell-centered E_K by weighted LS.
    static void reconstructCellCenteredField(
        const QuadtreeMesh& mesh,
        const FaceSegmentGraph& faceGraph,
        const Solvers::Poisson::CutCellMetrics* metrics,
        const std::vector<double>& EnBySegment,
        const std::vector<double>* EnGammaByCellIndex,
        const PicEmbeddedBC* embeddedBC,
        const PicEmbeddedBCMap* embeddedBCById,
        CellField<Vector2D>& Ecell,
        const Options& opt = defaultOptions());

    // Convenience: face-normal -> (optional EB normal) -> cell-centered.
    static void evaluateCellCenteredField(
        const QuadtreeMesh& mesh,
        const FaceSegmentGraph& faceGraph,
        const Solvers::Poisson::CutCellMetrics* metrics,
        const PicDomainBC& domainBC,
        const PicEmbeddedBC* embeddedBC,
        const PicEmbeddedBCMap* embeddedBCById,
        const CellField<double>& phi,
        CellField<Vector2D>& Ecell,
        std::vector<double>* EnBySegmentOut = nullptr,
        std::vector<double>* EnGammaByCellIndexOut = nullptr,
        const Options& opt = defaultOptions());
};

} // namespace Field
} // namespace ExplicitPIC
