#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "Core/Field/CellField.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Solvers/Poisson/PoissonBC.h"

class QuadtreeMesh;
class ParticleType;
class NormalizationContext;

namespace Solvers { namespace Poisson { class CutCellMetrics; } }

namespace ExplicitPIC
{
namespace Diagnostics
{

// PicDiagnostics
// ---------------------------------------------------------------------------
// Computes discrete Gauss residual r_{G,K} for the cell-centered FV Poisson
// operator (including optional EB Robin/Dirichlet/Neumann closure), plus a few
// global scalar diagnostics.
//
// Conventions:
//   Operator is assembled for:  -Delta(phi) = rhs
//   FV system is (A * phi) = b, where b includes:
//     - cell-integrated rhs: rhs_K * (alpha_K * area_K)
//     - domain BC contributions
//     - EB Robin closure contributions (if enabled)
//
// The Gauss residual is reported in the *assembled* balance form:
//   r_G = A*phi - b
// which is consistent with the solver's discrete FV balance (includes BC/EB).
// ---------------------------------------------------------------------------
class PicDiagnostics
{
public:
    struct Result
    {
        // Active unknowns count (FaceSegmentGraph-based)
        int nActive = 0;

        // Gauss residual norms over active unknown set
        double rG_L2 = 0.0;      // sqrt(sum r^2)
        double rG_Linf = 0.0;    // max |r|
        double rG_L2_per_active = 0.0; // sqrt(sum r^2 / nActive)
        double rG_L2_rel_rhs = 0.0;    // sqrt(sum r^2) / sqrt(sum b^2)
        double rG_Linf_rel_rhs = 0.0;  // max |r| / max |b|

        // Gauss residual spike indicators (max |r_G,c| over selected bands)
        //  - S_cf: coarse-fine composite band (refinement interfaces)
        //  - S_eb: embedded-boundary cut-cell band
        double S_cf = 0.0;
        double S_eb = 0.0;
        double S_eb_rel_rhs = 0.0;     // max |r| on cut-cell band / max |b|

        // Coarse-fine / embedded-boundary band quantiles of |r_G,c|.
        // These are more robust than a pure max and are useful for tracking
        // band-wide residual behavior across AMR events.
        double S_cf_p50 = 0.0;
        double S_cf_p90 = 0.0;
        double S_cf_p99 = 0.0;
        double S_eb_p50 = 0.0;
        double S_eb_p90 = 0.0;
        double S_eb_p99 = 0.0;

        // Electrostatic energy (discrete): W = 0.5 * sum_K Q_K * phi_K
        double energyW = 0.0;

        // Global charge totals
        double Qgrid = 0.0;      // sum_K Q_K
        double Qmacro = 0.0;     // sum_s (q_s * w_s * N_s)
        double dQ = 0.0;         // Qgrid - Qmacro
        double dQrel = 0.0;      // (Qgrid-Qmacro)/Qmacro (if |Qmacro|>0)

        // Optional physical summaries for cases/diagnostics.
        // Mean number densities are volume-averaged over the whole domain.
        std::vector<std::string> speciesNames;     // in config order
        std::vector<double> speciesNmean_1pm3;     // same size as speciesNames
        std::vector<long long> speciesMacroCount; // macro particle count per species
        std::vector<double> speciesPhysCount;     // physical count sum(w) per species
        long long macro_total = 0;
        double phys_total = 0.0;

        // Electron kinetic summaries (particle-weighted).
        // Emean is mean kinetic energy per electron [eV] based on 3D velocity-space.
        double electron_mean_energy_eV = 0.0;
        // Te is derived from Emean by Te = (2/3) * Emean for 3D Maxwellian convention.
        double Te_e_mean_eV = 0.0;
    };

    struct ExplicitSample
    {
        Result values;
        double phiMin, phiMax; // Code units; original all-cell scan, including NaN handling.
    };
    // Read-only explicit field-stage sample: Gauss balance, active-volume species
    // densities and 3D electron energy. Name matching precedes negative charge;
    // particle/default weight fallback and serial reduction order are preserved.
    // All fields and optional metrics must describe the current mesh.
    static ExplicitSample sampleExplicit(const QuadtreeMesh& mesh,
        const FaceSegmentGraph& graph, const PicDomainBC& domainBC,
        const CellField<double>& rhs, const CellField<double>& phi, const CellField<double>& Q,
        double totalMacroCharge, const std::vector<ParticleType>& species,
        const std::vector<double>& macroWeight, const NormalizationContext& normalization,
        const Solvers::Poisson::CutCellMetrics* cutMetrics,
        const PicEmbeddedBC& embeddedBC, const PicEmbeddedBCMap& embeddedBCById);

    // Current implicit diagnostic layer: retain variable coefficients and the
    // legacy all-four-children volume test. Population/energy arithmetic is
    // shared with the explicit sampler; no field extrema scan or observer here.
    static Result sampleImplicit(const QuadtreeMesh& mesh,
        const FaceSegmentGraph& graph, const PicDomainBC& domainBC,
        const CellField<double>& rhs, const CellField<double>& phi, const CellField<double>& Q,
        double totalMacroCharge, const std::vector<ParticleType>& species,
        const std::vector<double>& macroWeight, const NormalizationContext& normalization,
        const Solvers::Poisson::CutCellMetrics* cutMetrics,
        const PicEmbeddedBC* embeddedBC, const PicEmbeddedBCMap* embeddedBCById,
        const std::vector<double>& kappaBySegment, const std::vector<double>& kappaByCellIndex);

    // Compute Gauss residual and scalar diagnostics.
    //
    // Inputs:
    //   - mesh/graph: quadtree and oriented leaf-face segments
    //   - domainBC  : domain boundary conditions
    //   - rhs/phi   : cell-centered RHS and potential
    //   - Q         : cell-centered charge amount (C). If null, energy uses
    //                rho = rhs*epsilon0 and integrates rho*phi*area.
    //   - cutMetrics/embeddedBC: optional EB metrics + BC for Robin closure.
    //
    // Output:
    //   - rGout: optional per-cell residual field (size mesh.cellIndexNum())
    static Result compute(const QuadtreeMesh& mesh,
                          const FaceSegmentGraph& graph,
                          const PicDomainBC& domainBC,
                          const CellField<double>& rhs,
                          const CellField<double>& phi,
                          const CellField<double>* Q,
                          double Qmacro,
                          CellField<double>* rGout = nullptr,
                          const Solvers::Poisson::CutCellMetrics* cutMetrics = nullptr,
                          const PicEmbeddedBC* embeddedBC = nullptr,
                          const PicEmbeddedBCMap* embeddedBCById = nullptr,
                          const std::vector<double>* kappaBySegment = nullptr,
                          const std::vector<double>* kappaByCellIndex = nullptr);


};

} // namespace Diagnostics
} // namespace ExplicitPIC
