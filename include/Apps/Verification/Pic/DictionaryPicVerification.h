#pragma once
#include "Core/IO/Dictionary.h"
#include "Modules/Plasma/ExplicitPIC/IO/DictionaryPIC.h"

// Settings consumed only by PicVerification, never by a PIC time integrator.
struct PicVerificationConfig
{
    // -----------------------------------------------------------------
    // Regression harness (app-level baseline record/check)
    //
    // Used by verification / regression Apps that compute a set of
    // scalar metrics and compare them against a saved baseline.
    //
    // mode:
    //   - "off"   : do not read/write baseline, always pass
    //   - "record": write baseline file and exit success
    //   - "check" : compare against baseline file; non-zero exit on mismatch
    //
    // baseline is stored under the case directory (recommended: ./cache/).
    // -----------------------------------------------------------------
    struct RegressionHarnessConfig
    {
        std::string mode = "off"; // off | record | check
        std::string baselinePath = "./cache/regression_baseline.dict";
        std::string reportPath = "./output/meta/regression_report.txt";
        double absTol = 1e-12;
        double relTol = 1e-10;
    };

    RegressionHarnessConfig regressionHarness;

    // Regression tests
    // Verification/regression (handled by selftest tools). Default off for released cases.
    bool regressionUnitBoxSinSin = false;

    // Same unit-box regression but with y-periodic wrapping (bottom/top Periodic).
    // This specifically exercises Periodic BC assembly and is useful for validating
    // periodic closure before enabling it in RF PIC cases.
    bool regressionUnitBoxSinCosPeriodicY = false;

    // Regression: ElectricFieldEvaluator on an analytic field with y-periodic wrap.
    // Uses phi = sin(pi x) cos(2 pi y) on a unit box with x-Dirichlet (phi=0) and
    // y-Periodic boundaries. This compares periodic vs non-periodic boundary handling
    // to catch boundary spikes in E_n.
    bool regressionElectricFieldSinCosPeriodicY = false;

    // Regression: cut-circle cut-cell metrics.
    struct RegressionCutCircleConfig
    {
        bool enabled = false;
        int levelSetIndex = 0;
        Vector2D center = Vector2D(0.0, 0.0);
        double radius = 0.25;
        bool dumpCSV = true;
        std::string cellsCSV = "./output/fields/regression_cut_circle_cells.csv";
        std::string segmentsCSV = "./output/fields/regression_cut_circle_segments.csv";

        // If true, treat inside as (psi <= 0). Otherwise inside is (psi >= 0).
        bool insideIsNegative = true;

        // Alpha sampling in cut cells (N x N grid)
        int alphaSamplingN = 4;

        // Clamp to avoid exactly 0/len cut-length on nearly-tangent edges.
        double cutLengthClampEps = 1e-12;
    };

    RegressionCutCircleConfig regressionCutCircle;

    // Regression: FV Poisson on a cut-circle domain using a
    // manufactured solution (q^2 * sin(pi X) sin(pi Y)).
    //
    // The outer box uses homogeneous Dirichlet (phi=0). The embedded boundary is
    // the circle given by (center,radius), with cut-cell metrics computed from
    // the analytic circle level-set.
    struct RegressionCutCirclePoissonConfig
    {
        bool enabled = false;

        // Geometry / level-set
        int levelSetIndex = 0;
        Vector2D center = Vector2D(0.0, 0.0);
        double radius = 0.25;
        bool insideIsNegative = true;
        int alphaSamplingN = 6;
        double cutLengthClampEps = 1e-12;

        // Manufactured solution scaling (phi = coeff * q^2 * sin(pi X) sin(pi Y)).
        double coeff = 1.0;

        // Solver controls are configured in ./system/solvers.dict (PoissonCellCentered).

        // If true and BC/embedded.enabled=true, pass embedded BC to Poisson solver.
        bool useEmbeddedBC = false;

        // If true, compute E-field error norms against the analytic gradient.
        bool computeEError = false;

        // Optional dumps
        bool dumpPhiCSV = false;
        std::string phiCSV = "./output/fields/regression_cut_circle_phi.csv";
        bool dumpMetricsCSV = true;
        std::string cellsCSV = "./output/fields/regression_cut_circle_poisson_cells.csv";
        std::string segmentsCSV = "./output/fields/regression_cut_circle_poisson_segments.csv";
    };

    RegressionCutCirclePoissonConfig regressionCutCirclePoisson;

    // Regression: FV polarization Poisson on a cut-circle domain using a
    // manufactured solution (q^2 * sin(pi X) sin(pi Y)) with either constant
    // or variable polarization coefficient kappa(x,y).
    //
    // Governing form:
    //   -div(kappa * grad(phi)) = rhs
    //
    // For constant kappa, rhs is simply scaled from the kappa=1 MMS.
    struct RegressionCutCirclePolarizationPoissonConfig
    {
        bool enabled = false;

        // Geometry / level-set
        int levelSetIndex = 0;
        Vector2D center = Vector2D(0.0, 0.0);
        double radius = 0.25;
        bool insideIsNegative = true;
        int alphaSamplingN = 6;
        double cutLengthClampEps = 1e-12;

        // Manufactured solution scaling (phi = coeff * q^2 * sin(pi X) sin(pi Y)).
        double coeff = 1.0;

        // Constant coefficient kappa in -div(kappa grad(phi)). Must be >= 0.
        double kappa = 1.0;
        // kappa mode:
        //   constant: use scalar kappa above
        //   sin2d:    kappa(x,y)=kappa0*(1+kappaAmp*sin(kappaKx*x)*sin(kappaKy*y))
        std::string kappaMode = "constant";
        double kappa0 = 1.0;
        double kappaAmp = 0.0;
        double kappaKx = 0.0;
        double kappaKy = 0.0;

        // If true and BC/embedded.enabled=true, pass embedded BC to solver.
        bool useEmbeddedBC = false;

        // If true, compute E-field error norms against the analytic gradient.
        bool computeEError = false;

        // Optional dumps
        bool dumpPhiCSV = false;
        std::string phiCSV = "./output/fields/regression_cut_circle_pol_poisson_phi.csv";
        bool dumpMetricsCSV = true;
        std::string cellsCSV = "./output/fields/regression_cut_circle_pol_poisson_cells.csv";
        std::string segmentsCSV = "./output/fields/regression_cut_circle_pol_poisson_segments.csv";
    };

    RegressionCutCirclePolarizationPoissonConfig regressionCutCirclePolarizationPoisson;

    // Regression: embedded conductor (Dirichlet) response with a localized charge.
    struct RegressionEmbeddedConductorConfig
    {
        bool enabled = false;

        // Geometry (analytic circle level-set).
        int levelSetIndex = 0;
        Vector2D center = Vector2D(0.0, 0.0);
        double radius = 0.25;
        bool insideIsNegative = true; // solid is psi < 0
        int alphaSamplingN = 6;
        double cutLengthClampEps = 1e-12;

        // Source: single localized charge (Gaussian if sigma > 0).
        Vector2D particle = Vector2D(0.35, 0.40);
        double charge = 1.0;
        double sigma = 0.0; // 0 -> deposit to host cell only
        // chargeUnit:
        //   code: charge is used directly in code units for rhs assembly.
        //   si:   reserved for future conversion path; currently treated as code.
        std::string chargeUnit = "code";
        // sourceShape:
        //   delta:    point-like source (matched cache deposition in matched mode).
        //   gaussian: Gaussian source around particle location.
        std::string sourceShape = "gaussian";
        // sampleMode:
        //   hostCell: legacy host-cell-center sampling.
        //   matched:  matched gather with coupling weights.
        std::string sampleMode = "matched";
        // sigma auto policy (for gaussian source shape).
        bool sigmaAutoFromFine = true;
        double sigmaFineFactor = 1.5;

        // Output files.
        bool dumpPhiCSV = true;
        std::string phiCSV = "./output/fields/test5_phi.csv";
        bool dumpEnCSV = true;
        std::string enCSV = "./output/fields/test5_en.csv";
        bool dumpMetricsCSV = true;
        std::string cellsCSV = "./output/fields/test5_cells.csv";
        std::string segmentsCSV = "./output/fields/test5_segments.csv";
    };

    RegressionEmbeddedConductorConfig regressionEmbeddedConductor;

    // -----------------------------------------------------------------
    // paperA Test 3: Single-particle numerical self-force (periodic neutrality)
    // -----------------------------------------------------------------
    struct SelfForceTest3Config
    {
        bool enabled = false;

        // Which subtests to run (empty -> all):
        //   uniform, traverse, traverseCut, cutCross, periodicCross
        std::vector<std::string> subtests;

        struct ParticleConfig
        {
            double q = 1.0; // particle charge [C]
            double m = 1.0; // particle mass [kg] (unused)
        };
        ParticleConfig particle;

        struct NeutralityConfig
        {
            bool enabled = true;          // explicit neutrality
            std::string mode = "explicit";
        };
        NeutralityConfig neutrality;

        struct GaugeFixConfig
        {
            bool enabled = true;          // mean-gauge fix
            std::string mode = "mean";
        };
        GaugeFixConfig gaugeFix;

        struct UniformConfig
        {
            std::vector<int> Nlist; // override BackgroundMesh divide for uniform scaling
            std::vector<Vector2D> samplePoints = { Vector2D(0.37, 0.41) };
        };
        UniformConfig uniform;

        struct PathConfig
        {
            Vector2D start = Vector2D(0.25, 0.5);
            Vector2D end   = Vector2D(0.75, 0.5);
            int numSamples = 401;
        };
        PathConfig traverse;
        PathConfig traverseCut;
        PathConfig cutCross;

        struct PeriodicCrossConfig
        {
            Vector2D start = Vector2D(0.37, 0.98);
            Vector2D end   = Vector2D(0.37, 1.02);
            int numSamples = 401;
        };
        PeriodicCrossConfig periodicCross;

        // Ablation list (A0/A1). Empty -> all.
        std::vector<std::string> ablations;

        // Output files
        std::string outputUniform = "./output/fields/self_force_uniform.csv";
        std::string outputTraverse = "./output/fields/self_force_traverse.csv";
        std::string outputPeriodicCross = "./output/fields/self_force_periodic_cross.csv";
        std::string outputSummary = "./output/meta/self_force_summary.json";
    };

    SelfForceTest3Config selfForceTest3;

};

class DictionaryPicVerification : private Dictionary
{
public:
    bool readDict(const std::string& path);
    const PicConfig& physics() const { return pic_.config(); }
    const PicVerificationConfig& config() const { return config_; }
private:
    DictionaryPIC pic_;
    PicVerificationConfig config_;
    bool readRecord(const std::string& b, const std::string& parent,
                    const std::string& line, const std::vector<std::string>& tokens);
};
