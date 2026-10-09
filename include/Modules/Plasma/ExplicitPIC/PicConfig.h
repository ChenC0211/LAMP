#pragma once

#include "Modules/Plasma/Diagnostics/PicParticleTrajectoryConfig.h"

#include <string>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Modules/Plasma/Particles/ParticleManagementConfig.h"
#include "Solvers/Poisson/PoissonBC.h"

// ----------------------------------------------------------------------------
// PicConfig
//
// Configuration PODs for the Explicit PIC module.
// Parsed from ./system/pic.dict via DictionaryPIC.
//
// Keep this header lightweight (no heavy solver/mesh includes).
// ----------------------------------------------------------------------------

struct DepositionConfig
{
    std::string shape = "CIC"; // NGP/CIC/TSC (future)
    // Runtime defaults only; no longer overridden by pic.dict.
    bool cutRenorm = true;
    // Runtime defaults only; no longer overridden by pic.dict.
    bool compositeBand = true;
    // Runtime defaults only; no longer overridden by pic.dict.
    bool momentExactCut = false;
    bool ebRegularization = false;
    // Runtime defaults only; no longer overridden by pic.dict.
    bool useMatchedWeights = true;

    // ---------------------------------------------------------------------
    // Debug controls
    //
    // If debugParticleId >= 0, CouplingWeightsCache::build() can dump that
    // particle's coupling weights after building the cache:
    //   (cell, center, alpha, Wraw, Wtilde)
    //
    // dumpWeightsOnce: dump only on the first build() call in a run.
    // dumpWeightsPath: output CSV path.
    // ---------------------------------------------------------------------
    int debugParticleId = -1;
    bool dumpWeightsOnce = true;
    std::string dumpWeightsPath = ""; // e.g. "output/particle_weights.csv"
};

struct OutputConfig
{
    int every = 10;

    // -----------------------------------------------------------------
    // Field selection for case outputs
    //
    // Empty list means "use module defaults" (keeps backward-compat).
    // Suggested usage in pic.dict:
    //   Output { fields ( phi rho n Te E ) }
    // -----------------------------------------------------------------
    std::vector<std::string> fields;

    // Console verbosity controls for ExplicitPIC diagnostics.
    // level: error|warn|info|debug|trace (interpretation is module-defined).
    // every: print physical summary once every N steps at info-level.
    struct ConsoleConfig
    {
        std::string level = "info";
        int every = 10;

        // ImplicitPIC event logging controls (console anti-spam for dQrel).
        // dQrelRecover < 0 means "auto = 0.7 * dQrelWarn".
        double dQrelWarn = 5e-2;
        double dQrelRecover = -1.0;
        int dQrelDebounceSteps = 2;
        int dQrelEventCooldownSteps = 50;
        bool dQrelPrintRecover = true;

        // Physics-oriented summary control (primarily used by ImplicitPIC).
        // detail: minimal | medium | high
        bool physicsEnabled = true;
        std::string physicsDetail = "medium";
        bool physicsShowCharge = true;
        bool physicsShowCurrent = true;
        bool physicsShowWall = true;
        bool physicsShowMCC = true;
        bool physicsShowPerf = true;
    };
    ConsoleConfig console;

    // -----------------------------------------------------------------
    // Structured series cadence used by ImplicitPIC/ExplicitPIC secondary
    // outputs. Public dict schema exposes this as Output.series.every.
    // -----------------------------------------------------------------
    struct V2Config
    {
        bool enabled = false;
        int seriesEvery = 0;
        int fieldsEvery = 500;
        int perfEvery = 100;
        int physicsEvery = 100;
        int speciesEvery = 500;
        int amountEvery = 0;
        int boundaryEvery = 500;
        int flushEvery = 100;

        struct ConsoleV2Config
        {
            // oneLine | multiLine
            std::string format = "oneLine";
        } console;
    } v2;

    // -----------------------------------------------------------------
    // Periodic diagnostics (optional)
    //
    // When periodic wrapping is enabled (x/y), it is useful to monitor
    // (i) how often particles wrap and (ii) whether the solved potential
    // is periodic to numerical tolerance.
    //
    // This block is intentionally lightweight: it prints a short summary
    // to the console and optionally appends a CSV line.
    // -----------------------------------------------------------------
    struct PeriodicDiagnosticsConfig
    {
        bool enabled = false;
        bool printConsole = true;
    } periodicDiagnostics;

    // -----------------------------------------------------------------
    // Current waveform diagnostics (optional)
    //
    // Writes a time series CSV for the driven electrode:
    //   step, t, V_drive, J_cond, J_disp, J_total
    //
    // side: left|right|bottom|top|auto
    //   auto: prefer the first Dirichlet side that has timeDirichlet.enabled;
    //         otherwise use the first Dirichlet side.
    // Cadence is grouped under Output.series.every.
    // -----------------------------------------------------------------
    struct CurrentWaveformConfig
    {
        bool enabled = false;
        std::string side = "auto";
    } currentWaveform;

    // -----------------------------------------------------------------
    // 1D profile diagnostics (optional)
    //
    // Produces a binned profile along x by volume-averaging over y.
    // This is intended to support 2D-with-periodic-y emulation of common
    // 1D CCP/RF benchmarks (e.g. Turner et al. 2013).
    //
    // The output is written as a long CSV table: one row per (step,bin).
    // -----------------------------------------------------------------
    struct Profile1DConfig
    {
        bool enabled = false;

        // Number of bins in x.
        int nbinsX = 256;
    } profile1D;

    // -----------------------------------------------------------------
    // Time-averaged 1D profile diagnostics (optional)
    //
    // This is intended for benchmark comparisons where a steady periodic
    // state is characterized by phase-averaged or period-averaged profiles.
    // The solver accumulates statistics over a trailing window of
    // avgWindowSteps and writes a single CSV row per x-bin.
    //
    // avgWindowSteps: number of steps to average over (<=0 disables).
    // Cadence is grouped under Output.fields.every.
    // -----------------------------------------------------------------
    struct Profile1DAvgConfig
    {
        bool enabled = false;
        int avgWindowSteps = 0;
    } profile1DAvg;

    // -----------------------------------------------------------------
    // Source/sink current-balance diagnostics (optional)
    //
    // Writes a CSV time series that tracks injected source current and
    // absorbed/sink currents (domain outflow + EB wall absorption).
    // Cadence is grouped under Output.series.every.
    // -----------------------------------------------------------------
    struct CurrentBalanceConfig
    {
        bool enabled = false;
    } currentBalance;

    struct AmountConfig
    {
        bool enabled = false;
    } amount;

    // -----------------------------------------------------------------
    // EB near-wall E-field RMS proxy (optional)
    //
    // Tracks RMS(|E|) in a thin EB-adjacent cell band.
    // bandLayers:
    //   0 -> cut cells only
    //   1 -> cut cells + 1-ring neighbors
    //   ...
    // Cadence is grouped under Output.series.every.
    // -----------------------------------------------------------------
    struct EBBandRMSConfig
    {
        bool enabled = false;
        int bandLayers = 1;
    } ebBandRMS;

    // Unified physics diagnostics written under output/<baseDir>/.
    struct DiagnosticsConfig
    {
        struct ElectrodeGroup
        {
            std::string name;
            std::vector<int> boundaryIds;
        };

        struct EnergySample
        {
            std::string name;
            std::string species = "all";
            std::string location = "bulk";
        };

        bool enabled = false;
        int every = 0;
        std::string baseDir = "diagnostics";

        struct ElectrodeCurrentConfig
        {
            bool enabled = false;
            int every = 0;
            bool perBoundary = true;
            bool includeDomain = false;
            std::vector<ElectrodeGroup> groups;
        } electrodeCurrent;

        struct ReactionStatsConfig
        {
            bool enabled = false;
            int every = 0;
            bool summary = true;
            bool perChannel = true;
            bool perProjectile = true;
            bool spawnedSpecies = true;
        } reactions;

        struct EnergyHistogramConfig
        {
            bool enabled = false;
            int every = 0;
            int bins = 80;
            double maxEnergyEV = 200.0;
            std::vector<EnergySample> samples;
        } energyHistogram;
    } diagnostics;

    using ParticleTrajectoryConfig = PicParticleTrajectoryConfig;
    ParticleTrajectoryConfig particleTrajectory;

};

// Checkpoint / restart controls for ExplicitPIC.
// Stores particle caches under checkpoint.path.
struct CheckpointConfig
{
    bool enabled = false;
    int every = 100;              // write every N steps (0 disables)
    int keep = 2;                 // keep last K checkpoints (<=0 keeps all)
    std::string path = "../cache"; // default: cases/cache when run from case dir

    bool restart = false;         // if true, resume from checkpoint on start
    int restartStep = -1;         // <0: latest; >=0: use that step directory

    bool saveRng = true;          // persist RNG state for reproducibility
};

// Adaptive mesh refinement (AMR) controls for ExplicitPIC.
// These are PIC-side policies; mesh core logic is unchanged.
struct AmrConfig
{
    bool enabled = false;
    int adaptEverySteps = 20;
    // How often to recompute expensive AMR indicators/targets (in steps).
    // 1 = every step (default, legacy behavior).
    int indicatorEverySteps = 1;
    int persistenceSteps = 3;
    int cellBudgetPerAdapt = 0; // 0 = no budget
    double hysteresis = 0.2;    // ratio-based hysteresis band

    struct DebyeConfig
    {
        bool enabled = true;
        double CD = 0.5;        // h <= CD * lambda_D
        double Te_eV = -1.0;    // <=0: fallback to electron species Te
        double neFloor = 1e12;  // 1/m^3
        int avgWindowSteps = 20;
        // Spatial smoothing radius (in leaf-cell widths). 0 disables.
        int smoothWindowCells = 0;
    };
    DebyeConfig debye;

    struct IndicatorConfig
    {
        bool useEField = true;
        double eFieldPercentile = 95.0;
        int eFieldAvgWindow = 20;
        // Spatial smoothing radius (in leaf-cell widths). 0 disables.
        int eFieldSmoothWindowCells = 0;
        // Absolute floor for E (SI, V/m). 0 disables.
        double eFieldAbsFloor = 0.0;
        // Clamp reference to a fraction of Emax (0 disables).
        double eFieldRefFrac = 0.0;
        // EMA window for Eref (>=1, 1 disables smoothing).
        int eFieldRefAvgWindow = 1;
        // NOTE: EB band refinement has been removed. Use mesh.dict GeometryRefinement instead.
    };
    IndicatorConfig indicator;

    struct ParticleRelabelConfig
    {
        std::string mode = "localMapThenFindCell";
        std::string fallback = "findCellLabel";
    };
    ParticleRelabelConfig particleRelabel;

    struct GeometryRefreshConfig
    {
        bool present = false;
        // full | adaptive | skip
        std::string mode = "full";
        double changeRel = 0.2;
        bool verbose = false;
        bool allowWithMovingEb = false;
        // full | exact_refresh (experimental, env-gated in ImplicitPIC)
        std::string movingEbRefreshMode = "full";
    };
    GeometryRefreshConfig geometryRefresh;

    struct ReportConfig
    {
        // Runtime defaults only; no longer overridden by pic.dict.
        bool enabled = true;
        int every = 20;
        std::string path = "./output/series/amr_report.csv";
    };
    ReportConfig report;
};

// Secondary electron emission (SEE) configuration.
//
// Minimal model: ions that impact a chosen cathode-side boundary emit
// gamma electrons per incident ion (physical), with a prescribed emission
// energy. The emission angular distribution is added in later steps.
struct SEEConfig
{
    bool enabled = false;
    // Emitting sides (left/right/bottom/top). If empty, falls back to cathodeSide.
    std::vector<std::string> cathodeSides;
    std::string cathodeSide = "left";        // legacy single-side key
    std::string emittedSpecies = "electron"; // species name to inject
    // -----------------------------------------------------------------
    // Emission yield model (gamma)
    //
    // Backward compatible:
    //   - If gammaModel is unset, we default to "constant" and use `gamma`.
    //
    // Supported models:
    //   - constant: gamma(E) = gamma
    //   - powerlaw: gamma(E) = min(gammaMax, gamma0*(max(E-Emin,0)/Eref)^alpha)
    //   - piecewise: linear interpolation of (E_i, gamma_i) points
    // -----------------------------------------------------------------
    std::string gammaModel = "constant";   // constant | powerlaw | piecewise

    // Constant model (legacy)
    double gamma = 0.0;                    // electrons per ion (physical)

    // Power-law model parameters (E in eV)
    double gamma0 = 0.0;
    double Eref_eV = 100.0;
    double alpha = 0.0;
    double gammaMax = 10.0;
    double Emin_eV = 0.0;

    // Piecewise-linear table (E in eV). Must contain >= 2 points.
    std::vector<double> piecewiseE_eV;
    std::vector<double> piecewiseGamma;

    // -----------------------------------------------------------------
    // Emitted electron energy model
    //
    // Backward compatible:
    //   - Default is "mono" and uses emittedEnergy_eV (legacy key: emittedEnergy).
    //
    // Supported models:
    //   - mono       : E = emittedEnergy_eV
    //   - maxwell    : E ~ Maxwell-Boltzmann (3D) with kT = emittedTe_eV
    //   - exponential: E ~ Exp(mean = emittedEmean_eV)
    // -----------------------------------------------------------------
    std::string emittedEnergyModel = "mono"; // mono | maxwell | exponential
    double emittedEnergy_eV = 2.0;           // mono energy [eV] (legacy)
    double emittedTe_eV = 2.0;               // maxwell kT [eV]
    double emittedEmean_eV = 2.0;            // exponential mean [eV]
    double injectEps_m = 1e-6;             // inward offset from wall [m]
    int maxEmitMacro = 100;                // cap macro particles per impact

    // Angular distribution of emitted electrons:
    //   - "cosine"  : 3D cosine-weighted hemisphere about inward normal (default)
    //   - "isotropic": 3D isotropic hemisphere about inward normal
    //   - "fixed"   : emit along inward normal only
    std::string angular = "cosine";
};

struct SpeciesConfig
{
    std::string name;
    // Physical charge [C] and mass [kg].
    // These are used for MCC and to compute the normalized push coefficient kappa_s.
    double q = 0.0;
    double m = 1.0;
    double weight = 1.0; // macro-particle weight w_p
    int numParticles = 0;
    double energyCutoff_eV = 0.0;

    // Initialization knobs (minimal set; can be extended later)
    std::string init = "uniform";
    // Species thermal energy parameter Te [eV].
    // Interpreted as an isotropic Maxwellian with kT = Te*e.
    double Te_eV = 0.0;

    // Drift velocity [m/s] (physical). Will be mapped to code by v0=L0/T0.
    Vector2D drift = Vector2D(0.0, 0.0);

    // Optional source emission (boundary injection).
    struct SourceConfig
    {
        bool enabled = false;
        int boundaryId = 0;               // geometry boundary id (1-based)
        double numberDensity_m3 = 0.0;    // source number density [1/m^3]
        double temperature_eV = -1.0;     // if <0, use SpeciesConfig.Te_eV
        Vector2D drift = Vector2D(0.0, 0.0);
        bool driftSet = false;
        double current_A = 0.0;           // optional current override [A]
        int every = 1;                    // emit cadence (steps)
        int maxEmitMacro = 0;             // 0 -> no cap
    };
    SourceConfig source;
};

struct PicConfig
{
    double dt = 1e-8;
    int steps = 10;
    unsigned int randomSeed = 0;
    // If >0, override OpenMP thread count for the ExplicitPIC module.
    int ompThreads = 0;

    // -----------------------------------------------------------------
    // ImplicitPIC performance knobs (safe defaults)
    //
    // segmenterMaxLevel:
    //   -1  : default/baseline behavior (segment predictor orbits on leaf cells)
    //   >=0 : segment on cells coarsened up to this maximum level to reduce
    //         segments/samples per particle after AMR.
    //
    // Notes:
    // - Currently used by ImplicitPIC predictor orbit segmenter.
    // - If embedded-boundary interception splitting is enabled, the solver
    //   ignores this and uses leaf segmentation.
    // -----------------------------------------------------------------
    int segmenterMaxLevel = -1;
    // Runtime defaults only; no longer overridden by pic.dict.
    int segmenterMaxSegmentsPerParticle = 8192;

    struct PolarizationPoissonConfig
    {
        // Runtime defaults only; no longer overridden by pic.dict.
        // auto | direct | ic | cg
        // auto: prefer iterative IC-CG in ImplicitPIC and fall back to direct.
        std::string mode = "direct";
        // When true, keep using the cached direct LU factors/backsolves until
        // the matrix pattern changes (typically AMR / topology change).
        bool freezeDirectMatrixUntilPatternChange = true;
        bool temporalGuess = true;
        double temporalOmega = 0.8;
        bool reuseIterativePattern = true;
    } polarizationPoisson;

    // -----------------------------------------------------------------
    // ImplicitPIC orbit-averaged weights knobs
    // -----------------------------------------------------------------

    // Runtime defaults only; no longer overridden by pic.dict.
    // If true, enable the fused orbit-avg weights build path.
    // (May still be overridden via env var QUANTUM_IMPLICIT_ORBITAVG_FUSED.)
    bool orbitAvgFused = false;

    struct OrbitAvgAdaptiveConfig
    {
        // Runtime defaults only; no longer overridden by pic.dict.
        // If true, build orbit-averaged overlaps by adaptive time quadrature
        // over the whole predictor interval [t_n, t_{n+1}] (instead of segmenting
        // on cell-boundary events with fixed Gauss2 per segment).
        bool enabled = true;

        // Absolute tolerance for the local error estimate based on |I4 - I2|.
        // (This is a numerical quadrature control; physical acceptance is still
        // validated via dQrel/rG in the case compare harness.)
        double tol = 5e-3;

        // Maximum bisection depth for adaptive quadrature.
        int maxRefineDepth = 3;

        // If true, use a cheap refinement criterion based on host-cell label
        // changes along the interval (instead of the |I4 - I2| estimator).
        // This is intended as a low-cost "补刀" path and is default-on.
        bool cheapRefine = true;

        // If true, when cheapRefine hits maxRefineDepth on a crossing interval,
        // try one extra label-only bisection split before accepting Gauss2.
        bool crossSplitAtMaxDepth = false;

        // Maximum bisection iterations for locating the crossing split point t*.
        int crossSplitMaxIters = 10;

        // Sheath gate pad (meters): enable cross-split only for intervals close
        // to walls/cut-cells to keep overhead bounded.
        double sheathPad_m = 2e-4;

        // Per-particle cap for extra cross-splits at max depth.
        int crossSplitMaxSplits = 1;

        // Quadrature rule selector (currently: gauss2_4).
        std::string rule = "gauss2_4";
    } orbitAvgAdaptive;

    // -----------------------------------------------------------------
    // ImplicitPIC orbit-avg hybrid: coarse segmentation + adaptive patching
    // -----------------------------------------------------------------
    struct OrbitAvgPatchConfig
    {
        // If true, build orbit-averaged overlaps on the (possibly coarsened)
        // predictor segments, using fixed Gauss2 on most segments and applying
        // adaptive refinement only on segments that appear to cross leaf cells.
        //
        // Intended usage: combine with segmenterMaxLevel >= 0 to reduce
        // predictor segments after AMR, while keeping accuracy via "补刀".
        bool enabled = false;

        // Absolute tolerance for patch refinement error estimate |I4 - I2|.
        double tol = 1e-3;

        // Maximum bisection depth for patch refinement.
        int maxRefineDepth = 8;
    } orbitAvgPatch;

    // Wall interaction model (currently SEE only).
    SEEConfig see;

    std::vector<SpeciesConfig> species;

    PicDomainBC domainBC;
    PicEmbeddedBC embeddedBC;

    // V3-style boundaryId BCs (parsed from BC/Field and BC/Particle).
    struct BoundaryIdFieldBC
    {
        int id = -1;
        PicBoundarySide bc;
    };
    struct BoundaryIdParticleBC
    {
        int id = -1;
        std::string model = "absorb";
        SEEConfig see; // optional SEE override for this boundary
        bool seeSpecified = false;
    };
    std::vector<BoundaryIdFieldBC> fieldBCById;
    std::vector<BoundaryIdParticleBC> particleBCById;

    // (Note) SEEConfig is declared above; keep a single instance here.

    DepositionConfig deposition;
    OutputConfig output;

    // Adaptive mesh refinement controls (PIC-side policy).
    AmrConfig amr;

    // Particle population control (merge/split/thinning).
    ParticleManagementConfig particleManagement;

    // Checkpoint / restart
    CheckpointConfig checkpoint;
};
