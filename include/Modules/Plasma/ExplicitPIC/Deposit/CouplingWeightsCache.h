#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "Core/Base/Vector2D.h"

class BackgroundMesh;

template<typename T>
class VertField;

// ----------------------------------------------------------------------------
// CouplingWeightsCache
//
// CSR storage for particle-to-cell coupling weights.
//
// This class builds per-particle coupling weights W_{pK} and stores them in CSR form,
// for charge-amount deposition and matched gather. Weights can be optionally cut-aware
// and renormalized to enforce sum_K W_{pK}=1 for each particle.
// ----------------------------------------------------------------------------

class CouplingWeightsCache
{
public:
    struct HpStats
    {
        double min = 0.0;
        double max = 0.0;
        double mean = 0.0;
        std::size_t count = 0;
    };

    // Build options for coupling weight construction.
    //
    // Notes:
    // - We index per-cell alpha by *cell index* (0..mesh.cellIndexNum()-1),
    //   consistent with CellField storage.
    // - Cut detection uses cell.isCut[levelSetI]. If the mesh has not been
    //   tagged for cuts, isCut is treated as false (alpha ignored).
    struct BuildOptions
    {
        int levelSetI = 0;

        // Periodic coupling support:
        // If enabled, particles whose CIC support overlaps the periodic boundary will
        // build coupling weights using periodic ghost images, then merge the entries
        // (by cell index). This makes deposition/gather consistent with periodic
        // particle motion (PB2) without having to replicate mesh topology.
        //
        // Notes:
        // - Periodicity is currently axis-aligned with the domain box.
        // - For quadtree/nonuniform meshes, correctness assumes the periodic faces are
        //   compatible (same refinement pattern). If not, weights remain well-defined
        //   but strict translational invariance is not guaranteed.
        bool periodicX = false;
        bool periodicY = false;

        // AMR coarse–fine composite band:
        // If enabled, when a particle's CIC support overlaps the coarse–fine neighborhood set B,
        // we recompute the stencil/weights using a reduced kernel width
        //   hp* = min_K h_K,
        // where K ranges over leaf cells intersecting the particle support.
        // This is a first-step AMR hardening against coarse–fine self-force artifacts.
        // Note: when enabled and the particle support overlaps the composite band,
        // hp is clamped to the minimum leaf-cell size in the stencil.
        bool enableCompositeBand = true;

        // Moment-exact cut-cell weighting:
        // If enabled and a level-set vertex field is provided, cut-cell weights
        // are computed by moment-exact geometric integration of the CIC kernel
        // over the fluid sub-region, instead of the crude alpha_K truncation.
        bool enableMomentExactCut = true;

        // Level-set samples at mesh vertices (typically iGeometry.levelSets[levelSetI]).
        // Required when enableMomentExactCut is true. If null, we fall back to alpha.
        const VertField<double>* levelSetVert = nullptr;

        // Sign convention for level set: insideIsNegative==true means solid is phi<0,
        // fluid is phi>=0 (default used throughout the EB Poisson solver).
        bool insideIsNegative = true;

        // Optional alpha_K by cell index. If null, alpha=1 is assumed.
        const std::vector<double>* alphaByCellIndex = nullptr;

        // If true and levelSetI is invalid, infer cut-cell status from alpha.
        bool cutByAlpha = false;
        double cutAlphaEps = 1e-12;

        // Optional host leaf cell labels for each particle (same length as positions).
        // If null, host cell is computed via mesh.findCellLabel(position).
        const std::vector<int>* hostCellLabels = nullptr;

        // Optional: print stencil summary for a specific particle.
        bool debugPrintStencil = false;
        std::size_t debugParticleId = std::numeric_limits<std::size_t>::max();

        // Debug: dump one particle's weights after build().
        // Effective only when debugParticleId is set.
        bool dumpWeightsOnce = true;
        std::string dumpWeightsPath; // e.g. "output/particle_weights.csv"

        // Performance controls
        // parallelBuild: enable OpenMP parallel construction (two-pass, prefix-sum CSR fill).
        // deterministicOrdering: use static scheduling for reproducible build across runs.
        // ompMinParticles: enable OpenMP only when numParticles >= this threshold.
        bool parallelBuild = true;
        bool deterministicOrdering = true;
        int ompMinParticles = 256;

        // Optional bookkeeping controls for hot paths.
        bool recordHpStats = true;
        bool recordPerParticleHp = true;
        bool recordPerParticleCompositeBand = true;

        // Zero-weight policy:
        // When a particle/sample support produces no valid overlap weights, the default
        // behavior is to fall back to the host cell with weight=1 to preserve a valid CSR.
        // Orbit-averaged quadrature samples can disable this so solid-side / clipped
        // samples contribute nothing instead of injecting a spurious host-cell weight.
        bool fallbackToHostOnZeroWeight = true;
    };

    // g++ requires avoiding `BuildOptions{}` as a default argument inside the enclosing class
    // when BuildOptions has default member initializers. Provide an explicit factory instead.
    static BuildOptions defaultBuildOptions() { return BuildOptions{}; }

    CouplingWeightsCache() = default;

    // Clear all cached entries.
    void clear();

    // Assign CSR weights directly.
    //
    // This is intended for advanced workflows (e.g., orbit-averaged overlaps)
    // where the caller builds/merges CSR entries externally and wants to reuse
    // the existing deposit/gather kernels without re-running the stencil build.
    //
    // Requirements:
    // - offsets.size() == numParticles + 1
    // - offsets is non-decreasing and offsets.back() == weights.size()
    // - cellIndices.size() == weights.size()
    void assignCSR(std::size_t numParticles,
                   std::vector<int> offsets,
                   std::vector<int> cellIndices,
                   std::vector<double> weights);

    // Compatibility alias for plan/docs wording: identical to assignCSR().
    void buildFromCSR(std::size_t numParticles,
                      std::vector<int> offsets,
                      std::vector<int> cellIndices,
                      std::vector<double> weights);

    // Reserve storage for an upcoming build.
    void reserve(std::size_t numParticles, std::size_t nnzGuess);

    // Placeholder build hook (kept for backwards compatibility).
    // Later tasks will replace this with the real CIC analytic integration + cut alpha
    // truncation + renormalization pipeline.
    void buildPlaceholder(std::size_t numParticles = 0);

    // Build coupling weights cache.
    //
    // Inputs:
    //   mesh     : BackgroundMesh / QuadtreeMesh (geometry + neighbor topology)
    //   positions: particle positions
    //   hp       : particle CIC half-widths. If hp.size()==1, it is treated as constant
    //              but may be overridden by the local host-cell size (paperA default).
    //   opt      : build options, including level-set index and optional alpha_K.
    //
    // Output (CSR):
    //   offsets_[p]..offsets_[p+1]-1 list entries for particle p.
    //   cellIndices_ store *cell index* (for direct CellField access).
    //   weights_ store renormalized tilde weights.
    //
    // Requirements enforced:
    //   - For each particle: sum(tildeW)=1 (with a small residual correction).
    //   - If sumW==0 for a particle: fallback to host cell with weight=1, and emit warning.
    void build(const BackgroundMesh& mesh,
               const std::vector<Vector2D>& positions,
               const std::vector<double>& hp,
               const BuildOptions& opt = defaultBuildOptions());

    // Accessors
    std::size_t numParticles() const { return numParticles_; }
    std::size_t nnz() const { return weights_.size(); }

    // CSR arrays:
    //   offsets[p]..offsets[p+1]-1 are entries for particle p.
    const std::vector<int>& offsets() const { return offsets_; }
    const std::vector<int>& cellIndices() const { return cellIndices_; }
    const std::vector<double>& weights() const { return weights_; }

    // Optional raw (pre-renorm) weights (debug only; off by default).
    void enableRawWeights(bool enabled);
    bool hasRawWeights() const { return keepRawWeights_; }
    const std::vector<double>& rawWeights() const { return rawWeights_; }

    // Last build's effective hp statistics (after any local clamping).
    bool hasHpStats() const { return hasHpStats_; }
    const HpStats& lastHpStats() const { return lastHpStats_; }

    // Per-particle effective hp used to build weights (size = numParticles).
    const std::vector<double>& effectiveHp() const { return effectiveHp_; }

    // Per-particle composite-band flag (size = numParticles, 0/1).
    const std::vector<unsigned char>& inCompositeBand() const { return inCompositeBand_; }

    // Last build's zero-weight fallback count (particles that fell back to host cell).
    // For assignCSR()/buildPlaceholder(), this is reset to 0.
    std::size_t lastBuildFallbackCount() const { return lastBuildFallbackCount_; }
    std::size_t lastBuildParticleCount() const { return lastBuildParticleCount_; }

    // Debug dump for a single particle.
    // Returns true if file written.
    bool dumpParticleWeightsCSV(const std::string& filePath, std::size_t particleId) const;

    // Detailed debug dump (cell label + center + alpha + Wraw + Wtilde).
    // Notes:
    // - alpha is taken from alphaByCellIndex (if provided) for cut cells; otherwise 1.
    // - Wraw is available only if raw weights were stored (enableRawWeights(true) or
    //   build() performed a debug dump).
    bool dumpParticleWeightsDetailedCSV(const BackgroundMesh& mesh,
                                        int levelSetI,
                                        const std::vector<double>* alphaByCellIndex,
                                        const std::string& filePath,
                                        std::size_t particleId) const;

private:
    std::size_t numParticles_ = 0;

    std::vector<int> offsets_;
    std::vector<int> cellIndices_;
    std::vector<double> weights_;

    bool keepRawWeights_ = false;
    std::vector<double> rawWeights_;

    // If dumpWeightsOnce=true, dump only once per run.
    bool dumpedWeightsOnce_ = false;

    HpStats lastHpStats_{};
    bool hasHpStats_ = false;

    std::vector<double> effectiveHp_;
    std::vector<unsigned char> inCompositeBand_;

    std::size_t lastBuildFallbackCount_ = 0;
    std::size_t lastBuildParticleCount_ = 0;
};
