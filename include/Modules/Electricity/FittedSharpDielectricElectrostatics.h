#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>
#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Modules/Electricity/FittedElectricPotential.h"
#include "Core/Geometry/InterfaceTraceGeometry.h"
#include "Core/Mesh/UniformLeafGrid.h"

class LevelSetInterfacePanels;
namespace LevelSetHermiteGeometry { class Snapshot; }

// Electric result on an explicitly face-fitted active Cartesian leaf domain.
// This entry reuses the sharp P3 jump/phase-current assembler and SuperLU. It
// does not supply momentum loads or advance charge. All quantities use the
// caller's consistent code units, including charge and conductivity.
// Potential includes unchanged FINITE inactive payload from the supplied field;
// displacementFlux is the full corrected outward D integral in graph order.
struct FittedSharpDielectricResult
{
    FittedElectricPotential potential;
    std::vector<double> epsilonBySegment;
    std::vector<PolarizationPoissonFVSolverCellCentered::SegmentFluxCorrection> corrections;
    std::vector<double> displacementFlux;
    std::vector<double> conductionCurrentBySegment; // Requested J, integrated outward, graph order.
    // Only net boundary accounting is meaningful here. The legacy rectangular
    // lower/upper groups in report are not supplied; callers sum the segment
    // arrays using their actual boundary IDs instead of interpreting zeros.
    DielectricReport report;
    // Actual linear-solver work, including exact-matrix factor reuse. This is
    // observational; it does not replace any dielectric flux/current gate.
    PolarizationPoissonFVSolverCellCentered::SolveReport linearSolve;
    bool geometryAssemblyReused = false; // Same retained source and exact geometric/material inputs.
    // Actual center-phase queries for diagnostic D/epsilon recovery only.
    // Zero on reuse of the privately owned, fully keyed fixed geometry data.
    // Observational; does not affect the FV operator or any acceptance gate.
    std::size_t diagnosticCellPhaseQueries = 0;
    // Wall-clock observations only; never used by any acceptance condition.
    // The caller may aggregate these without changing field/diagnostic files.
    struct Timings
    {
        double input_s = 0.0, geometry_s = 0.0, linear_s = 0.0;
        double flux_s = 0.0, current_s = 0.0;
    } timings;
};

// The graph, retained scalar and charge panels must select the SAME active
// cells on a complete equal-level leaf lattice (not hanging-node AMR). Real
// internal/outer walls use explicit fitted BCs, never inactive cell values.
// Contact exemptions are only the panels' declared exact zero vertices; an
// unresolved/mixed boundary cell still rejects. Sampling filters solid cells
// from the fixed P3 windows, with the original fit gates and no lower order
// fallback. This active-patch path uses the original same-cell Hermite D3,
// not the full-domain reconstructed-D3 stencil (which crosses solid cells).
// No q=0 maximum-principle gate is imposed on a charged Poisson solution.
// Throws with stage context on failure; inputs and their storage are unchanged.
FittedSharpDielectricResult solveFittedSharpDielectricElectrostatics(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, CartesianGridScope scope,
    const DielectricConfig& config, const FittedPotentialInput& previousPotential,
    const SurfaceChargeSource& surfaceCharge,
    const DielectricConductivity* conductivity = nullptr);

// Repeated charged solves, e.g. the electric subproblem of an under-relaxed
// charge iteration. Retains the immutable source, owned P3 fits/sample windows
// and SuperLU factors, never a charge callback or another evolving geometry.
// Every call validates its complete inputs. Identical geometry/material data
// reuse the fits and cell-center diagnostic permittivity; affine charge terms
// are evaluated with the original formulas. The latter is not the FV face
// permittivity and never substitutes for the retained phase-face partitions.
// A new geometric fit uses the full validated linear assembly. Repeated fits
// change only affine constants of the privately owned stencil, and update the
// RHS of that owned linear operator without rebuilding its matrix. The live
// geometry/material/boundary checks above still run on EVERY public call.
// Unchanged full matrix values/row identities reuse the factorization. Any
// input, algebraic or physical flux/current failure clears BOTH caches.
// Use one instance per candidate/time step; do not call it concurrently.
// Charge time integration remains separate from this electric subproblem.
class FittedSharpDielectricSolver
{
public:
    FittedSharpDielectricSolver();
    ~FittedSharpDielectricSolver();
    FittedSharpDielectricResult solve(
        const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
        const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
        const LevelSetInterfacePanels& panels, CartesianGridScope scope,
        const DielectricConfig& config, const FittedPotentialInput& previousPotential,
        const SurfaceChargeSource& surfaceCharge,
        const DielectricConductivity* conductivity = nullptr);

    void clear();

    // Same charged P3/FV/SuperLU solve on actual Cartesian 2:1 leaves. Uses
    // retained native phase intervals and active physical sampling windows.
    // Pure-phase coarse/fine fluxes share Hydro's virtual-child gradient;
    // cut faces use Q(P3)-TPFA(P3) with the actual normal-projected distance.
    // Coarse/fine interpolation samples must be entirely active. Mixed-phase
    // point support uses the same P3 jump fit. No solid closure or time step is
    // supplied here. All existing fit, linear and current gates still apply.
    // Native wall phase certification covers the full face and half-link;
    // their zeros may only be exact declared contact endpoints. The remaining
    // part of a wall cell may be mixed. The old uniform whole-cell gate stays.
    FittedSharpDielectricResult solveNative(
        const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
        const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
        const LevelSetInterfacePanels& panels, const DielectricConfig& config,
        const FittedPotentialInput& previousPotential, const SurfaceChargeSource& surfaceCharge,
        const DielectricConductivity* conductivity = nullptr);

    // Complete native interface collection on one shared scalar: ONE global
    // Poisson/current solve with all component charge included in surfaceCharge.
    // Panel borrows last for this synchronous call. The collection must cover
    // every retained curve, without duplicates. An unresolved projection tie
    // between different components rejects; no nearest-child charge transfer
    // or independent per-drop global solve is implied by this entry.
    FittedSharpDielectricResult solveNative(
        const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
        const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
        const std::vector<const LevelSetInterfacePanels*>& surfaces, const DielectricConfig& config,
        const FittedPotentialInput& previousPotential, const SurfaceChargeSource& surfaceCharge,
        const DielectricConductivity* conductivity = nullptr);

private:
    FittedSharpDielectricResult solveImpl(
        const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
        const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
        const std::vector<const LevelSetInterfacePanels*>& surfaces, std::optional<CartesianGridScope> scope,
        const DielectricConfig& config, const FittedPotentialInput& previousPotential,
        const SurfaceChargeSource& surfaceCharge, const DielectricConductivity* conductivity,
        bool completeNativeCollection = false);
    struct PreparedAssembly;
    std::unique_ptr<PreparedAssembly> prepared_;
    PolarizationPoissonFVSolverCellCentered solver_;
};

// Physical Maxwell jumps in electric code units, on the ordered retained
// interface queries. No native velocity field, pressure gradient, dual volume
// or mechanical normalization is part of this electric result.
struct FittedSharpDielectricTractions
{
    std::shared_ptr<const InterfaceTraceGeometry> geometry;
    std::vector<DielectricInterfaceTraceReport> interfaceTraces;
    struct Diagnostics
    {
        std::size_t samplePoints = 0, localTraceFits = 0;
        double maxTraceCondition = 0.0, maxTraceReproductionError = 0.0;
        double maxNormalTraction = 0.0, maxTangentialTraction = 0.0;
        std::size_t axisCapExteriorSamples = 0;
        double maxAxisCapRelativeResidual = 0.0;
        std::size_t axisNeckHarmonicSamples = 0;
        double maxAxisNeckRelativeResidual = 0.0;
    } diagnostics;
    // Derived (normal,tangential) pairs for a receiving traction loader.
    // The authoritative physical values remain in interfaceTraces.
    std::vector<Vector2D> tractionComponents() const;
};

// Reconstruct the full Maxwell jump, including q*E_t, from the ACTUAL ABSOLUTE
// solved potential. Same P3 windows, charge jump, Cartesian stress identities
// and local D3 policy as before. Inactive potential payload is never read.
// Sampling contains source-local frames/raw meridional curvature derivatives;
// their retained panel incidence and local-jet provenance are checked. Native
// velocity support, pressure metrics and loading belong to the receiver.
// Throws with stage/sample context; inputs remain unchanged on failure.
// Independent sampling/factorizations may use OpenMP. Charge callbacks,
// ordered failure selection and publication remain on the calling thread.
FittedSharpDielectricTractions reconstructFittedSharpDielectricTractions(
    const QuadtreeMesh& mesh, std::shared_ptr<const InterfaceTraceGeometry> geometry,
    const LevelSetInterfacePanels& panels, CartesianGridScope scope,
    const DielectricConfig& config, const FittedElectricPotential& absolutePotential,
    const SurfaceChargeSource& surfaceCharge);

// Actual 2:1-leaf variant of the SAME physical trace implementation. Uses
// fixed physical sampling boxes and actual active centers, retaining original
// P3 rank/conditioning/jump rejection. Does not create pressure-side metrics.
// Explicit local closure for a declared unresolved RZ axis cap. The caller
// owns the event/model decision and supplies only its current axis-touching
// cells, bound to this same scalar. In those cells, fit the solved exterior
// potential with regular axisymmetric harmonic polynomials through degree 4;
// impose [Dn]=q and common Et to recover the liquid-side Maxwell stress.
// No 1/r Taylor coefficients, invented liquid samples, current update or
// change to the electric FV solve. Residual/model coverage are reported;
// physical event-model accuracy still requires grid/scale validation.
struct FittedAxisCapTractionModel
{
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    std::vector<int> cells;
};
FittedSharpDielectricTractions reconstructNativeFittedSharpDielectricTractions(
    const QuadtreeMesh& mesh, std::shared_ptr<const InterfaceTraceGeometry> geometry,
    const LevelSetInterfacePanels& panels, const DielectricConfig& config,
    const FittedElectricPotential& absolutePotential, const SurfaceChargeSource& surfaceCharge,
    const FittedAxisCapTractionModel* axisCapModel = nullptr);
