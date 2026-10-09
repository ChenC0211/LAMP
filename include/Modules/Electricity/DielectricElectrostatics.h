#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Geometry/InterfaceTraceGeometry.h"
#include "Modules/Electricity/DielectricPolarization.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Modules/Electricity/DielectricFluxChecks.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"

class QuadtreeMesh;
class InterfaceGeometry;
struct DielectricSharpAssemblyDiagnostics;
struct SurfaceChargeSource;

// All material, voltage, length and force quantities use the caller's
// consistent code units. This is not the legacy Electricity MMS/IIM model.
struct DielectricConfig
{
    std::string model = "perfectDielectric";
    // Explicit numerical choice, independent of the required SuperLU backend.
    // The default retains the original diffuse model and its output contract.
    std::string interfaceDiscretization = "diffuseTanh";
    double permittivityNegative = 2.0;
    double permittivityPositive = 1.0;
    std::string electrodeAxis = "y";
    double lowerPotential = 0.0;
    double upperPotential = 0.1;
    double linearResidualTolerance = 1e-9;
    // Physical checks have independent defaults. Only the legacy input adapter
    // expands a shared linear tolerance; unused explicit controls are rejected.
    // Diffuse boundary flux: 1e-9; sharp reconstruction/Gauss: 1e-9 each;
    // ordinary leakyDielectric charge/electrode closure: 16e-9.
    std::optional<double> boundaryFluxTolerance;
    std::optional<DielectricFluxChecks> sharpFluxChecks;
    std::optional<double> chargeClosureTolerance;
    // Explicit leakyDielectric App inputs; not legacy electrospray constants.
    double conductivityNegative = 0.0;
    double conductivityPositive = 0.0;
    double initialSurfaceChargeDensity = 0.0;
    bool chargeConvection = false;
    // Bound by the independent electric-domain reader or by EHD from the
    // canonical Hydro coordinateSystem. Must match the retained scalar geometry.
    bool axisymmetric = false;
    // Explicit native thin-column closure. Requires a certified liquid axis
    // throughout the fixed sampling box and radially single-valued interface.
    // Uses regular interior/logarithmic exterior harmonics when no liquid
    // center is resolved. Not a generic single-phase P3 fallback.
    bool nativeAxisNeckHarmonics = false;
};

// Numerical capability/material validation, also used by the input adapter.
void validateDielectricConfig(const DielectricConfig& config);

struct DielectricConductivity
{
    double negative = 0.0;
    double positive = 0.0;
};

// Currents from this solve's actual shared FV fluxes, integrated per unit
// depth in planar mode, or over revolved faces in RZ. Cell entries include
// BOTH interface and bulk cells, in solver order;
// the caller must not silently remap bulk current divergence to the interface.
struct DielectricConductionReport
{
    std::vector<double> cellIntegratedCurrent; // -sum(outward J), positive accumulation.
    double lowerBoundaryCurrent = 0.0, upperBoundaryCurrent = 0.0;
    double netBoundaryCurrent = 0.0, totalCellCurrent = 0.0;
    double closureResidual = 0.0, relativeClosure = 0.0;
    double maxDisplacementSplitResidual = 0.0;
    double absoluteFaceCurrent = 0.0; // Sum |J| once per physical face, including electrodes.
};

struct DielectricReport
{
    bool success = false;
    bool usedSuperLU = false;
    std::string stage;
    std::string reason;
    double relativeResidual = 0.0;
    // Outward D=epsilon*E flux through the low/high electrodes. All other
    // sides are homogeneous Neumann. Flux is per unit depth in planar mode,
    // and integrated over the revolved face (2*pi*r*ds) in RZ mode.
    double lowerBoundaryFlux = 0.0;
    double upperBoundaryFlux = 0.0;
    double netBoundaryFlux = 0.0;
    double boundaryFluxImbalance = 0.0;
    double fieldEnergy = 0.0;
    double maxField = 0.0;
    double minPermittivity = 0.0;
    double maxPermittivity = 0.0;
    double minPotential = 0.0;
    double maxPotential = 0.0;
    // Sharp-path audit only; existing diffuse diagnostics/CSV are unchanged.
    bool usedSharpInterface = false;
    std::size_t correctedPhysicalFaces = 0;
    std::size_t certifiedBoundaryFaces = 0;
    // Zero-area symmetry faces have regularity/metric checks, NOT a no-contact certificate.
    std::size_t zeroAreaAxisFaces = 0;
    std::size_t zeroNeumannFaces = 0;
    std::size_t localTraceFits = 0;
    double cellIntegratedFluxResidual = 0.0;
    double maxTraceCondition = 0.0;
    double maxTraceReproductionError = 0.0;
    std::size_t axisNeckHarmonicFaces = 0;
    bool usedSurfaceCharge = false;
    double prescribedSurfaceCharge = 0.0;
    double surfaceChargeClosureResidual = 0.0;
    // Immutable result of THIS successful solve only; failure or no request
    // leaves nullptr. There is no model-level last-current cache to go stale.
    std::shared_ptr<const DielectricConductionReport> conductionCurrent;
};

// Actual local traces in the retained InterfaceTraceGeometry query order. These
// diagnostics contain no reconstruction implementation or reference types.
struct DielectricInterfaceTraceReport
{
    double interfacePotential = 0.0;
    double tangentialElectricField = 0.0;
    double normalDisplacement = 0.0;
    double tangentialFieldDerivative = 0.0;
    double normalDisplacementDerivative = 0.0;
    double tangentialFieldSecondDerivative = 0.0;
    double normalDisplacementSecondDerivative = 0.0;
    double normalTraction = 0.0;
    double condition = 0.0;
    double reproduction = 0.0;
    std::vector<int> sampleCellIndices;
    // With prescribed charge, normalDisplacement and its derivatives above
    // are negative-side values. These additional fields are sharp diagnostics.
    double positiveNormalDisplacement = 0.0;
    double surfaceChargeDensity = 0.0;
    double tangentialTraction = 0.0;
    // Explicit sub-resolution axis-cap closure, never an automatic P3
    // fallback. Fits the actual exterior potential with regular RZ harmonics.
    bool axisCapExteriorModel = false;
    bool axisNeckHarmonicModel = false;
    double modelRelativeResidual = 0.0;
};

// Owned publication of one successful electrostatic solve. No legacy Field
// pointer escapes: const vector elements are read-only and retaining this
// snapshot never observes a later solve. Native indices belong to meshIdentity.
struct DielectricFields
{
    NativeMeshIdentity meshIdentity;
    std::vector<double> potential, permittivity, ex, ey;
    std::vector<double> levelSet; // Actual nodal input, retained for coherent output.
    bool axisymmetric = false;
    std::shared_ptr<const DielectricPolarization> polarization;
    std::shared_ptr<const InterfaceTraceGeometry> traceGeometry;
    std::vector<DielectricInterfaceTraceReport> interfaceTraces;
    explicit DielectricFields(const QuadtreeMesh& mesh) : meshIdentity(mesh) {}
    std::vector<Vector2D> interfaceTractions() const;
};
struct DielectricSolveResult
{
    DielectricReport report;
    // Present only when report.success is true. Failure never returns the
    // previous solve's fields; retained successful results remain unchanged.
    std::shared_ptr<const DielectricFields> fields;
};

// Quasi-static linear dielectric on a fixed uniform Cartesian mesh, both
// phases included. Sharp mode additionally supports RZ interfaces, including
// smooth axis contact with explicitly parity-enabled shared geometry and
// resolved positive-radius traces. Prescribed charge must match the area measure.
// No EB, AMR, periodic BC or emission. Conductive time evolution is owned by
// Electricity's evolution services, not this electrostatic solve.
//
// Diffuse mode publishes cell permittivity and the physical polarization
// coefficient; sharp mode publishes actual-potential two-sided P3 traces and
// Maxwell traction. Receiving force discretization belongs to Hydro.
// Complete-domain scalar geometry has no mechanical support information.
// Optional trace queries are retained on that exact source; omitting them
// solves potential/flux/current without constructing a momentum load band.
// An optional SurfaceChargeSource adds prescribed free surface charge; this
// solve does not evolve conduction or apply momentum loads.
class DielectricElectrostatics
{
public:
    explicit DielectricElectrostatics(const QuadtreeMesh& mesh);

    // The mesh must be the unchanged mesh supplied to the constructor.
    // On failure report.success=false and fields=nullptr. All output belongs
    // to the returned result, never to mutable public solver storage.
    // Optional caller-owned assembly observations are sharp-only, allocated
    // before commit and replaced only on success; nullptr collects nothing.
    // Optional conductivity requests conservative currents from this sharp
    // solve. Both phase values must be finite and nonnegative; an insulating
    // phase may be exactly zero. The returned report owns immutable currents.
    DielectricSolveResult solve(const QuadtreeMesh& mesh,
                           const VertField<double>& levelSet,
                           double thickness,
                           const DielectricConfig& config,
                           const InterfaceGeometry* sharedGeometry = nullptr,
                           DielectricSharpAssemblyDiagnostics* assemblyDiagnostics = nullptr,
                           const SurfaceChargeSource* surfaceCharge = nullptr,
                           const DielectricConductivity* conductivity = nullptr,
                           std::shared_ptr<const InterfaceTraceGeometry> traceQueries = nullptr);


private:
    const QuadtreeMesh* mesh_ = nullptr;
    std::uint64_t meshSignature_ = 0;
    FaceSegmentGraph graph_;
    PolarizationPoissonFVSolverCellCentered solver_;
    // Called only after common config/mesh/input validation.
    // Uses the existing graph and solver, and publishes a complete immutable result.
    void solveSharp(const QuadtreeMesh& mesh, const VertField<double>& levelSet,
                    const DielectricConfig& config,
                    const InterfaceGeometry& sharedGeometry,
                    DielectricReport& report,
                    DielectricSharpAssemblyDiagnostics* assemblyDiagnostics,
                    const SurfaceChargeSource* surfaceCharge,
                    const DielectricConductivity* conductivity,
                    const std::shared_ptr<const InterfaceTraceGeometry>& traceQueries,
                    std::shared_ptr<const DielectricFields>& publication);
};
