#pragma once

#include "Modules/Electricity/FittedElectricPotential.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Solvers/Poisson/PolarizationPoissonFVSolverCellCentered.h"
#include <memory>
#include <vector>

namespace Quantum::Electricity {
// q=0 cell-material harmonic FV problem. Bounds use an explicitly supplied
// reference scale: preserve the caller's comparison arithmetic without making
// the numerical module depend on SI, dictionaries or an App material table.
struct FittedCellHarmonicControls
{
    double linearTolerance = 0.0;
    double potentialBoundsTolerance = 0.0;
    double potentialScale = 1.0;
    double lowerPotential = 0.0, upperPotential = 0.0;
};
struct FittedElectricFaceValue { double voltageDropCode, distanceCode, chargeCode; };

class FittedCellHarmonicState
{
public:
    FittedCellHarmonicState() = default;
    explicit operator bool() const noexcept { return bool(data_); }
    // A binding prepared for a charged bootstrap/AMR seed is explicitly
    // unsolved. Only solve success publishes a Solved stage and report.
    enum class Stage { Prepared, Solved };
    Stage stage() const;
    const FittedElectricPotential& potential() const;
    const NativeLeafDomain& geometry() const;
    const FaceSegmentGraph& graph() const;
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundaries() const;
    const std::vector<double>& epsilonByCell() const;
    const std::vector<double>& epsilonBySegment() const;
    const PolarizationPoissonFVSolverCellCentered::SolveReport& report() const;
    // Physical stamp of a Solved state; Prepared seeds carry no solved time.
    double time() const;
    void validateMesh(const QuadtreeMesh&) const;
    // One validated batch, with the same fitted owner distance and RZ
    // multiplication order as the original diagnostic. No per-face mesh scan.
    std::vector<FittedElectricFaceValue> faceValues(const QuadtreeMesh&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    friend struct FittedCellHarmonicAccess;
};

// Own the fixed graph/BCs and actual materials. Empty graph rows select the
// inactive cells; their finite payload is set to inactivePotential. No phase
// ID, voltage conversion or physical boundary name is inferred here.
// With controls, solve before publishing (no extra initial-field copy).
// Without controls, publish Prepared only, for a subsequent charged solve.
FittedCellHarmonicState prepareFittedCellHarmonic(const QuadtreeMesh&, FaceSegmentGraph,
    std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>,
    std::vector<double> epsilonByCell, bool axisymmetric, double inactivePotential,
    const FittedPotentialInput* seed = nullptr, const FittedCellHarmonicControls* controls = nullptr);

// Same fixed domain/BCs, new cell materials and optional retained scalar.
// Shares fixed problem data; copies only the original working potential.
// Failure leaves all retained states usable. time is a caller-owned stamp,
// not a second clock; same-time synchronization is allowed, reversal is not.
FittedCellHarmonicState updateFittedCellHarmonic(const QuadtreeMesh&, const FittedCellHarmonicState&,
    std::vector<double> epsilonByCell, const FittedCellHarmonicControls&, double time,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source = nullptr);
}
