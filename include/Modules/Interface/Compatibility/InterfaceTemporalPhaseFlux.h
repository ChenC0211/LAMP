#pragma once
#include "Modules/Interface/Compatibility/InterfacePhaseFlux.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/ActiveHermitePatch.h"
#include <array>

namespace InterfaceTemporalPhaseFlux
{
struct Options
{
    double relativeTolerance=1e-8;
    double capacityTolerance=1e-13;
    int maxDepth=20, maxEvaluations=4096;
};
struct Integral
{
    // Time-integrated physical volume, NOT instantaneous volume/time.
    InterfacePhaseFlux::FluxBudget volume;
    std::vector<long double> absoluteTransportByCell, errorByNativeFace;
    // Estimated GL4/GL8 error plus floating accumulation uncertainty, summed
    // over incident shared faces. This is not a certified interval enclosure.
    std::vector<long double> integrationUncertaintyByCell;
    int evaluations=0, intervals=0;
};
template<class Sensitivity> struct PhaseResult
{
    Integral integral;
    // d(integrated volume flux)/d(current independent scalar data).
    Sensitivity derivative;
};
using Result=PhaseResult<InterfacePhaseFlux::QuarticFluxSensitivity>;
using CubicResult=PhaseResult<InterfacePhaseFlux::FluxSensitivity>;
template<class Patch> struct PhaseTrajectory
{
    // Frozen physical predictor at theta=0,1/4,1/2,3/4,1. Correction at an
    // interior time is theta*(candidate-predictorEnd), so its Jacobian uses
    // theta. Data describe one old state and one velocity time pair.
    std::array<Patch,5> samples;
};
using Trajectory=PhaseTrajectory<Quantum::Geometry::ActiveQuarticHermitePatch>;
using CubicTrajectory=PhaseTrajectory<Quantum::Geometry::ActiveHermitePatch>;
// Interpolate the SAME Q4 independent data and frozen velocity traces in time.
// Split at endpoint phase events, integrate with adaptive GL4/GL8, and rebuild
// phase cuts for every sample. The new-data derivative includes its theta
// weight. Old/new patches, sources and traces must describe the same native
// mesh and time pair; this service does not choose dt or commit a state.
Result integrate(const QuadtreeMesh& mesh,const InterfaceTransportDomain& domain,
    const Quantum::Geometry::ActiveQuarticHermitePatch& oldPatch,
    const Quantum::Geometry::ActiveQuarticHermitePatch& currentPatch,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& oldSource,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& currentSource,
    const InterfaceNormalFluxView& oldVelocity,const InterfaceNormalFluxView& currentVelocity,
    double dt,const Options& options = {},const Trajectory* trajectory = nullptr);
// Same integration and shared-face ledger, retaining cubic nodal identities.
// Declared sector identities stay fixed; their jets follow the prescribed
// predictor in time and remain fixed under correction variation. The caller
// still enforces pinned values and shared trace validity.
CubicResult integrate(const QuadtreeMesh& mesh,const InterfaceTransportDomain& domain,
    const Quantum::Geometry::ActiveHermitePatch& oldPatch,
    const Quantum::Geometry::ActiveHermitePatch& currentPatch,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& oldSource,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& currentSource,
    const InterfaceNormalFluxView& oldVelocity,const InterfaceNormalFluxView& currentVelocity,
    double dt,const Options& options = {},const CubicTrajectory* trajectory = nullptr);
}
