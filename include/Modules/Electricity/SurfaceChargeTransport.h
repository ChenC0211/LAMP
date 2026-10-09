#pragma once

#include <cstddef>
#include <string>
#include <vector>
#include <optional>

namespace Quantum::Electricity
{
// External upwind charge densities at the two ends of an OPEN chain. Face
// fluxes point along the chain: start F>0 and end F<0 are inflows. A missing
// density is legal only when that endpoint is not an inflow. No default q=0.
struct SurfaceChargeBoundaryState
{
    std::optional<double> startDensity;
    std::optional<double> endDensity;
};

struct SurfaceChargeStepReport
{
    bool success = false;
    std::string stage;
    std::string reason;
    double totalChargeBefore = 0.0;
    double totalChargeAfter = 0.0;
    double integratedNormalCurrent = 0.0;
    // Charge/time, integrated in space. In/out refer to the direction of
    // surface-measure transport; signed charge makes either value negative.
    double boundaryChargeInflow = 0.0;
    double boundaryChargeOutflow = 0.0;
    double netBoundaryCurrent = 0.0; // inflow - outflow, separate from normal I.
    double conservationResidual = 0.0;
    double relativeConservationResidual = 0.0;
    double conservationScale = 0.0; // sum |Q_old| + dt*(sum |I| + |Fq| at both ends).
    double maxAdvectiveCfl = 0.0;
    double maxAbsDensity = 0.0;
    // Closed: N faces from panel i into (i+1)%N. Open: N+1 boundaries,
    // boundary i points from i-1 into i. End fluxes are zero for the old API.
    std::vector<double> advectiveChargeFlux;
};

// Conservative charge state on one closed loop or open chain. Open endpoints
// have zero flux by default; nonzero flux requires the explicit boundary API.
// The geometry owner supplies panel measures and stable panel identities.
// No private mesh, level set or clock lives here. A geometry owner may supply
// an orientation-preserving periodic coordinate map for conservative remapping.
// Stored unknowns are INTEGRATED charges, not interpolated charge densities.
class SurfaceChargeTransport
{
public:
    struct PanelWeight {std::size_t targetPanel=0;long double share=0;};
    // Geometry-owner supplied positive conservative weights, one row per old
    // panel. Geometry coverage/orientation is the owner's responsibility;
    // this kernel checks stochastic rows and that every target is represented.
    // No normalization, current, extra tangential flux or density interpolation.
    SurfaceChargeStepReport remapWeights(const std::vector<std::vector<PanelWeight>>&,
        const std::vector<double>& newPanelMeasure);
    struct PanelCoordinate
    {
        std::size_t panelIndex = 0;
        double localMeasure = 0.0; // ds or RZ area measured from this NEW panel's begin.
    };
    // Requires at least three panels, positive finite measures, finite signed
    // densities and representable integrated charges. Throws on invalid input.
    SurfaceChargeTransport(const std::vector<double>& panelMeasure,
                           const std::vector<double>& initialDensity, bool closed = true);

    // Copies integrated charges directly, without a Q->density->Q round trip.
    // The implied densities must still be finite and representable.
    static SurfaceChargeTransport fromIntegratedCharge(
        const std::vector<double>& panelMeasure,
        const std::vector<double>& integratedCharge, bool closed = true);

    // Piecewise-constant conservative remapping on the NEW arclength axis.
    // mappedOldBoundaries has oldCount+1 strictly increasing entries spanning
    // one period, with an arbitrary periodic origin. newBoundaries starts at
    // exact zero and ends at the positive period; newCount may differ (>=3).
    // Each old Q is apportioned by overlap/mappedOldLength. There is no global
    // charge correction or normalization. Winding/coverage are checked to
    // fixed roundoff; unresolved reductions, folding and underflow fail loudly.
    // On failure both existing state arrays remain unchanged.
    // Open chains instead require both boundary arrays to start at exact zero
    // and end at the same value; overlap is nonperiodic and never wraps.
    SurfaceChargeStepReport remap(const std::vector<double>& mappedOldBoundaries,
                                 const std::vector<double>& newBoundaries);

    // Explicit OPEN-chain coordinate adapter for the same overlap rule. New
    // measures are supplied directly, never differenced from a whole-chain
    // cumulative double axis. Each old boundary is located in one new panel;
    // (j,newPanelMeasure[j]) and (j+1,0) are the SAME shared boundary identity.
    // Both final representations (N-1,A_last) and (N,0) are accepted. No other
    // coordinate with panelIndex==N is valid. Exact first/final coverage and
    // strict order are checked after identity normalization, without snapping.
    // Donor lengths are compensated FLOATING sums of local overlaps; every
    // positive overlap is retained separately for Q_old*overlap/donorLength.
    // No wider long-double precision, exact denominator, root merging, global
    // charge normalization or relaxed coverage/conservation gate is assumed.
    // Closed chains reject this overload. Failure preserves both state arrays.
    // This kernel alone does not implement moving geometry or Field/App remap.
    SurfaceChargeStepReport remap(const std::vector<PanelCoordinate>& mappedOldBoundaries,
                                 const std::vector<double>& newPanelMeasure);

    // Explicit upwind finite-volume transport on a moving interface partition:
    //   Q_i(new) = Q_i + dt*(I_i + F_(i-1) - F_i).
    // Face i points from i to (i+1)%N. relativeAdvectiveMeasureFlux is surface
    // measure/time relative to the panel boundary motion (planar: velocity).
    // I_i is normal current INTO panel i, integrated over its measure; the
    // coupled electric solver must supply it. Normal points negative->positive,
    // hence I_i = integral(J_negative.n - J_positive.n) on that panel.
    // This kernel does not solve the electric field or make stiff conduction
    // stable. The caller owns conduction time integration and its restrictions.
    // Stretching changes density via Q/newMeasure, without changing Q itself.
    // Panel count/order must not change. Failure leaves both state arrays intact.
    // Open chains use N+1 fluxes, exact zero at both ends, and update with
    // Q_i(new)=Q_i+dt*(I_i+F_i-F_(i+1)); no periodic neighbors are used.
    SurfaceChargeStepReport advance(
        const std::vector<double>& nextPanelMeasure,
        const std::vector<double>& relativeAdvectiveMeasureFlux,
        const std::vector<double>& integratedNormalCurrent,
        double dt);

    // Explicit open-boundary variant. Nonzero relative end measure fluxes
    // use an external supplied density for inflow, the old endpoint-panel
    // density for outflow. The CFL includes losses through either end.
    // Closed chains reject this overload. The old overload still requires
    // exact zero end fluxes; remap and material endpoint coverage are unchanged.
    SurfaceChargeStepReport advance(
        const std::vector<double>& nextPanelMeasure,
        const std::vector<double>& relativeAdvectiveMeasureFlux,
        const std::vector<double>& integratedNormalCurrent,
        double dt, const SurfaceChargeBoundaryState& boundary);

    const std::vector<double>& panelMeasure() const noexcept { return measure_; }
    const std::vector<double>& panelCharge() const noexcept { return charge_; }
    bool closed() const noexcept { return closed_; }
    std::vector<double> density() const;

private:
    SurfaceChargeStepReport advanceImpl(
        const std::vector<double>& nextPanelMeasure,
        const std::vector<double>& relativeAdvectiveMeasureFlux,
        const std::vector<double>& integratedNormalCurrent,
        double dt, const SurfaceChargeBoundaryState* boundary);
    struct IntegratedChargeTag {};
    SurfaceChargeTransport(const std::vector<double>& panelMeasure,
        const std::vector<double>& integratedCharge, IntegratedChargeTag, bool closed);
    std::vector<double> measure_;
    std::vector<double> charge_;
    bool closed_ = true;
};
}
