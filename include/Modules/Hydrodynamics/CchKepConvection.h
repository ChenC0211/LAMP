#pragma once

#include "Modules/Hydrodynamics/CchDualMaterial.h"
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
struct CchConvectionWeight { int velocity=-1; double weight=0; };
// Prepared final-operator correction: J_ij=coefficient, J_ji=-coefficient.
// Not another material face flux; both indices belong to one local patch.
struct CchConvectionCoupling { int first=-1, second=-1; long double coefficient=0; };
// Explicit exterior segment relation: second translated into first's image.
// Geometry/port data only; no extra velocity unknown or material flux.
struct CchConvectionPeriodicPair
{ int firstSegment=-1, secondSegment=-1; Vector2D translation{}; };
struct CchConvectionMomentReport
{
    int interiorRows=0, constraints=0, iterations=0;
    long double constantError=0, affineError=0, quadraticError=0;
    int boundedRows=0;
    long double coefficientRatio=0, quadraticBudgetRatio=0, maximumSupportRatio=0;
    int projectionCalls=0;
    long double quadraticBoundError=0;
};
struct CchKepComponentReport
{
    long double momentumRate=0, boundaryMomentumRate=0, convectivePower=0;
    long double boundaryKineticPower=0, materialKineticPower=0, kineticRate=0;
    long double momentumIdentity=0, kineticIdentity=0;
    long double dissipationPower=0, maximumSensor=0;
    int dissipationConnections=0, activeDissipationConnections=0;
};
struct CchKepConvectionReport { std::array<CchKepComponentReport,2> components{}; };

// Explicit local stabilization. Density contrast alone never activates it.
// The default build overloads remain the unchanged center-only operator.
struct CchConvectionDissipationConfig
{
    double strength=1, activation=0.2, saturation=0.5;
};
struct CchConvectionDissipationSample
{
    int first=-1, second=-1, firstSegment=-1, secondSegment=-1;
    long double sensor=0, chi=0, halfJump=0, momentumRate=0, power=0;
};

// Explicit transported-component traces on external dual segments. Internal
// entries must be zero. Bound to the actual H09 owner, identities and clock;
// no default boundary law, hidden ghosts or accepted state.
class CchConvectionBoundary
{
public:
    static CchConvectionBoundary bind(const QuadtreeMesh&, CchDualMaterial,
        const std::vector<double>& segmentTraces, CchMaterialSources,
        Interface::HermiteGalsTimeLayer);
    const CchDualMaterial& material() const noexcept;
    const std::vector<double>& traces() const noexcept;
    void validateFor(const QuadtreeMesh&, const CchDualMaterial&, CchMaterialSources,
        Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit CchConvectionBoundary(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

class CchKepConvectionResult
{
public:
    const CchDualMaterial& material() const noexcept;
    const CchConvectionBoundary& boundary() const noexcept;
    // Raw F*t^T*u diagnostic, NOT the final compatible native scatter.
    const std::vector<long double>& centerMomentumRates() const noexcept;
    // Positive outward convection N, in native velocity order.
    const std::vector<long double>& momentumRates() const noexcept;
    // -N + u_mid*R_dual, excluding pressure/viscosity/external loads.
    const std::vector<long double>& momentumLoads() const noexcept;
    // Actual native stabilization L, including periodic mass redistribution.
    const std::vector<long double>& dissipationRates() const noexcept;
    const std::vector<CchConvectionDissipationSample>& dissipationSamples() const noexcept;
    const CchKepConvectionReport& report() const noexcept;
private:
    friend class CchKepConvection;
    struct Data;
    explicit CchKepConvectionResult(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

// Common-position diagnostic plus final native spatial/skew moment constraints.
// Consumes the unique H09 flux; optional explicit local nonnegative dissipation.
// No mass/velocity update, IO or clock ownership.
class CchKepConvection
{
public:
    // Raw phase-volume Galerkin A_ij = integral B_i dot(F dot grad B_j).
    // Uses the unique flux's physical view and the same full [U,q] source.
    // No forced skew projection or R/boundary cancellation. General temporal,
    // boundary and continuity pairing remains the caller's obligation.
    static CchKineticHodge::Matrix modalBulkTransportMatrix(const QuadtreeMesh&,
        const CchKineticHodge&,const CchMaterialFlux&);
    static CchKepConvection build(const QuadtreeMesh&, CchDualMaterial);
    static CchKepConvection build(const QuadtreeMesh&, CchDualMaterial,
        const std::vector<CchConvectionPeriodicPair>&);
    static CchKepConvection buildControlled(const QuadtreeMesh&, CchDualMaterial,
        CchConvectionDissipationConfig);
    static CchKepConvection buildControlled(const QuadtreeMesh&, CchDualMaterial,
        const std::vector<CchConvectionPeriodicPair>&, CchConvectionDissipationConfig);
    const CchDualMaterial& material() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const std::vector<Vector2D>& points() const noexcept;
    const std::vector<std::vector<CchConvectionWeight>>& weights() const noexcept;
    const std::vector<CchConvectionCoupling>& couplings() const noexcept;
    const CchConvectionMomentReport& momentReport() const noexcept;
    const std::vector<CchConvectionPeriodicPair>& periodicPairs() const noexcept;
    const std::vector<int>& representatives() const noexcept;
    CchKepConvectionResult apply(const QuadtreeMesh&, const CchNativeValues&,
        const CchConvectionBoundary&, CchMaterialSources actualSources,
        Interface::HermiteGalsTimeLayer actualMidpointClock) const;
    CchKepConvectionResult applyPaired(const QuadtreeMesh&, const CchNativeValues&,
        CchMaterialSources actualSources, Interface::HermiteGalsTimeLayer actualMidpointClock) const;
private:
    struct Data;
    explicit CchKepConvection(std::shared_ptr<const Data>);
    static CchKepConvection buildWithPorts(const QuadtreeMesh&, CchDualMaterial,
        const std::vector<CchConvectionPeriodicPair>*, const CchConvectionDissipationConfig* = nullptr);
    std::shared_ptr<const Data> data_;
};
}
