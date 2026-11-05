#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Interface
{
struct HermiteGalsFaceVelocity
{
    enum class Boundary { Interior, Impermeable, Open };
    int nativeFace = -1; // VEdge indices followed by HEdge indices
    Boundary boundary = Boundary::Interior;
    // Physical +x/+y normal velocity, ascending powers of the ORIGINAL
    // half-face t in [0,1]. This is a local supplied reconstruction, never
    // a new transported/Hydro unknown. XY/RZ metric is added by this module.
    std::array<std::array<double,11>,2> halfPolynomial{};
};
class HermiteGalsFluxVelocity
{
public:
    HermiteGalsFluxVelocity(const QuadtreeMesh& mesh,const HermiteGalsState& midpoint,
        const std::vector<HermiteGalsFaceVelocity>& faces);
    HermiteGalsIdentity geometryIdentity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const std::vector<HermiteGalsFaceVelocity>& faces() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};
struct HermiteGalsPhaseFluxConfig
{
    double relativeIntegrationTolerance = 1e-13;
    double relativeDivergenceTolerance = 1e-12; // |dt sum s Q|/capacity
    double relativeConservationTolerance = 1e-12; // global closure / total capacity
};
struct HermiteGalsPhaseFaceFlux
{
    int nativeFace = -1, negativeCell = -1, positiveCell = -1;
    long double totalPositiveRate = 0, totalNegativeRate = 0;
    long double liquidPositiveRate = 0, liquidNegativeRate = 0;
    // Independent subset quadrature before error-budgeted reconciliation.
    // Preserve the original estimate for diagnostics.
    long double integratedLiquidPositiveRate = 0, integratedLiquidNegativeRate = 0;
    // After integration consistency, before conservative inventory repair.
    long double rawLiquidPositiveRate = 0, rawLiquidNegativeRate = 0;
    long double totalRate = 0, liquidRate = 0;
    long double liquidVolume = 0; // dt*liquidRate, used identically by both owners
    long double integrationErrorEstimate = 0;
};
struct HermiteGalsPhaseFluxReport
{
    std::size_t activeFaces=0,internalFaces=0,counterflowFaces=0,phaseRanges=0;
    std::size_t roundoffRepairSweeps=0,repairedFaces=0;
    std::size_t subsetReconciliationCount=0,directionalIdentityCount=0;
    long double maxRelativeSubsetAdjustment=0; // / half-face total gross rate
    long double liquidBefore=0,liquidAfter=0,boundaryOutward=0,closure=0;
    long double maxOutflowCfl=0,maxRelativeDivergence=0;
    long double maxRelativeFluxCorrection=0,integrationErrorEstimate=0;
};
class HermiteGalsPhaseFluxResult
{
public:
    HermiteGalsIdentity previousIdentity() const;
    HermiteGalsIdentity geometryIdentity() const;
    HermiteGalsTimeLayer timeLayer() const; // midpoint interval defining Q
    HermiteGalsGeometryView geometry() const;
    const std::vector<HermiteGalsPhaseFaceFlux>& faces() const;
    const std::vector<HermiteGalsCellVolume>& phaseVolumes() const;
    const HermiteGalsPhaseFluxReport& report() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsPhaseFluxResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsPhaseFluxResult advanceHermiteGalsPhaseVolumes(const QuadtreeMesh&,
        const HermiteGalsState&,const HermiteGalsState&,const HermiteGalsFluxVelocity&,
        const HermiteGalsPhaseFluxConfig&);
};
// One midpoint FV interval. Source inventories and geometry are never changed.
// Returns unique rates and updated M, NOT geometry, a new accepted State, or
// a legacy upwind/FCT fallback. HermiteGalsReflux exposes sums of these SAME
// final actual subfaces; it must not apply a second inventory correction.
HermiteGalsPhaseFluxResult advanceHermiteGalsPhaseVolumes(const QuadtreeMesh& mesh,
    const HermiteGalsState& previous,const HermiteGalsState& midpoint,
    const HermiteGalsFluxVelocity& velocity,const HermiteGalsPhaseFluxConfig& config = {});
}
