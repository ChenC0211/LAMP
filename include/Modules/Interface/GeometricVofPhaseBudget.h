#pragma once
#include "Modules/Interface/GeometricVofTransport.h"
#include "Modules/Interface/VofPhaseFraction.h"
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace LAMP::Interface::Detail
{
// Stateless arithmetic kernel. It consumes private candidate face data; it
// owns neither a mesh/geometry nor an accepted Interface or Hydro state.
struct GeometricVofPhaseBudgetCell
{
    long double capacity=0;
    double alpha=0; // Legacy input, or exact rounded view when fraction is set.
    std::optional<VofPhaseFraction> fraction; // Value only, no State/geometry owner.
};
struct GeometricVofDirectionalWorkMoments
{
    // [positive/negative branch][1,2*s-1] spatial work basis. The negative
    // branch coefficients themselves are NOT negated. Duration stays outside.
    std::array<std::array<long double,2>,2> gross{},liquid{};
    // liquid is the actual HALF-DEPTH back-line sharp midpoint, not a theta
    // AVF integral. The retained sharp finite extent fixes its Gonzalez work.
    std::array<long double,2> zeroTraceLiquid{};
    long double mean=0,halfVariation=0,duration=0;
};
struct GeometricVofAcceptedWorkMoments
{
    std::array<long double,2> value{}; // Exact branch sum, rounded only once.
    std::array<std::array<long double,2>,2> branch{};
    long double rawWorkRelativeResidual=0,roundingWorkResidual=0;
    bool zeroTrace=false;
};
struct GeometricVofPhaseBudgetResult;
// Immutable acceptance arithmetic only. Geometry and authentic owner/clock
// binding remain with the Transport snapshot; raw arguments are actual old-fan
// phase integrals, never locally normalized or accepted phase diagnostics.
class GeometricVofPhaseBudgetRestriction final
{
public:
    GeometricVofIntegratedFlux measure(std::size_t row,
        const GeometricVofIntegratedFlux& raw) const;
    // Same retained high/raw geometry, donor and limiter. Optional cleanup
    // phases are actual independent [branch][liquid,gas] values; preserve their
    // minor and restore the dominant from the retained exact branch gross.
    // roundingWorkResidual is signed volume: duration*dot(native trace,value)
    // minus the exact accepted signed liquid extent, not a new flux/stock.
    GeometricVofAcceptedWorkMoments liquidWorkMoments(std::size_t row,
        const GeometricVofDirectionalWorkMoments&,
        std::optional<std::array<std::array<long double,2>,2>> acceptedP0Phases=std::nullopt) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit GeometricVofPhaseBudgetRestriction(std::shared_ptr<const Data>);
    friend GeometricVofPhaseBudgetResult prepareGeometricVofPhaseBudget(
        const std::vector<GeometricVofPhaseBudgetCell>&,
        const std::vector<GeometricVofTransportFace>&,bool,long double,
        const std::vector<GeometricVofTransportFace>*,long double,const std::vector<double>*,bool);
};
struct GeometricVofPhaseBudgetResult
{
    // Independent phase inventories; both positive values remain observable
    // even when a later double alpha cannot represent their ratio.
    std::vector<std::array<long double,2>> phaseInventory; // liquid, gas
    std::vector<GeometricVofIntegratedFlux> flux;
    std::shared_ptr<const GeometricVofPhaseBudgetRestriction> restriction;
    std::size_t limitedFaceExtents=0;
    long double maximumOutflowCfl=0,maximumTotalDivergence=0;
    long double maximumCapacityReconciliation=0;
};
// Opt-in arithmetic for a signed total measured independently of two rounded
// counterflow gross diagnostics. Keeps each independent branch minor phase;
// the dominant phase is its exact complement before signed phase export.
// No owner, density, limiter or geometry is created; mass is left zero.
GeometricVofIntegratedFlux reconcileGeometricVofNativeSignedTotal(
    const GeometricVofIntegratedFlux&);
GeometricVofPhaseBudgetResult prepareGeometricVofPhaseBudget(
    const std::vector<GeometricVofPhaseBudgetCell>&,
    const std::vector<GeometricVofTransportFace>&,bool half,
    long double maximumOutflowCfl=.9L,
    const std::vector<GeometricVofTransportFace>* geometry=nullptr,
    long double relativeConservationTolerance=1e-12L,
    // Explicit boundary INPUT only, never an exterior inventory. Null retains
    // the original impermeable contract. Callers authenticate inflow support.
    const std::vector<double>* exteriorDonorAlpha=nullptr,
    // Explicit typed native-mean law only. False preserves the old arithmetic.
    bool preserveNativeSignedTotal=false);
}
