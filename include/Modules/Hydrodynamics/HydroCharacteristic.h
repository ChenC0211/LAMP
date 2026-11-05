#pragma once

#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include <optional>

namespace LAMP::Hydrodynamics
{
namespace Detail
{
// Six physical value/scaled-J traces, sampled at t=0,1/4,1/2,3/4,1.
// Complete for the actual piecewise Q4 field. Exposed internally for negative controls.
using HydroQuarticTrace=std::array<std::array<double,5>,6>;
struct HydroTraceAgreement { bool continuous=false;double relativeJump=0; };
HydroTraceAgreement hydroQuarticTraceAgreement(const HydroQuarticTrace&,const HydroQuarticTrace&);
}
class HydroCharacteristicSource;
// Actual original native points on one immutable endpoint source. Only the
// start locator and lazily used endpoint jets are shared across trace clocks.
class HydroCharacteristicNativeStarts
{
public:
    std::size_t pointCount() const;
    bool usesSource(const HydroCharacteristicSource&) const;
private:
    struct Data;
    explicit HydroCharacteristicNativeStarts(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroCharacteristicSource;
};
// Static search index of the SAME actual reconstruction pieces, prepared once.
class HydroCharacteristicGeometry
{
public:
    static HydroCharacteristicGeometry prepare(const QuadtreeMesh&, HydroAdvectorReconstruction);
    const HydroAdvectorReconstruction& reconstruction() const noexcept;
    std::vector<int> locatePieces(Vector2D) const;
    bool sameIdentity(const HydroCharacteristicGeometry&) const noexcept;
    void validateFor(const QuadtreeMesh&) const;
    HydroCharacteristicSource linear(const QuadtreeMesh&, HydroSpatialAdvector previous,
        HydroSpatialAdvector arrival) const;
private:
    struct Data;
    explicit HydroCharacteristicGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroCharacteristicSource;
};
struct HydroCharacteristicConfig
{
    double maxStep=0; // Explicit positive time budget, independent of the point.
    int maxSubsteps=4096, maxEvents=256;
};
struct HydroCharacteristicJet
{
    LevelSetBoundaryTrace::VelocityJet spatial;
    Vector2D timeDerivative{0,0};
};
struct HydroCharacteristicBoundary
{
    int normalAxis=0, outwardSign=0, piece=-1;
};
struct HydroCharacteristicContinuityReport
{
    std::size_t connections=0, transparentConnections=0;
    double maxTransparentRelativeJump=0;
};
class HydroCharacteristicResult
{
public:
    Vector2D point() const noexcept { return point_; }
    double time() const noexcept { return time_; }
    const std::optional<LevelSetBoundaryTrace::Jacobian>& jacobian() const noexcept { return derivative_; }
    const std::optional<HydroCharacteristicBoundary>& boundary() const noexcept { return boundary_; }
    bool oneSided() const noexcept { return oneSided_; }
    int substeps() const noexcept { return substeps_; }
    int events() const noexcept { return events_; }
    int splitEvents() const noexcept { return splitEvents_; }
    int initialPiece() const noexcept { return initialPiece_; }
    int finalPiece() const noexcept { return finalPiece_; }
private:
    Vector2D point_{0,0};double time_=0;
    std::optional<LevelSetBoundaryTrace::Jacobian> derivative_;
    std::optional<HydroCharacteristicBoundary> boundary_;
    bool oneSided_=false;
    int substeps_=0,events_=0,splitEvents_=0,initialPiece_=-1,finalPiece_=-1;
    friend class HydroCharacteristicSource;
};
// Immutable actual endpoint source. No arbitrary callback, geometry revision
// fabrication, field smoothing, accepted state or material inventory.
class HydroCharacteristicSource
{
public:
    const HydroCharacteristicGeometry& geometry() const noexcept;
    const HydroSpatialAdvector& previous() const noexcept;
    const HydroSpatialAdvector& arrival() const noexcept;
    const HydroCharacteristicContinuityReport& continuityReport() const noexcept;
    bool sameIdentity(const HydroCharacteristicSource&) const noexcept;
    void validateFor(const QuadtreeMesh&, HydroAdvectorIdentity previous,
        HydroAdvectorIdentity arrival) const;
    HydroCharacteristicJet sample(int piece, Vector2D, double physicalTime) const;
    HydroCharacteristicNativeStarts prepareNativeStarts(const QuadtreeMesh&,
        const Interface::HermiteGalsPublication&) const;
    // Derivative is of the SAME order-three piecewise Picard interior map
    // and cubic dense event map, at fixed end times. Quadrature nodes, widths
    // and crossing times are differentiated with the actual update.
    // A boundary hit has no interior derivative. Starting/ending on a piece
    // boundary is explicitly one-sided, not a globally differentiable view.
    HydroCharacteristicResult traceBackward(Vector2D arrivalPoint, double arrivalTime,
        double departureTime, HydroCharacteristicConfig, int initialPiece=-1) const;
    HydroCharacteristicResult traceBackward(Vector2D arrivalPoint, double arrivalTime,
        double departureTime, HydroCharacteristicConfig, const HydroCharacteristicNativeStarts&,
        int initialPiece=-1) const;
    // Same actual source and order-three piecewise/event kernel, progressing
    // from departureTime to arrivalTime within the immutable endpoint clocks.
    HydroCharacteristicResult traceForward(Vector2D departurePoint, double departureTime,
        double arrivalTime, HydroCharacteristicConfig, int initialPiece=-1) const;
private:
    HydroCharacteristicResult traceImpl(Vector2D startPoint, double startTime,
        double endTime, HydroCharacteristicConfig, const HydroCharacteristicNativeStarts*,
        int initialPiece, bool forward) const;
    struct Data;
    explicit HydroCharacteristicSource(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class HydroCharacteristicGeometry;
};
}
