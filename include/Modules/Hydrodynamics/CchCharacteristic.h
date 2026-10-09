#pragma once

#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include <optional>

namespace Quantum::Hydrodynamics
{
namespace Detail
{
// Six physical value/scaled-J traces, sampled at t=0,1/4,1/2,3/4,1.
// Complete for the actual piecewise Q4 field. Exposed internally for negative controls.
using CchQuarticTrace=std::array<std::array<double,5>,6>;
struct CchTraceAgreement { bool continuous=false;double relativeJump=0; };
CchTraceAgreement cchQuarticTraceAgreement(const CchQuarticTrace&,const CchQuarticTrace&);
}
class CchCharacteristicSource;
// Static search index of the SAME actual reconstruction pieces, prepared once.
class CchCharacteristicGeometry
{
public:
    static CchCharacteristicGeometry prepare(const QuadtreeMesh&, CchAdvectorReconstruction);
    const CchAdvectorReconstruction& reconstruction() const noexcept;
    std::vector<int> locatePieces(Vector2D) const;
    bool sameIdentity(const CchCharacteristicGeometry&) const noexcept;
    void validateFor(const QuadtreeMesh&) const;
    CchCharacteristicSource linear(const QuadtreeMesh&, CchSpatialAdvector previous,
        CchSpatialAdvector arrival) const;
private:
    struct Data;
    explicit CchCharacteristicGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class CchCharacteristicSource;
};
struct CchCharacteristicConfig
{
    double maxStep=0; // Explicit positive time budget, independent of the point.
    int maxSubsteps=4096, maxEvents=256;
};
struct CchCharacteristicJet
{
    LevelSetBoundaryTrace::VelocityJet spatial;
    Vector2D timeDerivative{0,0};
};
struct CchCharacteristicBoundary
{
    int normalAxis=0, outwardSign=0, piece=-1;
};
struct CchCharacteristicContinuityReport
{
    std::size_t connections=0, transparentConnections=0;
    double maxTransparentRelativeJump=0;
};
class CchCharacteristicResult
{
public:
    Vector2D point() const noexcept { return point_; }
    double time() const noexcept { return time_; }
    const std::optional<LevelSetBoundaryTrace::Jacobian>& jacobian() const noexcept { return derivative_; }
    const std::optional<CchCharacteristicBoundary>& boundary() const noexcept { return boundary_; }
    bool oneSided() const noexcept { return oneSided_; }
    int substeps() const noexcept { return substeps_; }
    int events() const noexcept { return events_; }
    int splitEvents() const noexcept { return splitEvents_; }
    int initialPiece() const noexcept { return initialPiece_; }
    int finalPiece() const noexcept { return finalPiece_; }
private:
    Vector2D point_{0,0};double time_=0;
    std::optional<LevelSetBoundaryTrace::Jacobian> derivative_;
    std::optional<CchCharacteristicBoundary> boundary_;
    bool oneSided_=false;
    int substeps_=0,events_=0,splitEvents_=0,initialPiece_=-1,finalPiece_=-1;
    friend class CchCharacteristicSource;
};
// Immutable actual endpoint source. No arbitrary callback, geometry revision
// fabrication, field smoothing, accepted state or material inventory.
class CchCharacteristicSource
{
public:
    const CchCharacteristicGeometry& geometry() const noexcept;
    const CchSpatialAdvector& previous() const noexcept;
    const CchSpatialAdvector& arrival() const noexcept;
    const CchCharacteristicContinuityReport& continuityReport() const noexcept;
    bool sameIdentity(const CchCharacteristicSource&) const noexcept;
    void validateFor(const QuadtreeMesh&, CchAdvectorIdentity previous,
        CchAdvectorIdentity arrival) const;
    CchCharacteristicJet sample(int piece, Vector2D, double physicalTime) const;
    // Derivative is of the SAME order-three piecewise Picard interior map
    // and cubic dense event map, at fixed end times. Quadrature nodes, widths
    // and crossing times are differentiated with the actual update.
    // A boundary hit has no interior derivative. Starting/ending on a piece
    // boundary is explicitly one-sided, not a globally differentiable view.
    CchCharacteristicResult traceBackward(Vector2D arrivalPoint, double arrivalTime,
        double departureTime, CchCharacteristicConfig, int initialPiece=-1) const;
private:
    struct Data;
    explicit CchCharacteristicSource(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class CchCharacteristicGeometry;
};
}
