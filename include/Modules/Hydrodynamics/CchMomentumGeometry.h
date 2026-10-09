#pragma once

#include "Modules/Hydrodynamics/CchPhaseGeometry.h"
#include <array>
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
struct CchMomentumOverlap
{
    int pressure = -1; // Same native pressure layout, never a virtual cell.
    Vector2D lower{0,0}, upper{0,0};
    long double volume = 0, cellFraction = 0;
};
struct CchMomentumVolume
{
    Vector2D lower{0,0}, upper{0,0};
    long double volume = 0;
    std::vector<CchMomentumOverlap> overlaps;
    std::array<std::vector<int>,4> sides; // -x/r,+x/r,-y/z,+y/z
};
struct CchMomentumSegment
{
    int negativeVelocity = -1, positiveVelocity = -1; // -1 outside selection.
    int normalAxis = 0;
    Vector2D begin{0,0}, end{0,0};
    long double area = 0;
    // Exactly one origin: an actual primal face or an active cell interior.
    int primalVelocity = -1, ownerPressure = -1;
};
struct CchMomentumGeometryReport
{
    long double maxRelativeCellPartition = 0, maxRelativeBoxClosure = 0;
    std::size_t overlaps = 0, sharedSegments = 0, boundarySegments = 0;
};
// Derived rectangles and common sides of the EXISTING native face velocities.
// Actual half-cell intersections; no global additional velocity unknowns,
// second primal divergence, pressure gauge, material state or boundary law.
class CchMomentumGeometry
{
public:
    static CchMomentumGeometry build(const QuadtreeMesh&, CchPhysicalDivergence);
    const CchPhysicalDivergence& metrics() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const std::vector<CchMomentumVolume>& volumes() const noexcept; // velocity DOFs
    const std::vector<CchMomentumSegment>& segments() const noexcept; // shared once
    const CchMomentumGeometryReport& report() const noexcept;
    bool sameIdentity(const CchMomentumGeometry&) const noexcept;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    explicit CchMomentumGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
struct CchMomentumMassReport
{
    std::array<CchPhaseMeasure,2> totalPhaseVolumes{}; // one partition per component
    std::array<long double,2> totalMass{};
    long double maxRelativeCellPartition = 0, maxRelativeLedgerDefect = 0;
    long double maxRelativeMomentEstimate = 0; // estimate, not a rigorous bound
    std::size_t cutOverlaps = 0;
};
// Read-only mass from the SAME retained FirstJet Q4 clipped into real dual
// rectangles. This is not an average of cell density or a conservative update.
class CchMomentumMass
{
public:
    static CchMomentumMass build(const QuadtreeMesh&, CchMomentumGeometry,
        CchPhaseGeometry, double liquidDensity, double gasDensity);
    const CchMomentumGeometry& geometry() const noexcept;
    const CchPhaseGeometry& phaseGeometry() const noexcept;
    const CchDofLayout& layout() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    double liquidDensity() const noexcept;
    double gasDensity() const noexcept;
    const std::vector<CchPhaseMeasure>& phaseVolumes() const noexcept;
    // Inner ordering matches geometry().volumes()[u].overlaps.
    const std::vector<std::vector<CchPhaseMeasure>>& overlapPhaseVolumes() const noexcept;
    const std::vector<double>& masses() const noexcept; // positive finite physical M_u
    const CchMomentumMassReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, Interface::HermiteGalsIdentity,
        Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit CchMomentumMass(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

class CchAdvectorReconstruction;
class CchSpatialAdvector;
struct CchKineticHodgePoint
{
    int velocityDual=-1, pressure=-1;
    Vector2D point{0,0};
    double physicalWeight=0;
    bool liquid=false;
};
struct CchKineticHodgeReport
{
    std::size_t markedDuals=0, points=0, sourceWeights=0, nonZeros=0,tangentialModes=0;
    long double maxRelativePhaseVolumeDifference=0, maxRelativeMomentEstimate=0;
    double symmetryDefect=0;
};
struct CchKineticHodgeEnergyChange
{
    long double beginEnergy=0, endEnergy=0, energyChange=0;
    long double midpointInertialWork=0, metricVariation=0;
    long double endpointQuadratureRemainder=0, identityResidual=0;
};
// Local kinetic metric of the supplied SAME prepared reconstruction, integrated
// on the authoritative phase-cut native duals. This assembly does not upgrade
// a C1 reconstruction into a phase-jump field or certify its physical accuracy.
// No velocity/geometry state, mass repair, clock owner or solver fallback.
class CchKineticHodge
{
public:
    using Matrix=Eigen::SparseMatrix<double,Eigen::RowMajor,Eigen::Index>;
    static CchKineticHodge build(const QuadtreeMesh&,CchMomentumMass,
        CchAdvectorReconstruction);
    const CchMomentumMass& materialMass() const noexcept;
    const CchAdvectorReconstruction& reconstruction() const noexcept;
    const CchDofLayout& layout() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    const Matrix& matrix() const noexcept;
    // Ordering [native normal U, independent tangential q], with ALL cross blocks.
    const Matrix& augmentedMatrix() const noexcept;
    const std::vector<bool>& markedDuals() const noexcept;
    const std::vector<CchKineticHodgePoint>& points() const noexcept;
    // One complete physical-volume partition, for same-source phase gradients.
    const std::vector<CchKineticHodgePoint>& strainQuadrature() const noexcept;
    const CchKineticHodgeReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&,Interface::HermiteGalsIdentity,
        Interface::HermiteGalsTimeLayer) const;
    std::vector<double> apply(const QuadtreeMesh&,const CchNativeValues&,
        Interface::HermiteGalsIdentity,Interface::HermiteGalsTimeLayer) const;
    long double kineticEnergy(const QuadtreeMesh&,const CchNativeValues&,
        Interface::HermiteGalsIdentity,Interface::HermiteGalsTimeLayer) const;
    std::vector<double> apply(const QuadtreeMesh&,const CchSpatialAdvector&) const;
    long double kineticEnergy(const QuadtreeMesh&,const CchSpatialAdvector&) const;
    // Force per physical volume at strainQuadrature(), tested by the SAME
    // full [U,q] velocity basis. No density multiplication or provider owner.
    std::vector<double> bodyForceLoad(const QuadtreeMesh&,const std::vector<Vector2D>&) const;
    // Exact endpoint algebra; the nonzero temporal remainder is retained.
    // This diagnostic does not replace the material F/R/uR energy terms.
    CchKineticHodgeEnergyChange energyChange(const QuadtreeMesh&,
        const CchKineticHodge& begin,const CchKineticHodge& end,
        const CchNativeValues& beginVelocity,const CchNativeValues& endVelocity) const;
    CchKineticHodgeEnergyChange energyChange(const QuadtreeMesh&,
        const CchKineticHodge& begin,const CchKineticHodge& end,
        const CchSpatialAdvector& beginVelocity,const CchSpatialAdvector& endVelocity) const;
private:
    struct Data;
    explicit CchKineticHodge(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
