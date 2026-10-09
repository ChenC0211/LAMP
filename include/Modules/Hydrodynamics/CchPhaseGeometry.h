#pragma once

#include "Modules/Hydrodynamics/CchPhysicalDivergence.h"
#include "Modules/Interface/HermiteGalsPublication.h"
#include <memory>
#include <vector>

namespace Quantum::Hydrodynamics
{
class CchSpatialAdvector;
class CchMomentumGeometry;
class CchAdvectorIdentity;
struct CchPhaseMeasure
{
    long double liquid = 0, gas = 0;
};
enum class CchPhaseFaceTraceKind { Ordinary, InteriorCoincident, BoundaryCoincident };
// Spatial geometry only. Cell availability belongs to the selected Hydro
// domain, even when the immutable publication has additional scalar support.
// A coincident face has no single strict negative-phi phase on the face.
struct CchPhaseFaceTrace
{
    CchPhaseFaceTraceKind kind = CchPhaseFaceTraceKind::Ordinary;
    int negativeCell = -1, positiveCell = -1;
    std::size_t coincidentHalfFaces = 0;
    CchPhaseMeasure minus, plus; // Whole physical trace on each available side.
    bool minusAvailable() const noexcept { return negativeCell>=0; }
    bool plusAvailable() const noexcept { return positiveCell>=0; }
};
struct CchPhaseGeometryReport
{
    long double maxRelativeFaceClosure = 0;
    std::size_t cutCells = 0, splitFaces = 0;
    std::size_t coincidentFaces = 0, boundaryCoincidentFaces = 0;
};
// Derived measurements of ONE retained FirstJet publication. No material
// inventory, scalar copy, clock, boundary condition or accepted transaction.
class CchPhaseGeometry
{
public:
    static CchPhaseGeometry build(const QuadtreeMesh&, CchPhysicalDivergence,
        Interface::HermiteGalsPublication);
    const CchPhysicalDivergence& metrics() const noexcept;
    const CchDofLayout& layout() const noexcept;
    const Interface::HermiteGalsPublication& publication() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    const std::vector<CchPhaseMeasure>& cellVolumes() const noexcept; // pressure DOFs
    // Ordinary strict phase; coincident: symmetric interior spatial trace or
    // available inside-domain boundary trace. Never a causal material rate.
    const std::vector<CchPhaseMeasure>& faceMeasures() const noexcept; // velocity DOFs
    const std::vector<CchPhaseFaceTrace>& faceTraces() const noexcept;
    const std::vector<double>& liquidFractions() const noexcept;
    const CchPhaseGeometryReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, Interface::HermiteGalsIdentity,
        Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit CchPhaseGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

struct CchSpatialPhaseReport
{
    std::size_t atoms=0, sourceQueries=0, ordinaryRanges=0, coincidentRanges=0;
    std::size_t events=0, refinements=0, tangentUnidentified=0;
    std::size_t tangentialAtoms=0;
    long double maxAbsTangentialRate=0;
    long double positiveIntegrationErrorEstimate=0;
    long double maxRelativeRefinementChange=0, maxRelativeAtomClosure=0;
};
// Spatial factory: the SAME actual compatible advector and FirstJet source.
// Instantaneous physical phase integrals, not dt-integrated transport or a
// conservative update of independently transported GALS geometry. Coincident
// traces consume causal inflow sides of the actual normal polynomial.
class CchPhaseRate
{
public:
    static CchPhaseRate prepareSpatial(const QuadtreeMesh&, CchPhaseGeometry,
        const CchSpatialAdvector&);
    const CchPhaseGeometry& geometry() const noexcept;
    const std::vector<CchPhaseMeasure>& faceRates() const noexcept;
    const std::vector<double>& totalFaceRates() const noexcept;
    const std::vector<CchPhaseMeasure>& cellOutwardRates() const noexcept;
    bool isSpatial() const noexcept;
    // Original immutable spatial owner; isSpatial remains true for compatibility.
    const CchSpatialAdvector& spatialAdvector() const;
    const CchMomentumGeometry& momentumGeometry() const;
    CchAdvectorIdentity advectorIdentity() const;
    const std::vector<CchPhaseMeasure>& atomRates() const;
    const std::vector<CchPhaseMeasure>& dualFaceRates() const;
    const CchSpatialPhaseReport& spatialReport() const;
private:
    struct Data;
    explicit CchPhaseRate(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
