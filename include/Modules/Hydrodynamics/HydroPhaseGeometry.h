#pragma once

#include "Modules/Hydrodynamics/HydroPhysicalDivergence.h"
#include "Modules/Interface/HermiteGalsPublication.h"
#include "Modules/Interface/GeometricVofPublication.h"
#include <memory>
#include <vector>

namespace LAMP::Interface { class HermiteGalsFirstJetSample; class HermiteGalsFirstJetLocalSample; class HermiteGalsFirstJetTraceSample; }

namespace LAMP::Hydrodynamics
{
class HydroSpatialAdvector;
class HydroMomentumGeometry;
class HydroAdvectorIdentity;
class HydroAdvectorReconstruction;
class HydroAdvectorTimeSlab;
class HydroAdvectorAtomSupport;
struct HydroPhaseMeasure
{
    long double liquid = 0, gas = 0;
};
enum class HydroPhaseFaceTraceKind { Ordinary, InteriorCoincident, BoundaryCoincident };
// Spatial geometry only. Cell availability belongs to the selected Hydro
// domain, even when the immutable publication has additional scalar support.
// A coincident face has no single strict negative-phi phase on the face.
struct HydroPhaseFaceTrace
{
    HydroPhaseFaceTraceKind kind = HydroPhaseFaceTraceKind::Ordinary;
    int negativeCell = -1, positiveCell = -1;
    std::size_t coincidentHalfFaces = 0; // FirstJet half partitions only; VOF: not applicable (0).
    HydroPhaseMeasure minus, plus; // Whole physical trace on each available side.
    bool minusAvailable() const noexcept { return negativeCell>=0; }
    bool plusAvailable() const noexcept { return positiveCell>=0; }
};
struct HydroPhaseGeometryReport
{
    long double maxRelativeFaceClosure = 0;
    long double maxRelativeCellClosure = 0;
    // VOF mixed cells use the original minor inventory normalization.
    long double maxRelativeInventoryMeasurement = 0;
    std::size_t cutCells = 0, splitFaces = 0;
    std::size_t coincidentFaces = 0, boundaryCoincidentFaces = 0;
};
struct HydroVofBox
{
    int nativeCell=-1;
    Geometry::PlanarCellCut::Rectangle rectangle;
};
// Derived measurements of ONE explicitly selected retained publication. No
// new material inventory, geometry owner, clock or accepted transaction.
class HydroPhaseGeometry
{
public:
    enum class Model { FirstJet, GeometricVof };
    static HydroPhaseGeometry build(const QuadtreeMesh&, HydroPhysicalDivergence,
        Interface::HermiteGalsPublication);
    // VOF metrics retain the source domain's exact native mesh binding.
    // Independently rebuilt/selected domains need an explicit transfer.
    static HydroPhaseGeometry build(const QuadtreeMesh&, HydroPhysicalDivergence,
        Interface::GeometricVofPublication);
    Model model() const noexcept;
    const HydroPhysicalDivergence& metrics() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    // Typed ports reject the other model; no FirstJet identity is fabricated.
    const Interface::HermiteGalsPublication& publication() const;
    const Interface::GeometricVofPublication& vofPublication() const;
    Interface::GeometricVofIdentity vofIdentity() const;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    const std::vector<HydroPhaseMeasure>& cellVolumes() const noexcept; // Actual two-phase geometry, pressure DOFs
    // Ordinary strict phase; coincident: symmetric interior spatial trace or
    // available inside-domain boundary trace. Never a causal material rate.
    const std::vector<HydroPhaseMeasure>& faceMeasures() const noexcept; // velocity DOFs
    const std::vector<HydroPhaseFaceTrace>& faceTraces() const noexcept;
    const std::vector<double>& liquidFractions() const noexcept; // VOF: original owner alpha, bitwise.
    const HydroPhaseGeometryReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, Interface::HermiteGalsIdentity,
        Interface::HermiteGalsTimeLayer) const;
    void validateVofFor(const QuadtreeMesh&, Interface::GeometricVofIdentity,
        Interface::GeometricVofTimeLayer) const;
    // One mesh/source authentication for the whole batch. Each query lies
    // strictly within its selected original native leaf; both actual fan
    // phases are integrated independently from original points and common frame.
    std::vector<HydroPhaseMeasure> measureVofBoxes(const QuadtreeMesh&,
        const std::vector<HydroVofBox>&) const;
    // Ordinary meridional dA of the SAME actual LocalPlic plane, in XY or
    // RZ. Both phases are cut independently; no 2*pi*r weight, inventory
    // remap, radius division, or cached physical whole-cell volume is used.
    std::vector<HydroPhaseMeasure> measureVofMeridionalBoxes(const QuadtreeMesh&,
        const std::vector<HydroVofBox>&) const;
private:
    std::vector<HydroPhaseMeasure> measureVofBoxesWithMeasure(const QuadtreeMesh&,
        const std::vector<HydroVofBox>&,bool meridional) const;
    struct Data;
    explicit HydroPhaseGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};

struct HydroSpatialPhaseReport
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
class HydroPhaseRate
{
public:
    static HydroPhaseRate prepareSpatial(const QuadtreeMesh&, HydroPhaseGeometry,
        const HydroSpatialAdvector&);
    // Internal temporal scalar rates only, retaining the actual Q4 owner and
    // clock. No publication or derived phase geometry is constructed.
    static std::vector<HydroPhaseMeasure> sampleAtomRates(const QuadtreeMesh&,
        const Interface::HermiteGalsFirstJetSample&, const HydroSpatialAdvector&);
    // Exact closed-box scalar support of requested actual atoms, including
    // coincident sides and endpoint candidates. IDs must be unique and valid.
    static std::vector<int> requiredNativeCells(const QuadtreeMesh&,
        const HydroAdvectorReconstruction&, const std::vector<int>& atomIds);
    // Results follow atomIds order. Missing support rejects even for zero flux.
    static std::vector<HydroPhaseMeasure> sampleAtomRates(const QuadtreeMesh&,
        const Interface::HermiteGalsFirstJetLocalSample&, const HydroSpatialAdvector&,
        const std::vector<int>& atomIds);
    // Retained-source local query: slab preparation already authenticated the
    // mesh and endpoint states. This overload never reads a mutable mesh.
    static std::vector<HydroPhaseMeasure> sampleAtomRates(
        const Interface::HermiteGalsFirstJetLocalSample&,const HydroAdvectorTimeSlab&,
        const std::vector<int>& atomIds);
    // Actual requested straight traces only; the slab retains zero modal
    // endpoint coefficients. No full scalar/state identity is manufactured.
    static std::vector<HydroPhaseMeasure> sampleAtomRates(
        const Interface::HermiteGalsFirstJetTraceSample&,const HydroAdvectorTimeSlab&,
        const std::vector<int>& atomIds);
    // Consume the SAME authenticated support prepared by the FirstJet caller.
    static std::vector<HydroPhaseMeasure> sampleAtomRates(
        const Interface::HermiteGalsFirstJetLocalSample&,const HydroAdvectorTimeSlab&,
        const HydroAdvectorAtomSupport&);
    static std::vector<HydroPhaseMeasure> sampleAtomRates(
        const Interface::HermiteGalsFirstJetTraceSample&,const HydroAdvectorTimeSlab&,
        const HydroAdvectorAtomSupport&);
    const HydroPhaseGeometry& geometry() const noexcept;
    const std::vector<HydroPhaseMeasure>& faceRates() const noexcept;
    const std::vector<double>& totalFaceRates() const noexcept;
    const std::vector<HydroPhaseMeasure>& cellOutwardRates() const noexcept;
    bool isSpatial() const noexcept;
    // Original immutable spatial owner; isSpatial remains true for compatibility.
    const HydroSpatialAdvector& spatialAdvector() const;
    const HydroMomentumGeometry& momentumGeometry() const;
    HydroAdvectorIdentity advectorIdentity() const;
    const std::vector<HydroPhaseMeasure>& atomRates() const;
    const std::vector<HydroPhaseMeasure>& dualFaceRates() const;
    const HydroSpatialPhaseReport& spatialReport() const;
private:
    struct Data;
    explicit HydroPhaseRate(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
