#pragma once

#include "Modules/Interface/HermiteGalsBubble.h"
#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace Quantum::Interface
{
class HermiteGalsMeshTransfer;
struct HermiteGalsPublicationConfig
{
    HermiteGalsPhaseVolumeConfig volume;
    double relativeInventoryTolerance = 1e-10; // |raw geometric M - conserved M| / cell capacity
    double relativeInterfaceResidualTolerance = 1e-9; // frame residual / min(cell width)
    double unitNormalTolerance = 1e-12; // |norm(n)-1|
};
// Local Q4 center values obtained from the same old scalar/actual trace as
// the first jet. These are used directly; there is no volume correction.
struct HermiteGalsCenterPrediction
{
    HermiteGalsIdentity primaryIdentity;
    std::vector<double> values;
};
struct HermiteGalsSurfaceSample
{
    HermiteGalsIdentity geometryIdentity;
    HermiteGalsTimeLayer timeLayer;
    LevelSetHermiteGeometry::Frame frame; // raw same-scalar frame; curvature is meridional
    double azimuthalCurvature = 0;
    double curvature = 0; // physical +div(n): XY meridional, RZ meridional + azimuthal
};
struct HermiteGalsSurfaceQuadrature
{
    HermiteGalsSurfaceSample sample;
    double weight = 0, physicalWeight = 0; // ds, and ds / 2*pi*r*ds
};
struct HermiteGalsPublicationReport
{
    std::size_t components = 0, graphs = 0, quadratureNodes = 0, halfFaces = 0, axisEndpoints = 0;
    long double maxRelativeInventoryResidual = 0; // NaN for FirstJet (no independent M)
    double maxRelativeInterfaceResidual = 0, maxUnitNormalDefect = 0, maxAbsCurvature = 0;
};
class HermiteGalsPublication
{
public:
    const HermiteGalsState& state() const;
    const HermiteGalsPhaseVolumeResult& volumes() const; // geometric measurement; never fed back as state M
    const NativeLeafDomain& domain() const;
    const LevelSetInterfacePanels::NativeComponents& topology() const;
    const InterfacePhaseGeometry::MappedFacePartitions& facePartitions() const;
    const std::vector<HermiteGalsSurfaceQuadrature>& quadrature(std::size_t component, std::size_t graph) const;
    HermiteGalsSurfaceSample sample(std::size_t component, std::size_t graph, double fraction) const;
    const HermiteGalsPublicationReport& report() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsPublication(std::shared_ptr<const Data>);
    const Data& checked() const;
    static HermiteGalsPublication publish(const QuadtreeMesh&, HermiteGalsCandidate&&,
        const HermiteGalsBubbleResult*, const HermiteGalsAmrClosureResult*,
        const HermiteGalsCenterPrediction*, const HermiteGalsPublicationConfig&);
    static HermiteGalsPublication publishScalar(const QuadtreeMesh&, HermiteGalsCandidate&&,
        const HermiteGalsClosure&, LevelSetHermiteGeometry::Snapshot&&, const HermiteGalsPublicationConfig&);
    friend class HermiteGalsMeshTransfer;
    friend HermiteGalsPublication publishHermiteGalsGeometry(const QuadtreeMesh&,
        HermiteGalsCandidate&&, const HermiteGalsBubbleResult&, const HermiteGalsPublicationConfig&);
    friend HermiteGalsPublication publishHermiteGalsFirstJetGeometry(const QuadtreeMesh&,
        HermiteGalsCandidate&&, const HermiteGalsAmrClosureResult&, const HermiteGalsCenterPrediction&,
        const HermiteGalsPublicationConfig&);
};
// Builds/checks ONE snapshot using the I10 zero-center cardinal gauge.
// Failure preserves the argument; success consumes it, preserving its identity
// in the returned immutable state. This is NOT an App/Coupling accepted commit.
HermiteGalsPublication publishHermiteGalsGeometry(const QuadtreeMesh& mesh,
    HermiteGalsCandidate&& candidate, const HermiteGalsBubbleResult& bubble,
    const HermiteGalsPublicationConfig& config = {});
// Explicit nonconservative publication. Requires FirstJet state with no M;
// all measured volumes are diagnostics from the unique returned snapshot.
// Same transactional ownership and geometric checks as the historical port.
HermiteGalsPublication publishHermiteGalsFirstJetGeometry(const QuadtreeMesh& mesh,
    HermiteGalsCandidate&& candidate, const HermiteGalsAmrClosureResult& mixed,
    const HermiteGalsCenterPrediction& centers, const HermiteGalsPublicationConfig& config = {});
}
