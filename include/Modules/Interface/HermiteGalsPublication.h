#pragma once

#include "Modules/Interface/HermiteGalsBubble.h"
#include "Modules/Interface/HermiteGalsPhaseVolume.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cstddef>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace LAMP::Interface
{
class HermiteGalsMeshTransfer;
enum class HermiteGalsScalarRepresentation { SingleQuartic, C2MacroQuartic };
struct HermiteGalsPublicationConfig
{
    HermiteGalsGradientUpdate gradientUpdate = HermiteGalsGradientUpdate::CharacteristicPullback;
    HermiteGalsScalarRepresentation scalarRepresentation = HermiteGalsScalarRepresentation::SingleQuartic;
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
// One immutable scalar preparation for full FirstJet publication and
// intermediate transport. Verifies primary/mixed/actual-center identity and
// full clock before constructing the selected representation. C2 requires
// complete uniform XY support; unsupported scopes reject without fallback.
LevelSetHermiteGeometry::Snapshot prepareHermiteGalsFirstJetScalar(const QuadtreeMesh&,
    const HermiteGalsCandidate&, const HermiteGalsAmrClosureResult&,
    const HermiteGalsCenterPrediction&, const HermiteGalsPublicationConfig&);
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
// Passive caller elapsed times for this publication. The total also includes
// small unclassified setup/ownership costs; the four stages need not sum to it.
struct HermiteGalsPublicationTimings
{
    double privateSnapshot = 0, topologyCuts = 0, phaseVolumes = 0, surfaceCertification = 0, total = 0;
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
    // Read-only GL8 on an original graph-fraction interval. Retains each
    // subarc frame/chart and applies this publication's actual checks.
    // [0,1] copies the accepted cache; an empty interval returns no nodes.
    std::vector<HermiteGalsSurfaceQuadrature> quadrature(std::size_t component, std::size_t graph,
        double begin, double end) const;
    HermiteGalsSurfaceSample sample(std::size_t component, std::size_t graph, double fraction) const;
    const HermiteGalsPublicationReport& report() const;
    const HermiteGalsPublicationTimings& timings() const;
    // Validate the caller's actual mesh against the state exact identity.
    // The private immutable domain/scalar layout retains its construction
    // proof; public independent partition queries keep their full checks.
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
        const HermiteGalsClosure&, LevelSetHermiteGeometry::Snapshot&&, const HermiteGalsPublicationConfig&,
        double precedingPrivateSeconds = 0);
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

namespace LAMP::Electric { class ElectricMaxwellTrace; }

namespace LAMP::Interface
{
// An authenticated, immutable vector trace on the ONE retained publication.
// This view contains no force discretization or accept operation. Its only
// producer currently is the actual coeval electric Maxwell trace.
class HermiteGalsSurfaceVectorView
{
public:
    const HermiteGalsPublication& publication() const { return checked().publication; }
    HermiteGalsIdentity geometryIdentity() const { return publication().state().identity(); }
    HermiteGalsTimeLayer timeLayer() const { return publication().state().timeLayer(); }
    Vector2D sample(std::size_t component,std::size_t graph,const HermiteGalsSurfaceSample& q) const
    {
        const auto& d=checked();const auto& curves=d.publication.topology().curves;
        require(component<curves.size()&&graph<curves[component].graphs().size(),"original graph is unavailable");
        require(q.geometryIdentity==geometryIdentity()&&sameClock(q.timeLayer,timeLayer()),
            "surface sample lost publication/full clock");
        const auto& frame=q.frame;
        require(frame.nativeCellIndex==curves[component].graphs()[graph].nativeCellIndex,
            "surface sample has a different native graph cell");
        // The endpoint-relative chart is authoritative; a rounded global
        // point cannot be used to recover an unresolved local position.
        const auto chart=LevelSetHermiteGeometry::cellChartPoint(frame);
        const double fraction=curves[component].graphFraction(static_cast<int>(graph),chart);
        require(std::isfinite(fraction)&&fraction>=0&&fraction<=1,"surface sample is outside its original graph");
        require(std::isfinite(frame.point.x)&&std::isfinite(frame.point.y)&&std::isfinite(frame.normal.x)
            &&std::isfinite(frame.normal.y)&&std::isfinite(q.curvature),"nonfinite surface sample");
        const auto value=d.provider(component,graph,q);
        require(std::isfinite(value.x)&&std::isfinite(value.y),"nonfinite source vector trace");
        return value;
    }
    void validateFor(const QuadtreeMesh& mesh,const HermiteGalsPublication& actual,
        HermiteGalsTimeLayer expectedClock) const
    {
        publication().validateMesh(mesh);actual.validateMesh(mesh);(void)expectedClock.time();
        require(geometryIdentity()==actual.state().identity(),"different surface publication identity");
        require(sameClock(timeLayer(),actual.state().timeLayer())&&sameClock(timeLayer(),expectedClock),
            "different surface full clock");
    }
    bool sameIdentity(const HermiteGalsSurfaceVectorView& other) const noexcept { return data_==other.data_; }
private:
    struct Data
    {
        Data(HermiteGalsPublication p,
            std::function<Vector2D(std::size_t,std::size_t,const HermiteGalsSurfaceSample&)> f)
            :publication(std::move(p)),provider(std::move(f)) {}
        HermiteGalsPublication publication;
        std::function<Vector2D(std::size_t,std::size_t,const HermiteGalsSurfaceSample&)> provider;
    };
    static void require(bool value,const char* reason)
    { if(!value)throw std::invalid_argument(std::string("[HermiteGalsSurfaceVectorView] ")+reason); }
    static bool sameClock(const HermiteGalsTimeLayer& a,const HermiteGalsTimeLayer& b)
    { return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position; }
    const Data& checked() const
    { if(!data_)throw std::logic_error("[HermiteGalsSurfaceVectorView] moved-from view");return *data_; }
    HermiteGalsSurfaceVectorView(HermiteGalsPublication publication,
        std::function<Vector2D(std::size_t,std::size_t,const HermiteGalsSurfaceSample&)> provider)
    {
        require(static_cast<bool>(provider),"missing immutable vector source");
        (void)publication.state().timeLayer().time();
        data_=std::make_shared<const Data>(std::move(publication),std::move(provider));
    }
    std::shared_ptr<const Data> data_;
    friend class LAMP::Electric::ElectricMaxwellTrace;
};
}
