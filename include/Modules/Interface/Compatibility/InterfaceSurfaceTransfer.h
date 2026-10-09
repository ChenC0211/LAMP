#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include <functional>

// Immutable geometric transfer plan for resolved motion without topology
// events on one fixed native mesh. Geometry supplies positive overlap weights;
// consuming modules retain ownership of their integrated material fields.
class InterfaceSurfaceTransfer
{
public:
    struct Boundary
    {
        LevelSetHermiteGeometry::Frame sourceFrame;
        Vector2D physicalPoint,physicalDisplacement,numericalDisplacement;
        LevelSetInterfacePanels::LocalProjection target;
        // Includes predictor interpolation and synchronization relocation,
        // not a second physical trajectory or a separately advected surface.
    };
    struct Overlap {std::size_t targetPanel=0;long double measure=0,share=0;};
    struct Row {long double mappedMeasure=0;std::vector<Overlap> overlaps;};
    struct Surface
    {
        std::uint64_t sourceId=0,targetId=0,liquidId=0;
        bool closed=false,axisymmetric=false;
        std::vector<double> sourceMeasures,targetMeasures;
        // Closed: N distinct old boundaries plus the same first boundary at N.
        // Open: N+1 distinct boundaries; ends are axes or explicitly retained
        // stationary wall contacts or sliding reflection-section endpoints.
        // No moving contact law is inferred.
        std::vector<Boundary> boundaries;
        // Piecewise constant density over each mapped old panel: the old
        // integrated quantity is distributed by new physical surface measure.
        // Conservative first-order remap, not higher-order charge advection.
        std::vector<Row> rows; // One per source panel; each row sums to one.
        double maxPhysicalDisplacement=0,maxNumericalDisplacement=0;
        // Observational: split GL8 measures are not an exactly additive metric.
        double maxRelativeMeasurePartitionError=0;
    };
    const std::vector<Surface>& surfaces() const noexcept {return surfaces_;}
    void validateFor(const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source,
        const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& target,const InterfaceStepTime&) const;
    const auto& predictor() const noexcept {return predictor_;}
private:
    InterfaceSurfaceTransfer()=default;
    static std::shared_ptr<const InterfaceSurfaceTransfer> create(const InterfaceState& source,
        const InterfaceState& target,const InterfaceStepTime&,
        const std::function<Vector2D(Vector2D)>& forwardMaterialPoint);
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_,predictor_,target_;
    InterfaceStepTime time_;
    std::vector<Surface> surfaces_;
    friend class InterfaceConservativeGeometry;
};
