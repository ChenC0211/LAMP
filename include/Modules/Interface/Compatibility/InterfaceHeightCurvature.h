#pragma once
#include "Core/Base/Vector2D.h"
#include <cstddef>
#include <memory>
#include <limits>
#include <vector>

class QuadtreeMesh;
class NativeLeafDomain;
class InterfaceState;

// Read-only height-function geometry on ONE immutable conservative state.
// Borrowed mesh/domain/state must remain alive and unchanged. Rebuild after
// accepting a step or changing the mesh. No scalar Hessian or old curvature
// is substituted when a height stencil is unavailable. Use five exact column
// moments for quartic heights; avoid overdetermined curvature stencils whose
// flat-interface high-frequency symbol reverses the restoring sign.
// Non-axis stencils require every sampled supporting normal to have
// axial component >= 1/2; the regular even axis stencil retains its sign guard.
// At the axis, try an even implicit conic matching the same three annular
// volume means; retain the polynomial chart when no regular fit is available.
// This read-only estimate does not replace the geometry synchronization chart.
// On mixed meshes each attempted column rectangle must contain actual leaves
// at the query leaf's level. Coarse/fine gaps are not walls or virtual pure
// cells; they return MixedResolution. Physical inactive spans retain their
// existing boundary handling. No conservative C is interpolated to a grid.
class InterfaceHeightCurvature
{
public:
    enum class Status { Valid, InsufficientColumns, MissingPureEnds, MultipleTransitions,
        MixedBranches, PoorOrientation, SingularFit, GeometryMismatch, PoorMomentFit, MixedResolution };
    struct Sample
    {
        Status status=Status::InsufficientColumns;
        Vector2D point,normal;
        double meridional=std::numeric_limits<double>::quiet_NaN();
        double azimuthal=std::numeric_limits<double>::quiet_NaN();
        double curvature=std::numeric_limits<double>::quiet_NaN(); // +div(n), liquid -> gas
        double heightMismatch=std::numeric_limits<double>::quiet_NaN();
        // Retained fitted graph H(t), H'(t), H''(t), and physical support box.
        // Diagnostics can distinguish column transport error from derivative
        // error without reevaluating or substituting a different inventory.
        double height=std::numeric_limits<double>::quiet_NaN();
        double first=std::numeric_limits<double>::quiet_NaN();
        double second=std::numeric_limits<double>::quiet_NaN();
        Vector2D supportLower,supportUpper;
        // Geometrically unresolved minority inventory in the accepted column
        // support. It remains included in the raw volume moments, not deleted.
        std::size_t unresolvedColumnCells=0;
        long double unresolvedColumnVolume=0;
        double scaledRankRatio=0,scaledMomentResidual=0;
        int normalAxis=-1,columns=0,halfHeight=0;
        bool evenAxisFit=false,evenSectionFit=false,leastSquaresMoments=false,implicitAxisFit=false;
        // One-sided columns anchored to a declared physical wall endpoint;
        // no ghost inventory, contact-angle law or endpoint motion is inferred.
        bool oneSidedBoundary=false;
        Vector2D boundaryAnchor;
        // Ordered graph IDs from a contiguous, consistently oriented arc of
        // this surface. Not a whole-liquid ID used as a branch substitute.
        std::vector<int> branchGraphs;
        bool valid() const noexcept { return status==Status::Valid; }
    };
    InterfaceHeightCurvature(const QuadtreeMesh&,const NativeLeafDomain&,const InterfaceState&);
    ~InterfaceHeightCurvature();
    InterfaceHeightCurvature(const InterfaceHeightCurvature&)=delete;
    InterfaceHeightCurvature& operator=(const InterfaceHeightCurvature&)=delete;
    // Actual retained native graph coordinate, endpoints included. Evaluates
    // at that state's canonical interface, never an offset level surface.
    Sample evaluate(std::size_t surface,int graph,double fraction=.5) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
