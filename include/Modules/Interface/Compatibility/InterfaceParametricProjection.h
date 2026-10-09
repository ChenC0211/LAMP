#pragma once
// S08_FROZEN_BASELINE_BEGIN
// Frozen experimental baseline: parametricGals global curve/H2 projection.
// Retained only for explicit A/B comparisons and legacy case regression.
// No V2 features or implicit fallback: V2 owns a first-jet Hermite skeleton
// and independent cell phase volume, not this global parametric curve state.
// Scope and source identities: docs/task/S08_Frozen_Baselines.json.
// S08_FROZEN_BASELINE_END
#include "Core/Geometry/ParametricBoundary.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include <memory>
#include <stdexcept>

// One oriented free boundary in a conforming global Hermite H2 space: both the
// affine origin and correction share node jets. This is an
// explicit parametric representation, not a second nodal level-set snapshot.
// XY loops, RZ loops/axis meridians and fixed-rim supplied meridians share
// the same projection. Essential endpoint conditions belong to that space.
// The same H2/volume projection consumes uniform or actual mixed leaf cells.
// Dynamic inventory transfer and topology change are separate operations.
namespace InterfaceParametricProjection {
// A trial geometry cannot meet the requested accuracy in the permitted
// material space. Coupled nonlinear solvers may reject that trial; the
// projection itself never changes its accuracy or resource limit.
class MaterialResolutionLimit : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
struct NeckSplit;
struct NeckAssessment;
struct NeckSplitOptions;
using Curve=Quantum::Geometry::Parametric::Curve;
using Boundary=Quantum::Geometry::Parametric::Boundary;
using Source=Quantum::Geometry::Parametric::PolynomialCurve<9>;
using Point=Quantum::Geometry::Parametric::Point;
using Arc=Quantum::Geometry::Parametric::Arc;
struct Grid {
    struct Rectangle {double left,right,bottom,top;};
    struct Face {double lower,upper;int negative=-1,positive=-1,nativeRow=-1;};
    struct Plane {double coordinate=0,capacity=0;int offset=0;std::vector<Face> faces;};
    // Immutable view/index of actual native leaves and shared face segments.
    // It introduces no auxiliary mesh or sampled geometry. Inventory slots
    // are explicit; a coarse face is represented by its native fine segments.
    struct Leaves {
        explicit Leaves(const QuadtreeMesh& mesh):locator(mesh){}
        CartesianLeafLocator locator;
        std::vector<Rectangle> cells;
        std::vector<int> slotByCell;
        std::array<std::vector<Plane>,2> planes;
    };
    int nx=0,ny=0;
    double left=0,bottom=0,width=1,height=1;
    bool axisymmetric=false;
    std::shared_ptr<const Leaves> leaves;
    bool sameLayout(const Grid& other) const {
        return nx==other.nx&&ny==other.ny&&left==other.left&&bottom==other.bottom
            &&width==other.width&&height==other.height&&axisymmetric==other.axisymmetric&&leaves==other.leaves;
    }
    std::size_t size() const {return leaves?leaves->cells.size():static_cast<std::size_t>(nx)*ny;}
    Rectangle rectangle(std::size_t slot) const {
        if(leaves)return leaves->cells.at(slot);
        const double hx=width/nx,hy=height/ny,x=left+(slot/ny)*hx,y=bottom+(slot%ny)*hy;
        return {x,x+hx,y,y+hy};
    }
    int locate(Point p) const;
    std::vector<int> overlapping(const Rectangle&) const;
    std::vector<double> lines(int axis) const;
    double capacity(std::size_t slot) const {
        if(leaves){const auto c=rectangle(slot);return (c.right-c.left)*(c.top-c.bottom)*(axisymmetric?std::acos(-1.)*(c.left+c.right):1.);}
        const double area=width*height/size();
        return axisymmetric?area*2*std::acos(-1.)*(left+(slot/ny+.5)*width/nx):area;
    }
    double totalCapacity() const {return width*height*(axisymmetric?2*std::acos(-1.)*(left+.5*width):1.);}
};
struct Options {
    int materialSegments=16;
    // Resource ceiling for the compact hierarchical approximation, separate
    // from H2/volume accuracy. Exhaustion rejects the whole candidate.
    int maximumMaterialSegments=8192;
    // Nested material partition supplied by a coupled outer solve. These
    // knots are retained; further refinement still obeys the same H2 and
    // inventory budgets. Empty preserves standalone adaptive reconstruction.
    std::vector<double> minimumMaterialKnots;
    double relativeH2=1e-3;
    double inventoryFractionTolerance=1e-10;
    // Declared numerical uncertainty of the independent FV data, distinct from
    // the acceptance ceiling. Geometry cannot spend the whole ceiling to lower
    // its H2 objective. This is a numerical budget, not a certified enclosure.
    double inventoryFractionUncertainty=1e-12;
    // Reserve for numerical clipping/integration, within the total budget.
    double measureFractionReserve=2e-12;
};
struct Report {
    std::size_t arcs=0,refinements=0,activeCellConstraints=0;
    std::size_t basisCoefficients=0,denseBasisCoefficients=0;
    double relativeH2=0,maximumInventoryFractionError=0;
    double totalGeometryInventoryDifference=0,maximumGeometryDisplacement=0;
    double minimumLongitudinalDerivative=1,positionJoin=0,tangentJoin=0;
    double stationarityResidual=0;
    std::size_t endpointIncidences=0;
    double maximumEndpointIncidenceWidth=0;
};
// An exited component owns no in-domain geometry. Its signed numerical FV
// remainder is retained, never clipped or reassigned to an active component.
// The time is the accepted endpoint of removal, not an interpolated exit time.
struct RetiredComponent {
    std::uint64_t id=0;
    double time=0;
    std::shared_ptr<const std::vector<double>> inventory;
};
struct MeshTransferReport {
    struct Overlap {int source,target;Grid::Rectangle box;double volume;};
    // Slots in the respective source/target grids. Shared by material and
    // momentum consumers, not an independently sampled donor geometry.
    std::vector<Overlap> overlaps;
    double inventoryBefore=0,inventoryAfter=0,inventoryBalance=0;
    double maximumFractionError=0,maximumInheritedFraction=0;
    std::size_t refinedSources=0,coarsenedTargets=0,retainedCells=0;
};
// Endpoint first-jet prediction, plus the area functional of its vector bubble.
// The final projection can move shared first jets through one continuous field;
// this helper also remains available to explicit arc-moment controls.
struct ArcPrediction {Arc free;Point areaNormal;double sourceAreaDifference=0;};
ArcPrediction predictArc(const Source& source,double begin,double end);

// Immutable paired geometry/inventory candidate. The given inventory is copied
// without clipping or reseeding. Preparation is functional: failure cannot
// mutate the caller's old geometry, inventory or clock. It is not a coupled
// physics commit, and supplies no claim that its caller advanced FV correctly.
class Candidate {
public:
    const Curve& geometry() const {return *geometryHandle();}
    const std::shared_ptr<const Curve>& geometryHandle() const {return boundary_->single();}
    const std::shared_ptr<const Boundary>& boundaryHandle() const noexcept {return boundary_;}
    std::size_t componentCount() const noexcept {return boundary_?boundary_->components().size():0;}
    Candidate component(std::size_t) const;
    const std::vector<double>& componentInventory(std::size_t i) const {
        if(componentCount()==1&&retired_.empty()){if(i)throw std::out_of_range("invalid material component");return inventory_;}return componentInventories_.at(i);
    }
    const std::vector<RetiredComponent>& retiredComponents() const noexcept {return retired_;}
    const Grid& grid() const noexcept {return grid_;}
    const std::vector<double>& inventory() const noexcept {return inventory_;}
    const Report& report() const noexcept {return report_;}
    double liquidArea(double left,double right,double bottom,double top) const;
    double liquidVolume(double left,double right,double bottom,double top) const;
    std::vector<double> geometricVolumes() const;
    Quantum::Geometry::Parametric::DistanceSample query(Point point) const;
private:
    Candidate()=default;
    std::shared_ptr<const Boundary> boundary_;
    std::vector<double> inventory_;
    // The public inventory sums active and retired ledgers. The single-active
    // optimization is valid only when there are no retired remainders.
    std::vector<std::vector<double>> componentInventories_;
    std::vector<Report> componentReports_;
    std::vector<RetiredComponent> retired_;
    Grid grid_;
    Report report_;
    friend Candidate prepare(const Source&,const Grid&,const std::vector<double>&,const Options&);
    friend Candidate initialize(Curve,const Grid&);
    friend Candidate combine(std::vector<Candidate>,const std::vector<std::uint64_t>&);
    friend Candidate combine(const Grid&,std::vector<Candidate>,const std::vector<std::uint64_t>&,std::vector<RetiredComponent>);
    friend Candidate transferMesh(const Candidate&,const Grid&,const Options&,MeshTransferReport&);
    friend NeckAssessment assessNeckSplit(const Candidate&,const NeckSplitOptions&);
};
Candidate initialize(Curve effectiveBoundary,const Grid& grid);
Candidate initialize(std::vector<Curve>,const Grid&);
Candidate combine(std::vector<Candidate>,const std::vector<std::uint64_t>& ids);
// Grid ownership is explicit even when the last active component has left.
Candidate combine(const Grid&,std::vector<Candidate>,const std::vector<std::uint64_t>& ids,
    std::vector<RetiredComponent> retired={});
// Reserve the frozen remainder's fraction error before sharing the original
// budget between active reconstructions. Removing geometry grants no budget.
Options activeOptions(const Candidate&,const Options&);
// Explicit absorbing-outlet operation, called only after the physical step's
// geometry/flux/energy acceptance. No fixed inlet component is removed, and
// enclosing the entire domain is not confused with absence of the phase.
Candidate retireExterior(const Candidate&,double acceptedTime,const Options&);
// Zero-time nested Cartesian remap: exact same immutable curve, independent
// inventory conserved on each donor. Coarsening sums inventories; refinement
// distributes the donor's geometric/inventory discrepancy in the volume norm.
// It never reseeds/clips inventory or changes geometry. Original fraction
// budgets still apply; an unresolved remap rejects the whole candidate.
Candidate transferMesh(const Candidate&,const Grid&,const Options&,MeshTransferReport&);
Candidate prepare(const Source& source,const Grid& grid,
    const std::vector<double>& independentInventory,const Options& options);
// Necessary geometric feasibility, before evaluating expensive conservative
// face fluxes. Uses the same unconstrained H2 minimizer as prepare(), returns
// its resolved nested space, and never constructs/substitutes FV inventory.
std::vector<double> resolveMaterialSpace(const Source&,const Grid&,const Options&);
} // namespace InterfaceParametricProjection
