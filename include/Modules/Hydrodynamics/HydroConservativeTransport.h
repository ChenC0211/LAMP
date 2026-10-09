#pragma once
#include "Modules/Hydrodynamics/Detail/HydroFittedMomentumDuals.h"
#include <functional>
#include <vector>
#include <memory>
#include <cstddef>

class InterfaceConservativePhase;
class InterfaceState;
struct InterfaceStepTime;
struct HydroStepTime;
enum class HydroForceTimeLevel;
class NativeMeshIdentity;
namespace Quantum::Geometry::Parametric {template<std::size_t> struct PolynomialCurve;}
template<class T> class CellField;

namespace Quantum::Hydrodynamics
{
// One accepted mass inventory step mapped onto the EXISTING staggered
// control volumes. All fluxes are time-integrated masses, not rates. This
// service does not advance Hydro or replace its pressure/viscous operators.
class HydroConservativeTransport
{
public:
    using Duals=HydroFittedDetail::MomentumDuals;
    struct MassFlux {long double positive=0,negative=0;};
    struct MomentumFlux
    {
        std::vector<long double> outward; // Native packed velocity rows.
        long double boundaryOutward[2]{}; // Separate component partitions.
    };
    // Independent incompressible single-phase source on a closed native mesh.
    // The supplied native velocities are FACE AVERAGES, frozen over this
    // interval; rho*dt*area*u supplies the same primal-to-dual flux map as an
    // Interface step. Cell mass stays rho*volume. Rejects nonzero wall flux
    // and a cell/dual mass defect above the common inventory-roundoff budget;
    // it does not project or repair the velocity.
    static HydroConservativeTransport singlePhase(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroStepTime&,double density,const std::vector<double>& nativeVelocity);
    bool isSinglePhase() const noexcept {return !singlePhaseVelocity_.empty();}
    void validateSinglePhaseMidpoint(const std::vector<double>& nativeVelocity) const;
    // Complete native XY/RZ leaves, including balanced coarse/fine faces.
    // Cell mass is rho_l*V_l+rho_g*(V-V_l).
    // Half-cell mass fractions include the RZ radial measure. Shared dual
    // fluxes commute with this mass mapping; gross counterflow is retained.
    static HydroConservativeTransport build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const InterfaceConservativePhase& previous,const InterfaceConservativePhase& next,
        const InterfaceStepTime& time,double liquidDensity,double gasDensity);
    // Same dual mass/flux mapping for the immutable parametric Interface step.
    // Reads its actual FV inventory and retained gross face volumes directly;
    // neither reseeds a legacy phase state nor runs a second phase limiter.
    static HydroConservativeTransport build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const InterfaceState& previous,const InterfaceState& next,
        const InterfaceStepTime& time,double liquidDensity,double gasDensity,
        const InterfaceState* midpoint=nullptr);
    struct MomentumRemap
    {
        // Unprojected candidate only: no pressure correction, boundary law,
        // accepted Hydro owner or physical clock is published by this service.
        std::vector<long double> mass,momentum;
        std::vector<double> velocity; // Existing native VEdge then HEdge rows.
        long double sourceMomentum[2]{},targetMomentum[2]{};
        long double sourceAbsoluteMomentum[2]{};
        long double sourceKineticEnergy[2]{},targetKineticEnergy[2]{};
        double maxRelativeMassClosure=0,maxRelativeMomentumClosure=0;
        std::size_t refinedParents=0,retainedParents=0;
    };
    struct ProjectedMomentumRemap
    {
        MomentumRemap unprojected;
        std::vector<double> velocity;
        long double momentum[2]{},kineticEnergy=0;
        // If the momentum-constrained mass projection exceeds the original
        // energy budget, contract only its zero-momentum solenoidal part
        // toward the minimum-energy feasible field. No global velocity scale.
        long double kineticEnergyBeforeConstraint=0,minimumKineticEnergy=0;
        double energyContraction=1;
        bool energyContractionApplied=false;
        double divergenceBefore=0,divergenceAfter=0,localContinuity=0;
        double momentumRelativeClosure=0,maxVelocityCorrection=0,linearRelativeResidual=0;
        int independentMomentumConstraints=0,projectionSolves=0;
        bool usedSuperLU=false;
    };
    struct StateRemapReport
    {
        // Each velocity layer uses its OWN source and transferred C masses.
        ProjectedMomentumRemap current,previous;
        long double sourcePressureIntegral=0,targetPressureIntegral=0;
    };
    // Zero-time refinement on the same native background lattice. The target
    // C MUST be the actual conservative transfer of source C. Refined parents
    // distribute each old half-piece mass share across their target mixture
    // masses; retained parents use physical half-piece intersections. This
    // positive first-order map preserves total component momentum and constant
    // velocity to measure roundoff, but does NOT guarantee divergence or wall
    // normal velocity. It cannot be accepted as an AMR NS state without a joint
    // constrained projection. No coarsening, momentum clipping or phase reseed.
    static MomentumRemap remapMomentum(const QuadtreeMesh& sourceMesh,const NativeLeafDomain& sourceDomain,
        const InterfaceConservativePhase& source,const std::vector<double>& sourceVelocity,
        const QuadtreeMesh& targetMesh,const NativeLeafDomain& targetDomain,
        const InterfaceConservativePhase& target,double liquidDensity,double gasDensity);
    const Duals& duals() const noexcept {return duals_;}
    const std::vector<long double>& massBefore() const noexcept {return before_;}
    const std::vector<long double>& massAfter() const noexcept {return after_;}
    const std::vector<MassFlux>& fluxes() const noexcept {return fluxes_;}
    double maxRelativeClosure() const noexcept {return maxRelativeClosure_;}
    std::size_t splitCellMappings() const noexcept {return splitCellMappings_;}
    // Bind the retained phase step to Hydro's current mesh, clock and source
    // mixture. Managed closed mixed leaves use the physical native NS strategy;
    // root-grid boundary patches retain their original operator. An unprojected
    // momentum remap alone is not an accepted Hydro state.
    // Validation precedes any Hydro field/history mutation.
    void validateFor(const QuadtreeMesh& mesh,bool axisymmetric,const HydroStepTime& time,
        const CellField<double>& sourceDensity) const;
    void validateParametricForce(const InterfaceState&,HydroForceTimeLevel) const;
    // First-order donor velocity for each GROSS direction, including temporal
    // reversal. A physical boundary callback is required only for inflow; it
    // receives the actual dual segment and transported velocity component.
    // Results are extensive momentum, with no additional AB2 extrapolation.
    MomentumFlux upwind(const std::vector<double>& nativeVelocity,
        const std::function<double(const Duals::Segment&,int)>& incoming={}) const;
    // Shared net mass flux times the arithmetic mean of adjacent velocities.
    // Closed uniform XY/RZ without swirl, resolved smooth baseline; no limiter.
    MomentumFlux centered(const std::vector<double>& nativeVelocity) const;
private:
    HydroConservativeTransport()=default;
    static HydroConservativeTransport fromCellMasses(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceStepTime&,const std::vector<long double>& before,const std::vector<long double>& after,
        const std::vector<MassFlux>& primal);
    Duals duals_;
    std::vector<long double> before_,after_;
    std::vector<MassFlux> fluxes_;
    int radialRows_=0;
    double maxRelativeClosure_=0;
    std::size_t splitCellMappings_=0;
    std::shared_ptr<const NativeMeshIdentity> identity_;
    std::shared_ptr<const Quantum::Geometry::Parametric::PolynomialCurve<5>> parameterBefore_,parameterAfter_,parameterMiddle_;
    std::vector<double> sourceDensity_;
    std::vector<double> singlePhaseVelocity_;
    bool axisymmetric_=false;
    bool uniformLeaves_=false;
    int step_=0;
    double begin_=0,dt_=0;
};
}
