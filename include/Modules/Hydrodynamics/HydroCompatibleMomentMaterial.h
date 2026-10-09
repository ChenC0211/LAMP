#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Core/Geometry/ParametricBoundary.h"
#include <optional>

namespace Quantum::Hydrodynamics {
// Frozen material on one immutable effective parametric boundary, or an
// explicitly prescribed single-phase viscosity field with no interface.
// This is a geometric integration consumer, NOT a moving-phase FV inventory.
// XY and meridional RZ use the same curve, mesh identity and physical measure
// for density, viscosity and kinetic mass. Properties remain fixed in time.
class HydroCompatibleMomentMaterial {
public:
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    using Curve=Geometry::Parametric::Curve;
    struct Properties { double liquidDensity=1,gasDensity=1,liquidViscosity=0,gasViscosity=0; };
    struct Cell {
        // Integral of (1,X,Y)*(1,X,Y)^T with each indicated measure.
        Eigen::Matrix3d volume,liquid,density,viscosity;
        std::vector<Eigen::Matrix3d> componentLiquid;
    };
    static HydroCompatibleMomentMaterial buildFrozen(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,std::shared_ptr<const Curve>,double geometryTime,Properties);
    static HydroCompatibleMomentMaterial buildFrozen(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,std::shared_ptr<const Geometry::Parametric::Boundary>,double geometryTime,Properties);
    // Physical viscosity, not a fit inferred from six old material moments.
    // Each matrix holds the tensor quadratic Bernstein coefficients on the
    // cell's [0,1]^2 chart (row=x/r degree, column=y/z degree). Nonnegative
    // coefficients define a nonnegative measure, including exact zero mu.
    using ViscosityField=std::vector<Eigen::Matrix3d>;
    static HydroCompatibleMomentMaterial buildSinglePhase(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,double density,const ViscosityField&,double time=0);
    static HydroCompatibleMomentMaterial buildSinglePhase(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroCompatibleMomentSpace&,double density,const std::vector<double>& cellViscosity,double time=0);
    void validate(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&) const;
    const std::vector<Cell>& cells() const noexcept{return cells_;}
    const Matrix& kineticMass() const noexcept{return mass_;}
    const std::shared_ptr<const Curve>& geometryHandle() const {return boundary_->single();}
    const std::shared_ptr<const Geometry::Parametric::Boundary>& boundaryHandle() const noexcept{return boundary_;}
    double geometryTime() const noexcept{return time_;}
    // The same immutable, fixed spatial law at a new evaluation time.
    // Does not advect geometry or synthesize a material flux.
    HydroCompatibleMomentMaterial atTime(double time) const;
    Properties properties() const noexcept{return properties_;}
    // A prescribed field replaces the phase-constant viscosities in Properties.
    const std::shared_ptr<const ViscosityField>& prescribedViscosity() const noexcept{return viscosity_;}
    double viscosityUpperBound() const noexcept;
    double viscosityAt(int cell,Vector2D unitCoordinate,bool liquid) const;
    std::optional<double> uniformDensity() const noexcept{return uniformDensity_;}
    std::vector<Eigen::Matrix3d> viscosityMoments() const;
private:
    std::shared_ptr<const NativeMeshIdentity> identity_;
    std::shared_ptr<const Geometry::Parametric::Boundary> boundary_;
    double time_=0;
    Properties properties_;
    std::optional<double> uniformDensity_;
    std::shared_ptr<const ViscosityField> viscosity_;
    double viscosityBound_=0;
    std::vector<Cell> cells_;
    Matrix mass_;
};
}
