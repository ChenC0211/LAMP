#pragma once
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/Detail/HydroNativeVelocityTrace.h"
#include <memory>
#include <functional>
class InterfaceState;

namespace Quantum::Hydrodynamics
{
class HydroConservativePressure;
// Integrated div(2*mu*D(u)) on exactly the pressure/transport physical duals.
// Geometry and traces are immutable; materials are supplied by the Hydro owner.
// Wall/axis normal rows are retained as reactions, not unconstrained DOFs.
class HydroConservativeStress
{
public:
    struct EnergyReport {double maximumMeasureFractionAdjustment=0;};
    // Integral of viscosity over a real strain-support rectangle. Includes
    // dV=dx*dy in XY, or dV=2*pi*r*dr*dz in RZ, with no strain factor.
    using ViscosityMeasure=std::function<double(double x0,double x1,double y0,double y1)>;
    // Common +K=E^T H_mu E on uniform XY/RZ native rows. Normal strains and
    // RZ hoop strain live in cells; shear lives at interior vertices.
    // Zero-normal/free-slip boundary; the caller constrains boundary normal
    // rows (including u_r=0 on the axis), retaining their reactions in K.
    static Eigen::SparseMatrix<double> assembleEnergy(const QuadtreeMesh&,const NativeLeafDomain&,
        const ViscosityMeasure&);
    // Same energy with the Hydro owner's cell viscosity as a piecewise
    // constant measure. Real cell/support intersections determine all strain
    // weights, including RZ shear supports; no second vertex mixture law.
    static Eigen::SparseMatrix<double> assembleCellEnergy(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::vector<double>& cellViscosity);
    // Positive integrated dissipation K=E^T H_mu E on uniform XY MAC rows.
    // Normal strains live in cells, shear at interior vertices; viscosity is
    // integrated on the same parameter geometry over their real rectangles.
    // Closed zero-normal/free-slip box. Returns +K, unlike assemble's div(tau).
    static Eigen::SparseMatrix<double> assembleParametricEnergy(const QuadtreeMesh&,const NativeLeafDomain&,
        const InterfaceState&,double liquidViscosity,double gasViscosity,EnergyReport* report=nullptr);
    static HydroConservativeStress build(const QuadtreeMesh&,const NativeLeafDomain&,
        const HydroConservativePressure&,const std::vector<HydroFittedDomain::Boundary>& boundaryByNative);
    Eigen::SparseMatrix<double> assemble(const QuadtreeMesh&,const NativeLeafDomain&,
        const std::vector<double>& cellViscosity,const std::vector<double>& vertexViscosity) const;
private:
    struct Connection {
        int negative=-1,positive=-1;
        double area=0;
        HydroFittedDetail::MomentumTrace strain;
        // Positive material interpolation on real cells or real vertices.
        HydroFittedDetail::MomentumTrace cellMaterial,vertexMaterial;
    };
    struct Hoop {
        int row=-1,owner=-1;
        double constant=0,slope=0;
        HydroFittedDetail::MomentumTrace derivative;
    };
    std::shared_ptr<const NativeMeshIdentity> identity_;
    bool axisymmetric_=false;
    int faces_=0,vertices_=0,cells_=0;
    std::vector<Connection> connections_;
    std::vector<Hoop> hoop_;
};
}
