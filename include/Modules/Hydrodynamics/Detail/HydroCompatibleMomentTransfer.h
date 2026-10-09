#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"

namespace Quantum::Hydrodynamics::MomentTransfer {
// Same finite-dimensional momentum functional as the material mass matrix.
// Constraints conserve its two components on each maximal old/new leaf patch
// and enforce target D and the old kinetic-energy upper bound, in one target
// kinetic-metric minimization. Incompatible mesh changes are rejected.
Eigen::VectorXd velocity(const QuadtreeMesh&,const NativeLeafDomain&,
    const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,const Eigen::VectorXd&,double,
    const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,const HydroCompatibleMomentMaterial&,
    const InterfaceParametricProjection::MeshTransferReport&,const std::vector<int>& sourceCells,const std::vector<int>& targetCells,
    HydroCompatibleMomentFlow::TransferLevel&,
    const HydroCompatibleMomentBoundary* sourceBoundary=nullptr,const HydroCompatibleMomentBoundary* targetBoundary=nullptr);
Eigen::VectorXd pressure(const QuadtreeMesh&,const Eigen::VectorXd&,const QuadtreeMesh&,const NativeLeafDomain&,
    const InterfaceParametricProjection::MeshTransferReport&,const std::vector<int>& sourceCells,const std::vector<int>& targetCells);
// Same constrained kinetic projection, on unchanged native cells. Retain
// each cell momentum and actual inhomogeneous boundary values. Surface energy
// change sets the kinetic upper bound through total-energy nonincrease.
enum class Feasibility {Admissible,Momentum,Energy};
// Complete-field counterpart: closest source in the target physical mass
// metric, with patch momentum and original boundary data. The only velocity
// unknown is q; native moments remain derived observations.
std::optional<HydroCompatibleMomentSpline::NativeVelocity> nativeVelocity(
    const HydroCompatibleMomentSpline::NativeVelocity& source,const HydroCompatibleMomentSpline::NativeMaterialForms& before,
    const HydroCompatibleMomentSpline::NativeVelocity& target,const HydroCompatibleMomentSpline::NativeProjectionForms& projection,
    const HydroCompatibleMomentSpace&,const HydroCompatibleMomentBoundary&,
    const InterfaceParametricProjection::MeshTransferReport&,const std::vector<int>& sourceCells,const std::vector<int>& targetCells,
    double surfaceEnergyChange,HydroCompatibleMomentFlow::TransferLevel&,Feasibility&);
Eigen::VectorXd eventVelocity(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,
    const HydroCompatibleMomentMaterial& before,const HydroCompatibleMomentMaterial& after,
    const Eigen::VectorXd&,double time,const HydroCompatibleMomentBoundary&,double surfaceEnergyChange,
    const std::vector<int>& cells,HydroCompatibleMomentFlow::TransferLevel&,Feasibility&);
}
