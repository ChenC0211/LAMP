#pragma once
#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include <array>
#include <vector>

// Explicit complete-grid XY reconstruction from native FV face averages.
// Tensor natural cubic interpolation of cumulative face fluxes gives C1
// velocity, preserves every native face integral, and commutes with discrete
// divergence. Natural primitive boundary conditions represent the closed
// zero-normal/free-slip box. No projection, smoothing of NS unknowns, or
// alternate velocity/flux pair is hidden in this immutable snapshot.
class HydroSplineVelocity {
public:
    HydroSplineVelocity(const QuadtreeMesh&,const NativeLeafDomain&,
        const VEdgeField<double>&,const HEdgeField<double>&);
    ~HydroSplineVelocity();
    HydroSplineVelocity(HydroSplineVelocity&&) noexcept;
    HydroSplineVelocity& operator=(HydroSplineVelocity&&) noexcept;
    HydroSplineVelocity(const HydroSplineVelocity&)=delete;
    HydroSplineVelocity& operator=(const HydroSplineVelocity&)=delete;
    HydroFittedVelocitySampler::FirstJet firstJet(Vector2D) const;
    const HydroFittedFaceFlux& faceFlux(int nativeRow) const;
    Vector2D lower() const;
    Vector2D upper() const;
    bool sameLayout(const HydroSplineVelocity&) const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Transpose of the SAME cumulative-flux/natural-spline reconstruction. Surface
// quadrature accumulates local primitive coefficient moments, then one reverse
// tensor interpolation and reverse cumulative sum returns integrated forces on
// native velocity DOFs. Consequently u.F equals integral traction.R(u) without
// reconstructing another velocity or approximating global basis functions.
class HydroSplineVelocityAdjoint {
public:
    using Moment=std::array<long double,32>;
    HydroSplineVelocityAdjoint(const QuadtreeMesh&,const NativeLeafDomain&);
    ~HydroSplineVelocityAdjoint();
    Moment moment(int nativeCell,Vector2D unitCoordinates,Vector2D forceMeasure) const;
    void accumulate(int nativeCell,const Moment&);
    // Includes constrained wall rows as reaction loads. Their velocity is
    // prescribed zero by Hydro, so they perform no work in the closed box.
    std::vector<double> integratedForce() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
