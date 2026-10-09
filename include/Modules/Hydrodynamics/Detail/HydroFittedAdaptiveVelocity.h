#pragma once

#include "Modules/Hydrodynamics/HydroFittedVelocity.h"
#include <memory>

namespace HydroFittedDetail
{
// Internal implementation of the explicit shared-vertex-trace constructor.
// No independent Hydro state/time step or numerical output is owned here.
class AdaptiveVelocity
{
  public:
    AdaptiveVelocity(const QuadtreeMesh &, const HydroFittedDomain &, const VEdgeField<double> &,
                     const HEdgeField<double> &,
                     const std::vector<HydroFittedVelocitySampler::VertexTrace> &);
    AdaptiveVelocity(const QuadtreeMesh &, const NativeLeafDomain &, const VEdgeField<double> &,
                     const HEdgeField<double> &,
                     const std::vector<HydroFittedVelocitySampler::VertexTrace> &);
    ~AdaptiveVelocity();
    Vector2D sample(Vector2D) const;
    HydroFittedVelocitySampler::FirstJet firstJet(Vector2D) const;
    HydroFittedVelocitySampler::FirstJet volumeFirstJet(Vector2D) const;
    double divergence(Vector2D) const;
    const HydroFittedFaceFlux &face(int nativeRow) const;

  private:
    struct Impl;
    std::unique_ptr<const Impl> impl_;
};
} // namespace HydroFittedDetail
