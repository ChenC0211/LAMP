#pragma once

#include "Modules/Hydrodynamics/HydroFittedVelocity.h"

namespace HydroFittedDetail
{
// Prepare one shared physical velocity and axial derivative per actual vertex
// from the accepted native NS face field. No scalar geometry, time advancement,
// solid ghost values or application-specific boundary model is introduced.
// Inactive vertices retain NaN; the adaptive sampler never consumes them.
std::vector<HydroFittedVelocitySampler::VertexTrace> reconstructVertexVelocity(
    const QuadtreeMesh &mesh, const HydroFittedDomain &domain, const VEdgeField<double> &ux,
    const HEdgeField<double> &uy);
} // namespace HydroFittedDetail
