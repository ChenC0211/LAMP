#pragma once

#include "Modules/Hydrodynamics/HydroStep.h"

struct HydroFittedCapillaryState;

// Completed capillary + external native force density on one fitted domain.
// Shares immutable owned fields and their producer's provenance; copying this
// handle does not copy arrays. No electric model, unit system, clock or IO.
class HydroFittedForce
{
public:
    HydroFittedForce() = default; // No published load yet.
    static HydroFittedForce combine(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const HydroFittedCapillaryState& capillary,
        const VEdgeField<double>& externalX, const HEdgeField<double>& externalY,
        double externalScale = 1.0);

    std::shared_ptr<const HydroForceSource> source() const noexcept;
    HydroReadOnlyField<VEdgeField<double>> x() const;
    HydroReadOnlyField<HEdgeField<double>> y() const;
    // Synchronous numerical borrow, valid while a handle retains this load.
    // Fields must not be modified. Evaluation time initially labels the source;
    // an App may evaluate a Frozen source later, subject to Hydro's step checks.
    HydroForceDensity density() const;
private:
    struct Data;
    const Data& data() const;
    std::shared_ptr<const Data> data_;
};
