#pragma once

#include <string>
#include "Modules/Hydrodynamics/HydroStep.h"
class QuadtreeMesh;
class HydroInterfaceLoadGeometry;
template<typename Mesh> class ImmersedGeometry;
// Complete-domain capillarity/traction adapter over the existing balanced-force
// and legacy reference quadrature kernels. No case, clock or electric callback.
namespace HydroInterfaceForce
{
struct Options
{
    double surfaceTension = 0.0;
    double normalTraction = 0.0;
    double tangentialTraction = 0.0;
    bool boundaryContact = false;
    bool axisContactSupported = false;
};

struct Report
{
    bool legacyReferenceAvailable = true;
    int nodeCount = 0;
    double surfaceMeasure = 0.0;
    double maxAbsCurvature = 0.0;
    double forceTime = 0.0;
    double forceX = 0.0, forceY = 0.0;
    double spreadForceX = 0.0, spreadForceY = 0.0;
    double maxForceDensity = 0.0;
    double activeCapillaryFaceX = 0.0, activeCapillaryFaceY = 0.0;
    double maxActiveCapillaryFaceDensity = 0.0;
    double maxForcePressureImbalance = 0.0, l2ForcePressureImbalance = 0.0;
};

// One complete-domain load assembly. Owns the combined fields; borrows the
// unchanged mesh for this force/solve phase. No provider, IO, clock or solve.
// Add one coupling cell result, then its optional face result. Cell and face
// densities are additive sources, not two representations to apply twice.
// An arithmetic/storage failure invalidates the assembly; no rollback/retry.
class Assembly
{
public:
    Assembly(QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
        int interfaceIndex, bool rz, bool axisContact, const Options& options,
        Report& diagnostic, std::string& stage,
        const HydroInterfaceLoadGeometry* loadGeometry = nullptr);
    Assembly(const Assembly&) = delete;
    Assembly& operator=(const Assembly&) = delete;
    void addCouplingCells(const CellField<double>& x, const CellField<double>& y);
    void addCouplingFaces(const VEdgeField<double>& x, const HEdgeField<double>& y, bool supplied);
    // Freeze once; returned fields borrow this assembly until destruction.
    // The caller attaches the existing Hydro force-source identity before
    // prepareStep; this object does not manufacture source time/provenance.
    HydroForceDensity publish(double evaluationTime);
    // Read-only diagnostic fields after publication, including absent loads.
    const VEdgeField<double>& faceX() const;
    const HEdgeField<double>& faceY() const;
private:
    void requireOpen() const;
    QuadtreeMesh& mesh_;
    CellField<double> cellX_, cellY_;
    VEdgeField<double> faceX_;
    HEdgeField<double> faceY_;
    bool interfaceEnabled_, hasFaces_;
    bool cellsAdded_ = false, facesAdded_ = false, failed_ = false, published_ = false;
};
}
