#pragma once

#include "Core/Field/CellField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/VertField.h"
#include "Modules/Hydrodynamics/HydroFittedNumerics.h"
#include "Modules/Hydrodynamics/HydroStep.h"

#include <memory>
#include <vector>

class QuadtreeMesh;
class Hydrodynamics;
class HydroFittedDomain;

// One accepted fitted CN/AB2 state, staged independently of its live solver.
// The source mesh must remain alive and unchanged through prepare(). Register
// these fields in EVERY level of the original QuadtreeMesh::adapt operation;
// a single updateWithMesh on the final grid is not a valid transfer.
//
// Scope: fixed physical RZ domain, unchanged wall/axis/reservoir boundary
// faces and boundary data, same-level boundary stencil bands. One-way pressure
// outlets are not yet supported here. This does not replace the fitted step's
// mesh/history locks and does not schedule refinement or advance time.
class HydroFittedMeshTransfer
{
public:
    struct VelocityReport
    {
        double maxDivergenceBefore = 0, maxDivergenceAfter = 0;
        double localContinuityResidual = 0, maxCorrection = 0;
        double sourceMomentumR = 0, sourceMomentumZ = 0;
        double targetMomentumR = 0, targetMomentumZ = 0;
        double sourceBoundaryFlux = 0, targetBoundaryFlux = 0;
        int projectionSolves = 0;
    };
    struct Report
    {
        VelocityReport current, previous;
        double sourcePressureIntegral = 0, targetPressureIntegral = 0;
        bool usedSuperLU = false;
    };
    struct Result
    {
        std::unique_ptr<Hydrodynamics> solver;
        Report report;
        Result();
        ~Result();
        Result(Result&&) noexcept;
        Result& operator=(Result&&) noexcept;
    };

    HydroFittedMeshTransfer(const QuadtreeMesh& sourceMesh, const Hydrodynamics& source);
    ~HydroFittedMeshTransfer();
    HydroFittedMeshTransfer(const HydroFittedMeshTransfer&) = delete;
    HydroFittedMeshTransfer& operator=(const HydroFittedMeshTransfer&) = delete;

    // Owned staging fields, never aliases of the accepted solver. Material
    // closures may be rebuilt here from the transferred interface before
    // prepare(); current/previous velocities and physical pressure retain
    // their original time meaning. Predictor fields are scratch, NOT U(n-2).
    CellField<double> pressure, density, viscosity;
    VertField<double> vertexDensity, vertexViscosity;
    VEdgeField<double> velocityR, previousVelocityR, predictorR;
    HEdgeField<double> velocityZ, previousVelocityZ, predictorZ;

    void appendTransferFields(std::vector<VertField<double>*>& vertices,
        std::vector<HEdgeField<double>*>& horizontal,
        std::vector<VEdgeField<double>*>& vertical,
        std::vector<CellField<double>*>& cells);

    // Uses the ORIGINAL native face density, C/D pressure projection and
    // checked SuperLU solve. Projects U(n) and U(n-1) separately with the same
    // rebuilt material closure, retaining the fixed physical dt. The algebraic
    // multiplier is discarded: it is NOT a physical pressure increment.
    // Original field interpolation is not a conservative momentum/pressure
    // remap; report its changes, do not claim conservation. Failure publishes
    // no solver and leaves source and staging fields unchanged. Target mesh
    // stencil caches may be refreshed; target must be a private candidate.
    Result prepare(QuadtreeMesh& targetMesh, const HydroFittedDomain& targetDomain,
        const HydroFittedNumerics& numerics) const;

private:
    static const Hydrodynamics& checkedSource(const Hydrodynamics& source);
    static std::vector<double> nativeFaceDensities(const QuadtreeMesh& mesh,
        const HydroFittedDomain& domain, const Hydrodynamics& solver);
    const QuadtreeMesh* sourceMesh_;
    std::shared_ptr<const HydroFittedDomain> sourceDomain_;
    double fixedDt_;
    HydroStepProgress stepProgress_;
    HydroMaterialState materialState_;
    Report sourceReport_;
};
