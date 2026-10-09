#pragma once

#include "Modules/Hydrodynamics/HydroInterfaceForce.h"
#include <vector>

class Hydrodynamics;
class NativeLeafDomain;
namespace Quantum::Hydrodynamics { class HydroConservativePressure; }

// Complete-domain native-field diagnostics, with their original accumulation
// order and acceptance gates. No files, main time loop or App configuration.
namespace HydroFlowDiagnostics
{
struct GeometryMaterialDiagnostics
{
    double interfaceMaxChange = 0.0;
    double negativeVolume = 0.0;
    double rhoMin = 0.0, rhoMax = 0.0;
    double muMin = 0.0, muMax = 0.0;
    int cutCells = 0;
    double levelSetMin = 0.0, levelSetMax = 0.0;
};

struct BoundaryFlowDiagnostics
{
    double inflow = 0.0, outflow = 0.0, net = 0.0, massNet = 0.0;
};
struct NativeConservationDiagnostics {double mass=0,kineticEnergy=0;};
// Uses the same physical half-cell mass measure as native momentum; no
// interpolated nodal velocity or manufactured time interval is involved.
NativeConservationDiagnostics nativeConservation(const QuadtreeMesh&,const Hydrodynamics&,const NativeLeafDomain&);

double axisVelocity(const QuadtreeMesh& mesh, const Hydrodynamics& hydro, bool rz);
// t^n load minus G(p^(n+1)), NOT a full momentum residual.
void collectForcePressureImbalance(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    const VEdgeField<double>& faceForceX, const HEdgeField<double>& faceForceY,
    bool rz, HydroInterfaceForce::Report& diagnostic);
// Same physical pressure operator and half-dual body-force measures as mixed
// conservative NS. Fixed boundary-normal rows are excluded from the norm.
void collectForcePressureImbalance(const QuadtreeMesh&,const Hydrodynamics&,
    const VEdgeField<double>&,const HEdgeField<double>&,const NativeLeafDomain&,
    const Quantum::Hydrodynamics::HydroConservativePressure&,HydroInterfaceForce::Report&);
// negativeVolume retains the original smoothed phase-indicator quadrature;
// it is not conservative VOF or the zero-contour geometric-volume diagnostic.
GeometryMaterialDiagnostics collectDiagnostics(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
    const std::vector<double>& initialLevelSet, bool rz, double interfaceThickness,
    int boundaryCutCells = -1);
BoundaryFlowDiagnostics boundaryFlow(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    bool axisymmetric, double tolerance, bool constantMaterial);
}
