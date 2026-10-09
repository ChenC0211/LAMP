#pragma once

#include "Modules/Hydrodynamics/HydroFlowDiagnostics.h"
#include <fstream>
#include <string>

struct HydrodynamicsConfig;
struct HydroSolveReport;

// Serialization of the existing complete-domain Hydro outputs. The App owns
// scheduling, supplies exact report/field times and decides acceptance.
// IO neither derives a second clock nor retains runtime configuration.
namespace HydroOutput
{
class DiagnosticsWriter
{
public:
    // The configuration must outlive this writer. No state/mesh is retained.
    DiagnosticsWriter(const HydrodynamicsConfig& config, bool scalarBoundary);
    void appendMidpoint(int step,double time,double dt,double pressureTime,const HydroSolveReport&,
        const HydroFlowDiagnostics::NativeConservationDiagnostics&,int iterations=0,double velocityDefect=0);
    void appendReport(int step, double time, double dt, const HydroSolveReport& report, double axisRadialSpeed,
        const HydroFlowDiagnostics::GeometryMaterialDiagnostics& geometryDiagnostic,
        const HydroInterfaceForce::Report& forceDiagnostic,
        const HydroFlowDiagnostics::BoundaryFlowDiagnostics& flowDiagnostic,
        double geometricVolume, double initialGeometricVolume, double integratedNegativeInflow);

private:
    const HydrodynamicsConfig& config_;
    const bool scalarBoundary_;
    std::ofstream diagnostics;
    std::ofstream energy;
    double initialMass_=0,initialEnergy_=0,dissipation_=0,externalWork_=0;
    int energyStep_=-1;
};

void writeFields(const QuadtreeMesh& mesh, const Hydrodynamics& hydro,
    const ImmersedGeometry<QuadtreeMesh>* geometry, int interfaceIndex,
    const HydrodynamicsConfig& config, int step, double time,
    const std::string& outputName, bool radialAxisSymmetry, bool scalarBoundary,double pressureTime=-1);
}
