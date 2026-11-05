#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "Core/IO/TecplotWriter.h"
#include <ios>
#include <string>
namespace LAMP::Hydrodynamics {class HydroVofStep;class HydroBalancedCsfLoad;}
namespace LAMP::Interface {class GeometricVofMeshTransfer;}
// Pure serialization of borrowed actual fluid receipts; no files, accepted owner or transaction.
namespace AppIO::HydroOutput {
struct Observation {
    std::string capillaryDiscretization="not_observed",scalarRepresentation="not_observed",gradientUpdate="not_observed";
    bool heightFunctionForce=false,incrementForce=false,faceCurvatureDefined=false;
    double faceCurvatureMin=0,faceCurvatureMax=0;
};
struct Totals {long double initialVolume=0,bodyWork=0,viscousLoss=0;};
struct Rows {
    std::string progress,conservation,solver,energyPrefix,energySuffix,capillaryFile,capillary;
    Totals next;long double workCorrected=0;
};
Observation observeVof(const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroVofStep*,const LAMP::Interface::GeometricVofMeshTransfer* = nullptr);
std::string vofBuildInfo(const HydrodynamicsInput&,const Observation&);
std::string vofManifestFields();
std::string vofNativeVelocityCsvHeader();
std::string vofNativeVelocityRow(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
    const LAMP::Hydrodynamics::HydroVofStep*,const LAMP::Interface::GeometricVofMeshTransfer* = nullptr);
Rows vofRows(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroVofStep*,
    const HydrodynamicsInput&,const Totals&,double,std::ios_base::fmtflags,std::streamsize,bool,const LAMP::Interface::GeometricVofMeshTransfer* = nullptr);
TecplotWriter::Zone vofVolumeZone(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroVofStep*,const LAMP::Interface::GeometricVofMeshTransfer* = nullptr);
// Diagnostic of the retained force owner on its original mesh, including
// through an accepted End remap. Never substitutes the accepted frame geometry.
TecplotWriter::Zone actualForceGeometryZone(const QuadtreeMesh& sourceMesh,const LAMP::Hydrodynamics::HydroBalancedCsfLoad&);
TecplotWriter::Zone vofInterfaceZone(const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroVofStep*,bool,const LAMP::Interface::GeometricVofMeshTransfer* = nullptr);
std::string vofPerformance(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroVofStep*,double,double,double,const LAMP::Interface::GeometricVofMeshTransfer* = nullptr);
Observation observe(const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroTwoPhaseStep*,const Observation&);
Rows rows(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroTwoPhaseStep*,
    const HydrodynamicsInput&,const Totals&,const Observation&,double compute,std::ios_base::fmtflags,std::streamsize,bool final);
std::string performance(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,const LAMP::Hydrodynamics::HydroTwoPhaseStep*,
    const HydrodynamicsInput&,const Observation&,const char* acceptedStage,double compute,double io,double elapsed);
std::string buildInfo(const HydrodynamicsInput&,const Observation&);
std::string manifestFields(const HydrodynamicsInput&,const Observation&);
std::string volumeXml(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&);
TecplotWriter::Zone interfaceZone(const LAMP::Hydrodynamics::HydroTwoPhaseState&,bool radial);
std::uint64_t geometryIdentityValue(const LAMP::Hydrodynamics::HydroTwoPhaseState&);
}

