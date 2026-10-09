#pragma once

#include <fstream>
#include <string>

namespace Quantum::Hydrodynamics { struct HydroUniformLoad; }
class QuadtreeMesh;
class HydroInterfaceLoadGeometry;
struct DielectricFields;
struct DielectricConfig;
struct DielectricReport;
struct DielectricSharpAssemblyDiagnostics;
template<typename T> class VertField;

// EHD App I/O: combine committed electric fields and Hydro loads from the same
// solve and its borrowed t^n geometry; no force/field is solved or fitted here.
// Separate files preserve the original electric CSV/DAT schemas.
class SharpDielectricDiagnostics
{
public:
    void initialize();
    void writeStep(const QuadtreeMesh& mesh, const VertField<double>& phi,
        const HydroInterfaceLoadGeometry& geometry, const DielectricFields& electric, const Quantum::Hydrodynamics::HydroUniformLoad& load,
        const DielectricConfig& config, const DielectricReport& report,
        int step, double forceTime, double dt, double materialThickness,
        const std::string& outputName, bool writeFields,
        const DielectricSharpAssemblyDiagnostics* assemblyDiagnostics);

private:
    std::ofstream summary_;
};
