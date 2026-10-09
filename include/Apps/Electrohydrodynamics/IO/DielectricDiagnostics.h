#pragma once

#include <fstream>
#include <string>
namespace Quantum::Hydrodynamics { struct HydroUniformLoad; }
struct DielectricReport;
struct DielectricFields;
class QuadtreeMesh;
template<typename T> class VertField;

// Code-unit EHD diagnostics only; time is the force source time, not the
// just-completed flow time. No solve, reconstruction or state advancement.
class DielectricDiagnostics
{
public:
    void initialize();
    void writeStep(const DielectricReport& report, const Quantum::Hydrodynamics::HydroUniformLoad& load, int step, double forceTime, double dt);
private:
    std::ofstream summary_;
};
void writeDielectricLoadFields(const QuadtreeMesh& mesh,
    int step, double forceTime, const DielectricFields& electric, const Quantum::Hydrodynamics::HydroUniformLoad& load,
    bool axisymmetric, const std::string& outputName);
