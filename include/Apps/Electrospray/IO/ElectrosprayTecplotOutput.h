#pragma once

struct PhaseMaterial;
enum class ElectrosprayProblem;

#include <fstream>
#include <string>
#include <vector>

class QuadtreeMesh;
class HydroFittedDomain;
class Hydrodynamics;
class UnitsContext;
struct ElectrosprayElectricState;
struct ElectrosprayInterfaceState;
struct ElectrosprayDielectricForces;
struct FittedSharpDielectricResult;
struct ElectrosprayConductingForces;
struct SurfaceChargeSource;
namespace Quantum::Electricity { class SurfaceChargeField; class FittedChargeState; }
template<typename T> class VEdgeField;
template<typename T> class HEdgeField;

// App-level SI/output binding. Numerical fields and the interface come from
// Modules/Core; the format writer contains no spray physics. Call only for
// the initial state or a completed coupled step, never a rejected candidate.
class ElectrosprayTecplotOutput
{
public:
    explicit ElectrosprayTecplotOutput(std::string name, bool diagnosticFields = true);
    void write(const QuadtreeMesh&, const HydroFittedDomain&, const Hydrodynamics&,
        const ElectrosprayElectricState&, const PhaseMaterial& liquid, double surfaceTensionSI,
        const UnitsContext&, const ElectrosprayInterfaceState*, const ElectrosprayDielectricForces*,
        const VEdgeField<double>* totalForceX, const HEdgeField<double>* totalForceY,
        int step, double timeCode, ElectrosprayProblem problem);
    // Only an accepted, same-time conducting state. The fixed binding supplies
    // actual graph/boundary IDs, never its old potential/material arrays.
    // Forces are rebuilt at timeCode for the NEXT NS step. Cell electric data
    // remain potential/integrated charge/current; two-sided E/D are explicitly
    // located at reconstructed load-interface traces, not invented cell E.
    void writeConducting(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        const Hydrodynamics& hydro, const std::vector<int>& phaseByCell,
        const ElectrosprayElectricState& binding,
        const Quantum::Electricity::FittedChargeState& electrical,
        const PhaseMaterial& liquid, const PhaseMaterial& gas, double surfaceTensionSI, const UnitsContext& units,
        const ElectrosprayInterfaceState& interfaceState,
        const ElectrosprayConductingForces& forces, int step, double timeCode);
private:
    std::string name_;
    bool diagnosticFields_ = true;
    std::ofstream series_;
    std::ofstream conductingSeries_;
};

// Separate t=0 conducting/P3 output; never relabel the old q=0 flow writer.
// Uses the binding's actual active graph and boundary IDs, the successful
// result's corrected D/phase-current integrals, and the supplied cell-integrated
// surface charge. No nodal potential interpolation, charge advancement, Maxwell
// load or emission model is computed. Nonzero initial charge is supported.
// SI/RZ only. Existing output files are not overwritten; a failed write can
// leave a partial new file and is reported by exception, not as completion.
void writeElectrosprayConductingInitialOutput(const QuadtreeMesh& mesh,
    const ElectrosprayElectricState& initialElectricBinding,
    const FittedSharpDielectricResult& result, const SurfaceChargeSource& charge,
    const UnitsContext& units, const std::string& outputName);

// Initial full Maxwell/capillary native-face data and reconstructed interface
// traces in SI. No flow or charge time advancement is inferred by this writer.
// Uses the same output/fields convention; refuses an existing target.
void writeElectrosprayConductingInitialLoads(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayConductingForces& forces,
    const PhaseMaterial& liquid, const PhaseMaterial& gas, const UnitsContext& units,
    const std::string& outputName);
