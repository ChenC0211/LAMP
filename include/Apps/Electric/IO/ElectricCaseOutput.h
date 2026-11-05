#pragma once

#include "Electric/IO/ElectricCaseInput.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include "Core/IO/TecplotWriter.h"
#include <array>
#include <memory>

// Passive App-stage wall clocks. Zero means this stage was not invoked in this row.
// totalSeconds excludes performance serialization and console formatting.
struct ElectricPerformanceReport {
    double totalSeconds=0, inputSeconds=0, outputStartupSeconds=0;
    double meshSeconds=0, publicationSeconds=0, geometrySeconds=0;
    double panelChargeSeconds=0, momentPrepareSeconds=0, initialSolveSeconds=0;
    double chargePrepareSeconds=0, initialCurrentSeconds=0, chargeAdvanceSeconds=0;
    double maxwellSeconds=0, commitSeconds=0, outputSeconds=0;
};

// Serialization of one accepted electric Q/V/I/Maxwell layer. No solve, charge
// update, geometry publication or physical-state acceptance is performed here.
class ElectricCaseOutput {
public:
    ElectricCaseOutput(const ElectricCaseInput&, int finalStep);
    void write(const QuadtreeMesh&, const LAMP::Electric::ElectricPhaseResult&,
        const LAMP::Electric::ElectricPanelCurrentOperator&,
        const LAMP::Electric::ElectricMaxwellTrace&, int step,
        const LAMP::Electric::ElectricChargeTimeReport* report=nullptr);
    void writePerformance(int step, double time, const ElectricPerformanceReport&,
        const LAMP::Electric::ElectricPhaseResult&,
        const LAMP::Electric::ElectricChargeTimeReport* report=nullptr);
private:
    struct OutputTiming {
        double validation=0, cloudBuild=0, cloudWrite=0, samples=0, series=0;
    } outputTiming_;
    RunControls run_;
    const int finalStep_;
    bool radial_=false;
    std::array<int,2> phaseIds_{{0,1}};
    std::array<double,2> epsilon_{{1,1}}, conductivity_{{0,0}};
    double sourceTolerance_=0, projectionTolerance_=0;
    int lastWrittenStep_=-1; // Output ordering only, not a physical clock.
    std::unique_ptr<TecplotWriter::History> history_;
};
