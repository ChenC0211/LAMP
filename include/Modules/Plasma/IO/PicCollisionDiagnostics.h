#pragma once
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include "Modules/Plasma/Physics/MCCCollisionOperator.h"

class PicSpeciesState;
class UnitsContext;
class QuadtreeMesh;
struct GasState;
struct NormalizationScales;
struct SpeciesConfig;

// Existing explicit-PIC startup diagnostics; does not modify operator, mesh,
// particles, random streams or the driver's integration schedule.
void writeExplicitPicCollisionInitialization(const MCCCollisionOperator& operation,
    const GasState& gas, const UnitsContext* units, const NormalizationScales& scales,
    const QuadtreeMesh& mesh, const std::vector<SpeciesConfig>& species,
    double dt, const std::string& reportPath);

// Own the explicit MCC/ion/neutral report files. The driver chooses sampling
// phases and converts time to SI; no live particle/operator state is retained.
// Ion/neutral files keep the legacy truncate-on-open behavior on restart.
class ExplicitPicCollisionReports
{
public:
    ExplicitPicCollisionReports(const std::filesystem::path& metaDir, int statsEvery);
    const std::string& reportPath() const { return reportPath_; }
    void writeStep(int step, double timeSI, const PicSpeciesState& particles,
        double velocityScale, const MCCCollisionOperator* operation,
        const MCCCollisionOperator::StepStats& stats);
    // Called at the driver's output/final-dump phase when MCC is enabled.
    void append(int step, double timeSI, const PicSpeciesState& particles,
        double velocityScale, const MCCCollisionOperator::StepStats& stats,
        const MCCCollisionOperator::StepStats& cumulative);
private:
    std::string reportPath_;
    int statsEvery_;
    std::ofstream ionOut_, neutralOut_;
};
