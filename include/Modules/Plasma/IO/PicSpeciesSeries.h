#pragma once

#include "Modules/Plasma/Diagnostics/PicStepDiagnostics.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// Owns species column layout and the legacy species/boundary CSV streams.
// Callers supply diagnostic values, cadence and the original resume step.
class PicSpeciesSeries
{
public:
    enum class Mode { Explicit, Implicit };
    struct Species { std::string name; double charge; };
    struct Moments
    {
        std::size_t available;
        const std::vector<long long>& macroCounts;
        const std::vector<double>& physicalCounts;
        const std::vector<double>& meanDensities;
        double electronMeanEnergy;
        double electronTemperature;
    };
    using Count = PicSpeciesCount;
    using BudgetRow = PicSpeciesBudget;
    template<class MacroCount> using FluxRow = PicBoundaryFlux<MacroCount>;
    using FloatingFlux = FluxRow<double>;
    using IntegerFlux = FluxRow<std::uint64_t>;
    using BoundaryEnergy = PicBoundaryEnergy;

    PicSpeciesSeries(Mode mode, const std::vector<Species>& species);
    PicSpeciesSeries(const PicSpeciesSeries&) = delete;
    PicSpeciesSeries& operator=(const PicSpeciesSeries&) = delete;

    void openPrimary(const std::filesystem::path& path, int resumeStep);
    void openDetail(const std::filesystem::path& path, int resumeStep);
    void openBudget(const std::filesystem::path& path, int resumeStep);
    bool primaryReady() const { return primary_.is_open() && static_cast<bool>(primary_); }
    bool detailReady() const { return detail_.is_open() && static_cast<bool>(detail_); }
    bool budgetReady() const { return budget_.is_open() && static_cast<bool>(budget_); }

    // Borrowed vectors must contain at least 'available' entries. Missing species
    // use the original zero fallback; electron summaries use the selected column.
    void writePrimary(int step, double time_s, const Moments& moments);
    void writeDetail(int step, double time_s, const Moments& moments);
    void flushDetail();

    // Stream one complete budget line in species order without another species
    // array. begin/add/end retain the driver's original calculation/write order.
    // No automatic row completion or retry is performed on destruction.
    void beginBudget(int step, double time_s);
    void addBudget(const BudgetRow& row);
    void endBudget();

    void openBoundary(const std::filesystem::path& path, int resumeStep);
    void openBoundaryEnergy(const std::filesystem::path& path, int resumeStep);
    void openFlux(const std::filesystem::path& path, int resumeStep, bool writeSI);
    bool boundaryReady() const { return boundary_.is_open() && static_cast<bool>(boundary_); }
    bool boundaryEnergyReady() const { return boundaryEnergy_.is_open() && static_cast<bool>(boundaryEnergy_); }
    bool fluxReady() const { return flux_.is_open() && static_cast<bool>(flux_); }
    // Wide rows: species order, then L/R/B/T; energy adds embedded after T.
    // Keep floating (explicit) and integer (implicit) count serialization distinct.
    void beginBoundary(int step, double time_s);
    void addBoundary(const FloatingFlux& row);
    void addBoundary(const IntegerFlux& row);
    void endBoundary();
    void beginBoundaryEnergy(int step, double time_s);
    void addBoundaryEnergy(const BoundaryEnergy& row);
    void endBoundaryEnergy();
    // Long-form values, time and raw species/side names are supplied as before.
    // Explicit uses defaultfloat; implicit uses fixed. Flush cadence stays outside.
    void writeFlux(int step, double time, const char* side, const std::string& species, const FloatingFlux& row);
    void writeFlux(int step, double time, const char* side, const std::string& species, const IntegerFlux& row);
    void flushFlux();

    // Explicit per-species weights/occupancy, in the original defaultfloat format.
    void openWeightStatistics(const std::filesystem::path& path, int resumeStep);
    bool weightStatisticsReady() const { return weights_.is_open() && weights_.good(); }
    void writeWeightStatistics(int step, double timeS, const std::string& species,
        const PicParticleWeightStatistics& row);
    void flushWeightStatistics() { if (weightStatisticsReady()) weights_.flush(); }

private:
    enum class Layout { Species, Budget, Boundary, BoundaryEnergy };
    void open(std::ofstream& stream, const std::filesystem::path& path, int resumeStep,
        Layout layout, const char* role);
    void writeMoments(std::ofstream& stream, int step, double time_s, const Moments& moments);
    const Mode mode_;
    int electronIndex_ = -1;
    std::vector<std::string> tokens_;
    std::ofstream primary_, detail_, budget_;
    std::ofstream boundary_, boundaryEnergy_, flux_;
    std::ofstream weights_;
};
