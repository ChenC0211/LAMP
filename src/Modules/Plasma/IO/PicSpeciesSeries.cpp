#include "Modules/Plasma/IO/PicSpeciesSeries.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include <cctype>
#include <iomanip>
#include <iostream>
#include <unordered_set>

namespace
{
constexpr const char* primaryFields[] = {
    "macro_count", "phys_count", "n_mean_1pm3", "mean_energy_eV", "Te_eV"
};
constexpr const char* budgetFields[] = {
    "macro_begin", "phys_begin", "macro_source", "phys_source",
    "macro_mcc_spawn", "phys_mcc_spawn", "macro_see_source", "phys_see_source",
    "macro_out_domain", "phys_out_domain", "macro_out_embedded", "phys_out_embedded",
    "macro_end", "phys_end", "macro_particle_management", "phys_particle_management",
    "macro_reaction_net", "phys_reaction_net", "macro_residual", "phys_residual"
};
constexpr const char* boundaryDomainTokens[] = {"l", "r", "b", "t"};
constexpr const char* boundaryAllTokens[] = {"l", "r", "b", "t", "embedded"};
constexpr const char* boundaryFields[] = {
    "macro_out", "phys_out", "flux_macro_1pm2s", "flux_phys_1pm2s", "current_Apm2"
};
constexpr const char* boundaryEnergyFields[] = {
    "macro_out", "phys_out", "absorbed_energy_J", "absorbed_power_W",
    "mean_energy_eV", "energy_flux_Wpm2", "surface_area_m2"
};

template<class MacroCount>
void writeFluxValues(std::ostream& out, const PicSpeciesSeries::FluxRow<MacroCount>& row)
{
    out << ',' << row.macro << ',' << row.physical
        << ',' << row.macroFlux << ',' << row.physicalFlux << ',' << row.current;
}

template<class MacroCount>
void writeFluxRecord(std::ostream& out, int step, double time, const char* side,
    const std::string& species, const PicSpeciesSeries::FluxRow<MacroCount>& row)
{
    out << step << ',' << std::setprecision(16) << time << ',' << side << ',' << species;
    writeFluxValues(out, row);
    out << "\n";
}
}

PicSpeciesSeries::PicSpeciesSeries(Mode mode, const std::vector<Species>& species) : mode_(mode)
{
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        std::string name = species[si].name;
        for (auto& ch : name)
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        if (name == "e" || name == "electron" || name.find("electron") != std::string::npos)
        {
            electronIndex_ = static_cast<int>(si);
            break;
        }
    }
    if (electronIndex_ < 0)
    {
        for (std::size_t si = 0; si < species.size(); ++si)
        {
            if (species[si].charge < 0.0)
            {
                electronIndex_ = static_cast<int>(si);
                break;
            }
        }
    }
    tokens_.resize(species.size());
    std::unordered_set<std::string> used;
    for (std::size_t si = 0; si < species.size(); ++si)
    {
        const std::string base = TimeSeriesIO::sanitizeCsvToken(species[si].name);
        std::string token = base;
        int suffix = 2;
        while (used.count(token) > 0u)
            token = base + "_" + std::to_string(suffix++);
        used.insert(token);
        tokens_[si] = token;
    }
}

void PicSpeciesSeries::openWeightStatistics(const std::filesystem::path& path, int resumeStep)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(weights_, path, resumeStep);
    if (weights_)
    {
        if (!hasExisting)
            weights_ << "step,time_s,species,n_macro,n_phys,w_min,w_mean,w_p50,w_p90,w_p99,w_max,"
                << "cell_macro_max,cell_phys_max,cell_index_macro_max,cell_index_phys_max\n";
        weights_ << std::setprecision(16);
    }
    else
        std::cerr << "[ExplicitPIC] WARNING: failed to create particle weight stats: " << path.string() << std::endl;
}
void PicSpeciesSeries::writeWeightStatistics(int step, double timeS, const std::string& species,
    const PicParticleWeightStatistics& row)
{
    weights_ << step << ',' << timeS << ',' << species << ',' << row.macro << ',' << row.physical << ','
        << row.minimum << ',' << row.mean << ',' << row.percentiles[0] << ',' << row.percentiles[1] << ','
        << row.percentiles[2] << ',' << row.maximum << ',' << row.cellMacroMaximum << ',' << row.cellPhysicalMaximum << ','
        << row.cellMacroIndex << ',' << row.cellPhysicalIndex << "\n";
}

void PicSpeciesSeries::open(std::ofstream& stream, const std::filesystem::path& path,
    int resumeStep, Layout layout, const char* role)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(stream, path, resumeStep);
    if (stream)
    {
        if (!hasExisting)
        {
            stream << "step,time_s";
            for (const auto& token : tokens_)
            {
                if (layout == Layout::Budget)
                    for (const char* field : budgetFields) stream << ',' << token << '_' << field;
                else if (layout == Layout::Boundary)
                    for (const char* side : boundaryDomainTokens)
                        for (const char* field : boundaryFields) stream << ',' << token << '_' << side << '_' << field;
                else if (layout == Layout::BoundaryEnergy)
                    for (const char* side : boundaryAllTokens)
                        for (const char* field : boundaryEnergyFields) stream << ',' << token << '_' << side << '_' << field;
                else
                    for (const char* field : primaryFields) stream << ',' << token << '_' << field;
            }
            stream << "\n";
        }
        stream << std::setprecision(16) << std::fixed;
        if (mode_ == Mode::Implicit) stream.flush();
    }
    else
    {
        std::cerr << '[' << (mode_ == Mode::Explicit ? "ExplicitPIC" : "ImplicitPIC")
                  << "] WARNING: failed to create " << role << " csv: " << path.string() << std::endl;
    }
}

void PicSpeciesSeries::openPrimary(const std::filesystem::path& path, int resumeStep)
{ open(primary_, path, resumeStep, Layout::Species, "species primary"); }
void PicSpeciesSeries::openDetail(const std::filesystem::path& path, int resumeStep)
{ open(detail_, path, resumeStep, Layout::Species, "species detail"); }
void PicSpeciesSeries::openBudget(const std::filesystem::path& path, int resumeStep)
{ open(budget_, path, resumeStep, Layout::Budget, "species budget"); }

void PicSpeciesSeries::writeMoments(std::ofstream& stream, int step, double time_s, const Moments& moments)
{
    stream << step << ',' << time_s;
    for (std::size_t si = 0; si < tokens_.size(); ++si)
    {
        const bool have = (si < moments.available);
        const bool isElectron = (electronIndex_ >= 0) && (static_cast<std::size_t>(electronIndex_) == si);
        stream << ',' << (have ? moments.macroCounts[si] : 0ll)
               << ',' << (have ? moments.physicalCounts[si] : 0.0)
               << ',' << (have ? moments.meanDensities[si] : 0.0)
               << ',' << (isElectron ? moments.electronMeanEnergy : 0.0)
               << ',' << (isElectron ? moments.electronTemperature : 0.0);
    }
    stream << "\n";
    stream.flush();
}

void PicSpeciesSeries::writePrimary(int step, double time_s, const Moments& moments)
{ writeMoments(primary_, step, time_s, moments); }
void PicSpeciesSeries::writeDetail(int step, double time_s, const Moments& moments)
{ writeMoments(detail_, step, time_s, moments); }
void PicSpeciesSeries::flushDetail() { if (detailReady()) detail_.flush(); }

void PicSpeciesSeries::beginBudget(int step, double time_s) { budget_ << step << ',' << time_s; }
void PicSpeciesSeries::addBudget(const BudgetRow& row)
{
    budget_ << ',' << row.begin.macro << ',' << row.begin.physical
            << ',' << row.source.macro << ',' << row.source.physical
            << ',' << row.mccSpawn.macro << ',' << row.mccSpawn.physical
            << ',' << row.seeSource.macro << ',' << row.seeSource.physical
            << ',' << row.outDomain.macro << ',' << row.outDomain.physical
            << ',' << row.outEmbedded.macro << ',' << row.outEmbedded.physical
            << ',' << row.end.macro << ',' << row.end.physical
            << ',' << row.particleManagement.macro << ',' << row.particleManagement.physical
            << ',' << row.reactionNet.macro << ',' << row.reactionNet.physical
            << ',' << row.residual.macro << ',' << row.residual.physical;
}
void PicSpeciesSeries::endBudget() { budget_ << "\n"; budget_.flush(); }

void PicSpeciesSeries::openBoundary(const std::filesystem::path& path, int resumeStep)
{ open(boundary_, path, resumeStep, Layout::Boundary, "boundary primary"); }
void PicSpeciesSeries::openBoundaryEnergy(const std::filesystem::path& path, int resumeStep)
{ open(boundaryEnergy_, path, resumeStep, Layout::BoundaryEnergy, "boundary energy"); }

void PicSpeciesSeries::openFlux(const std::filesystem::path& path, int resumeStep, bool writeSI)
{
    const bool hasExisting = TimeSeriesIO::openPreparedStepCsv(flux_, path, resumeStep);
    if (flux_)
    {
        if (!hasExisting)
        {
            if (mode_ == Mode::Implicit || writeSI)
                flux_ << "step,time_s,side,species,macro_out,phys_out,flux_macro_1pm2s,flux_phys_1pm2s,current_Apm2\n";
            else
                flux_ << "step,time,side,species,macro_out,phys_out,flux_macro,flux_phys,current\n";
        }
        flux_ << std::setprecision(16);
        if (mode_ == Mode::Implicit) { flux_ << std::fixed; flux_.flush(); }
    }
    else if (mode_ == Mode::Explicit)
        std::cerr << "[ExplicitPIC] WARNING: failed to create flux report: " << path.string() << std::endl;
    else
        std::cerr << "[ImplicitPIC] WARNING: failed to create boundary flux csv: " << path.string() << std::endl;
}

void PicSpeciesSeries::beginBoundary(int step, double time_s)
{ boundary_ << step << ',' << std::setprecision(16) << time_s; }
void PicSpeciesSeries::addBoundary(const FloatingFlux& row) { writeFluxValues(boundary_, row); }
void PicSpeciesSeries::addBoundary(const IntegerFlux& row) { writeFluxValues(boundary_, row); }
void PicSpeciesSeries::endBoundary() { boundary_ << "\n"; boundary_.flush(); }
void PicSpeciesSeries::beginBoundaryEnergy(int step, double time_s)
{ boundaryEnergy_ << step << ',' << std::setprecision(16) << time_s; }
void PicSpeciesSeries::addBoundaryEnergy(const BoundaryEnergy& row)
{
    boundaryEnergy_ << ',' << row.out.macro << ',' << row.out.physical
        << ',' << row.absorbedEnergyJ << ',' << row.absorbedPowerW << ',' << row.meanEnergyEV
        << ',' << row.energyFluxWpm2 << ',' << row.surfaceAreaM2;
}
void PicSpeciesSeries::endBoundaryEnergy() { boundaryEnergy_ << "\n"; boundaryEnergy_.flush(); }
void PicSpeciesSeries::writeFlux(int step, double time, const char* side,
    const std::string& species, const FloatingFlux& row)
{ writeFluxRecord(flux_, step, time, side, species, row); }
void PicSpeciesSeries::writeFlux(int step, double time, const char* side,
    const std::string& species, const IntegerFlux& row)
{ writeFluxRecord(flux_, step, time, side, species, row); }
void PicSpeciesSeries::flushFlux() { flux_.flush(); }
