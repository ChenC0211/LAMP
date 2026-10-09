#include "Modules/Plasma/IO/PicOutputManifest.h"
#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/Diagnostics/PhysicsDiagnostics.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <sstream>
#include <utility>
#include <vector>

namespace
{
constexpr const char* explicitObsolete[] = {
    "reports/collisions_report.txt", "reports/see_report.txt", "reports/ion_report.txt", "reports/neutral_report.txt",
    "series/impact_report.csv", "series/particle_management_report.csv", "series/flux_report.csv",
    "series/rf_phase_report.csv", "series/particle_weight_stats.csv", "series/periodic_diagnostics.csv",
    "series/current_waveform.csv", "series/current_balance.csv", "series/eb_band_rms.csv",
    "series/timing.csv", "series/pic_diagnostics.csv"
};
constexpr const char* implicitObsolete[] = {
    "series/coupling_fallbacks.csv", "series/current_balance.csv", "series/geometry_metrics.csv",
    "series/eb_band_rms.csv", "series/periodic_diagnostics.csv", "series/current_waveform.csv"
};
struct File
{
    const char* path;
    const char* category;
    const char* role;
    const char* cadenceKey;
    const char* columns;
    int schemaOrder; // Preserve the legacy schema order independently of run order.
};
constexpr const char* explicitOverviewColumns = "step,time_s,rG_L2,rG_Linf,rG_L2_per_active,rG_L2_rel_rhs,rG_Linf_rel_rhs,S_eb_rel_rhs,field_energy_J,kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,total_energy_J,J_cond_Apm2,J_disp_Apm2,J_total_Apm2,macro_total,phys_total,charge_balance_rel";
constexpr const char* processColumns = "step,time_s,drive_energy_J,injection_energy_J,reaction_energy_loss_J,domain_absorption_J,embedded_absorption_J,see_emission_energy_J,ion_impact_energy_J,displacement_current_fraction,dominant_source,dominant_source_fraction,dominant_boundary,dominant_boundary_absorbed_power_W,dominant_boundary_absorbed_fraction,dominant_reaction_channel,dominant_reaction_power_W,dominant_reaction_fraction";
constexpr const char* geometryColumns = "step,time_s,geometry_updated,active_cut_cell_count,embedded_surface_area_m2,moving_region_count,moving_levelset_count,moving_region_total_area_m2,moving_region_total_flux_1pm2s,moving_region_mean_energy_eV,moving_region_mean_cos,moving_region_mean_abs_speed_mps,moving_region_max_abs_speed_mps,geometry_gsa_mm,geometry_gad_mm,geometry_da_mm,field_energy_J,domain_absorption_J,embedded_absorption_J,reaction_energy_loss_J";
constexpr const char* fluxColumns = "step,time_s,side,species,macro_out,phys_out,flux_macro_1pm2s,flux_phys_1pm2s,current_Apm2";
constexpr const char* reactionEnergyColumns = "step,time_s,projectile,channel,macro_collisions,phys_collisions,projectile_energy_loss_J_total,projectile_energy_loss_power_W,projectile_energy_loss_fraction_step";
constexpr const char* implicitOverviewColumns = "step,time_s,rG_L2,rG_Linf,field_energy_J,kinetic_e_J,kinetic_ion_J,kinetic_neu_J,kinetic_total_J,total_energy_J,J_cond_Apm2,J_disp_Apm2,J_total_Apm2,macro_total,phys_total,charge_balance_rel";
constexpr const char* managementColumns = "step,time_s,kind,name,macro_before,phys_before,macro_after,phys_after,merged_pairs,split_count,split_emergency_count,split_count_based_count,spawn_in,spawn_kept,spawn_drop,keep_prob,weight_p50,weight_p90,weight_p99,weight_max,weight_top1_phys_frac";
const File explicitFiles[] = {
    {"fields/profile1d.csv", "fields", "profile", "fieldsEvery", nullptr, 0},
    {"fields/profiles_avg.csv", "fields", "profile", "fieldsEvery", nullptr, 1},
    {"series/overview_step.csv", "series", "primary", "seriesEvery", explicitOverviewColumns, 2},
    {"series/overview_detail_step.csv", "series", "secondary", "seriesEvery", nullptr, 3},
    {"series/species_step.csv", "series", "primary", "seriesEvery", nullptr, 4},
    {"series/species_budget_step.csv", "series", "primary", "seriesEvery", nullptr, 5},
    {"series/boundary_step.csv", "series", "primary", "seriesEvery", nullptr, 6},
    {"series/boundary_energy_step.csv", "series", "primary", "seriesEvery", nullptr, 7},
    {"series/process_step.csv", "series", "primary", "seriesEvery", processColumns, 8},
    {"series/geometry_coupling_step.csv", "series", "primary", "seriesEvery", geometryColumns, 9},
    {"series/perf_step.csv", "series", "primary", "seriesEvery", nullptr, 10},
    {"series/boundary_flux_step.csv", "series", "secondary", "seriesEvery", fluxColumns, 11},
    {"series/power_step.csv", "series", "primary", "fieldsEvery", nullptr, 12},
    {"series/impact_step.csv", "series", "secondary", "seriesEvery", nullptr, 13},
    {"series/particle_management_step.csv", "series", "secondary", "seriesEvery", nullptr, 14},
    {"series/particle_weight_stats_step.csv", "series", "secondary", "seriesEvery", nullptr, 15},
    {"series/rf_phase_step.csv", "series", "secondary", "seriesEvery", nullptr, 16},
    {"series/periodic_step.csv", "series", "secondary", "seriesEvery", nullptr, 17},
    {"series/waveform_step.csv", "series", "secondary", "seriesEvery", nullptr, 18},
    {"series/current_balance_step.csv", "series", "secondary", "seriesEvery", nullptr, 19},
    {"series/eb_band_rms_step.csv", "series", "secondary", "seriesEvery", nullptr, 20},
    {"series/diagnostics/electrodes/current_boundaries_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 21},
    {"series/diagnostics/electrodes/current_groups_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 22},
    {"series/diagnostics/reactions/reaction_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 23},
    {"series/diagnostics/reactions/reaction_channel_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 24},
    {"series/diagnostics/reactions/reaction_energy_step.csv", "series", "diagnostic", "diagnosticsEvery", reactionEnergyColumns, 25},
    {"series/diagnostics/reactions/reaction_projectile_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 26},
    {"series/diagnostics/reactions/reaction_spawn_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 27},
    {"series/diagnostics/energy/histogram_step.csv", "series", "diagnostic", "diagnosticsEvery", nullptr, 28},
    {"meta/run_manifest.json", "meta", "manifest", "startup", nullptr, 29},
    {"meta/schema_manifest.json", "meta", "manifest", "startup", nullptr, 30},
    {"meta/collisions_report.txt", "meta", "report", "startup+append", nullptr, 31},
    {"meta/see_report.txt", "meta", "report", "startup+append", nullptr, 32},
    {"meta/ion_report.txt", "meta", "report", "startup+append", nullptr, 33},
    {"meta/neutral_report.txt", "meta", "report", "startup+append", nullptr, 34},
};
const File implicitFiles[] = {
    {"fields/profile1d.csv", "fields", "profile", "fieldsEvery", nullptr, 0},
    {"fields/profiles_avg.csv", "fields", "profile", "fieldsEvery", nullptr, 1},
    {"series/overview_step.csv", "series", "primary", "physicsEvery", implicitOverviewColumns, 2},
    {"series/perf_step.csv", "series", "primary", "perfEvery", nullptr, 3},
    {"series/physics_step.csv", "series", "secondary", "physicsEvery", nullptr, 4},
    {"series/species_step.csv", "series", "primary", "speciesEvery", nullptr, 5},
    {"series/species_budget_step.csv", "series", "primary", "speciesEvery", nullptr, 8},
    {"series/particle_management_step.csv", "series", "secondary", "speciesEvery", managementColumns, 9},
    {"series/species_detail_step.csv", "series", "secondary", "speciesEvery", nullptr, 10},
    {"series/process_step.csv", "series", "primary", "physicsEvery", processColumns, 6},
    {"series/geometry_coupling_step.csv", "series", "primary", "physicsEvery", geometryColumns, 7},
    {"series/amount_step.csv", "series", "secondary", "amountEvery", nullptr, 11},
    {"series/amount_species_step.csv", "series", "secondary", "amountEvery", nullptr, 12},
    {"series/boundary_step.csv", "series", "primary", "boundaryEvery", nullptr, 13},
    {"series/boundary_energy_step.csv", "series", "primary", "boundaryEvery", nullptr, 14},
    {"series/boundary_flux_step.csv", "series", "secondary", "boundaryEvery", fluxColumns, 15},
    {"series/power_step.csv", "series", "primary", "physicsEvery", nullptr, 16},
    {"series/coupling_fallback_step.csv", "series", "secondary", "perfEvery", nullptr, 17},
    {"series/current_balance_step.csv", "series", "secondary", "boundaryEvery", nullptr, 18},
    {"series/geometry_step.csv", "series", "secondary", "startup+geometryUpdate", nullptr, 19},
    {"series/eb_band_rms_step.csv", "series", "secondary", "physicsEvery", nullptr, 20},
    {"series/periodic_step.csv", "series", "secondary", "physicsEvery", nullptr, 21},
    {"series/waveform_step.csv", "series", "secondary", "physicsEvery", nullptr, 22},
    {"series/amr_adapt_step.csv", "series", "secondary", "geometryUpdate", nullptr, 23},
    {"series/amr_levelset_diag.csv", "series", "secondary", "geometryUpdate", nullptr, 24},
    {"series/diagnostics/electrodes/current_boundaries_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 25},
    {"series/diagnostics/electrodes/current_groups_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 26},
    {"series/diagnostics/reactions/reaction_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 27},
    {"series/diagnostics/reactions/reaction_channel_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 28},
    {"series/diagnostics/reactions/reaction_energy_step.csv", "series", "diagnostic", "physicsEvery", reactionEnergyColumns, 29},
    {"series/diagnostics/reactions/reaction_projectile_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 30},
    {"series/diagnostics/reactions/reaction_spawn_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 31},
    {"series/diagnostics/energy/histogram_step.csv", "series", "diagnostic", "physicsEvery", nullptr, 32},
    {"meta/run_manifest.json", "meta", "manifest", "startup", nullptr, 33},
    {"meta/schema_manifest.json", "meta", "manifest", "startup", nullptr, 34},
};

std::string jsonEscape(const std::string& text)
{
    std::string result;
    result.reserve(text.size() + 16);
    for (const unsigned char ch : text)
    {
        switch (ch)
        {
        case '\\': result += "\\\\"; break;
        case '"': result += "\\\""; break;
        case '\n': result += "\\n"; break;
        case '\r': result += "\\r"; break;
        case '\t': result += "\\t"; break;
        default:
            if (ch < 0x20)
            {
                constexpr char hex[] = "0123456789abcdef";
                result += "\\u00";
                result += hex[ch >> 4];
                result += hex[ch & 15];
            }
            else result.push_back(static_cast<char>(ch));
        }
    }
    return result;
}

void writeManifest(const PicOutputManifest::Input& input, const char* module,
    const File* files, std::size_t fileCount,
    std::initializer_list<std::pair<const char*, int>> cadence)
{
    const auto& cfg = input.config;
    const auto& outputRoot = input.outputRoot;
    const auto fieldsDir = outputRoot / "fields";
    const auto seriesDir = outputRoot / "series";
    const auto metaDir = outputRoot / "meta";
    const auto manifestPath = metaDir / "run_manifest.json";
    std::ofstream manifest(manifestPath, std::ios::out | std::ios::trunc);
    if (manifest)
    {
        manifest << "{\n";
        manifest << "  \"schemaVersion\": 4,\n";
        manifest << "  \"module\": \"" << module << "\",\n";
        manifest << "  \"casePath\": \"" << jsonEscape(std::filesystem::current_path().string()) << "\",\n";
        manifest << "  \"dictPath\": \"" << jsonEscape(input.dictionaryPath) << "\",\n";
        manifest << "  \"output\": {\n";
        manifest << "    \"root\": \"" << jsonEscape(outputRoot.string()) << "\",\n";
        manifest << "    \"fieldsDir\": \"" << jsonEscape(fieldsDir.string()) << "\",\n";
        manifest << "    \"seriesDir\": \"" << jsonEscape(seriesDir.string()) << "\",\n";
        manifest << "    \"metaDir\": \"" << jsonEscape(metaDir.string()) << "\"\n";
        manifest << "  },\n";
        manifest << "  \"cadence\": {\n";
        std::size_t cadenceIndex = 0;
        for (const auto& entry : cadence)
            manifest << "    \"" << entry.first << "\": " << entry.second
                     << (++cadenceIndex < cadence.size() ? ",\n" : "\n");
        manifest << "  },\n";
        if (input.units != nullptr)
        {
            const auto& us = input.units->scales();
            manifest << "  \"units\": {\n";
            manifest << "    \"L0_m\": " << us.L0_m << ",\n";
            manifest << "    \"T0_s\": " << us.T0_s << ",\n";
            manifest << "    \"Phi0_V\": " << us.Phi0_V << ",\n";
            manifest << "    \"n0_m3\": " << us.n0_m3 << ",\n";
            manifest << "    \"v0_mps\": " << us.v0_mps << ",\n";
            manifest << "    \"E0_Vpm\": " << us.E0_Vpm << ",\n";
            manifest << "    \"rho0_Cpm3\": " << us.rho0_Cpm3 << "\n";
            manifest << "  },\n";
        }
        if (input.normalization != nullptr)
        {
            const auto& nc = *input.normalization;
            manifest << "  \"normalization\": {\n";
            manifest << "    \"dimension\": " << nc.dimension << ",\n";
            manifest << "    \"slabThickness_m\": " << nc.slabThickness_m << "\n";
            manifest << "  },\n";
        }
        manifest << "  \"species\": [\n";
        for (std::size_t si = 0; si < cfg.species.size(); ++si)
        {
            const auto& sp = cfg.species[si];
            manifest << "    {\"name\":\"" << jsonEscape(sp.name) << "\",\"q_C\":" << sp.q
                     << ",\"m_kg\":" << sp.m << ",\"weight\":" << sp.weight << "}";
            if (si + 1 < cfg.species.size())
                manifest << ",";
            manifest << "\n";
        }
        manifest << "  ],\n";
        manifest << "  \"fieldTokenToTecplot\": {\n";
        manifest << "    \"phi\": \"phi_SI_V\",\n";
        manifest << "    \"rho\": \"rho_SI_Cpm3\",\n";
        manifest << "    \"E\": [\"Ex_SI_Vpm\",\"Ey_SI_Vpm\"],\n";
        manifest << "    \"n\": \"n_<species>_SI_1pm3\",\n";
        manifest << "    \"Te_e\": \"Te_e_eV\",\n";
        manifest << "    \"electron_mean_energy\": \"electron_mean_energy_eV\",\n";
        manifest << "    \"psi\": \"psi_SI_m\"\n";
        manifest << "  },\n";
        manifest << "  \"selectedFieldTokens\": [";
        for (std::size_t i = 0; i < cfg.output.fields.size(); ++i)
        {
            manifest << "\"" << jsonEscape(cfg.output.fields[i]) << "\"";
            if (i + 1 < cfg.output.fields.size())
                manifest << ",";
        }
        manifest << "],\n";
        manifest << "  \"electrodeCurrentDiagnostics\": {\n";
        manifest << "    \"enabled\": " << (cfg.output.diagnostics.electrodeCurrent.enabled ? "true" : "false") << ",\n";
        manifest << "    \"perBoundary\": " << (cfg.output.diagnostics.electrodeCurrent.perBoundary ? "true" : "false") << ",\n";
        manifest << "    \"includeDomain\": " << (cfg.output.diagnostics.electrodeCurrent.includeDomain ? "true" : "false") << ",\n";
        manifest << "    \"groups\": [\n";
        for (std::size_t gi = 0; gi < cfg.output.diagnostics.electrodeCurrent.groups.size(); ++gi)
        {
            const auto& group = cfg.output.diagnostics.electrodeCurrent.groups[gi];
            manifest << "      {\"name\":\"" << jsonEscape(group.name) << "\",\"boundaryIds\":[";
            for (std::size_t bi = 0; bi < group.boundaryIds.size(); ++bi)
            {
                manifest << group.boundaryIds[bi];
                if (bi + 1 < group.boundaryIds.size())
                    manifest << ",";
            }
            manifest << "]}";
            if (gi + 1 < cfg.output.diagnostics.electrodeCurrent.groups.size())
                manifest << ",";
            manifest << "\n";
        }
        manifest << "    ]\n";
        manifest << "  },\n";
        manifest << "  \"outputFiles\": [\n";
        for (std::size_t n = 0; n < fileCount; ++n)
        {
            const auto& file = files[n];
            if (n != 0) manifest << ",\n";
            manifest << "    {\"path\":\"" << file.path << "\",\"category\":\"" << file.category
                     << "\",\"role\":\"" << file.role << "\",\"cadenceKey\":\"" << file.cadenceKey << "\"}";
        }
        manifest << "\n  ]\n";
        manifest << "}\n";
        manifest.flush();
    }
    else
        std::cerr << '[' << module << "] WARNING: failed to create run manifest: "
                  << manifestPath.string() << std::endl;

    const auto schemaManifestPath = metaDir / "schema_manifest.json";
    std::ofstream schema(schemaManifestPath, std::ios::out | std::ios::trunc);
    if (schema)
    {
        schema << "{\n";
        schema << "  \"schemaVersion\": 2,\n";
        schema << "  \"module\": \"" << module << "\",\n";
        schema << "  \"files\": [\n";
        std::vector<const File*> ordered;
        ordered.reserve(fileCount);
        for (std::size_t n = 0; n < fileCount; ++n)
            ordered.push_back(files + n);
        std::sort(ordered.begin(), ordered.end(), [](const File* a, const File* b) {
            return a->schemaOrder < b->schemaOrder;
        });
        for (std::size_t n = 0; n < ordered.size(); ++n)
        {
            const auto& file = *ordered[n];
            if (n != 0) schema << ",\n";
            schema << "    {\"path\":\"" << file.path << "\",\"category\":\""
                   << file.category << "\",\"role\":\"" << file.role << "\"";
            if (file.columns)
            {
                schema << ",\"columns\":[";
                std::istringstream columns(file.columns);
                std::string column;
                bool first = true;
                while (std::getline(columns, column, ','))
                {
                    if (!first) schema << ',';
                    first = false;
                    schema << "{\"name\":\"" << column << "\"}";
                }
                schema << ']';
            }
            schema << '}';
        }
        schema << "\n  ]\n}\n";
        schema.flush();
    }
    else
        std::cerr << '[' << module << "] WARNING: failed to create schema manifest: "
                  << schemaManifestPath.string() << std::endl;
}
}

PicOutputDirectory::PicOutputDirectory(Mode mode, const std::filesystem::path& diagnosticsBase,
    std::filesystem::path outputRoot)
    : root(std::move(outputRoot)), fields(root / "fields"), series(root / "series"),
      meta(root / "meta"), diagnostics(series / diagnosticsBase), mode_(mode)
{}

void PicOutputDirectory::create() const
{
    if (mode_ == Mode::Explicit) std::filesystem::create_directories(root);
    std::filesystem::create_directories(fields);
    std::filesystem::create_directories(series);
    std::filesystem::create_directories(meta);
}
void PicOutputDirectory::prepareFields(int resumeStep) const
{ TimeSeriesIO::prepareStepTaggedDirectory(fields, resumeStep); }
void PicOutputDirectory::removeObsoleteReports() const
{
    const auto* names = mode_ == Mode::Explicit ? explicitObsolete : implicitObsolete;
    const auto count = mode_ == Mode::Explicit ? std::size(explicitObsolete) : std::size(implicitObsolete);
    for (std::size_t i = 0; i < count; ++i)
    {
        std::error_code ec;
        std::filesystem::remove(root / names[i], ec);
    }
}

PicOutputCadence::PicOutputCadence(Mode mode, const PicConfig& config, int finalStep)
    : mode_(mode), finalStep_(finalStep), fieldsEvery_(std::max(1, config.output.every)),
      seriesEvery_(std::max(1, config.output.v2.seriesEvery)),
      diagnosticsEvery_(PhysicsDiagnostics::resolveEvery(seriesEvery_, config.output.diagnostics.every)),
      amountEvery_(config.output.v2.amountEvery > 0 ? std::max(1, config.output.v2.amountEvery) : 0),
      fieldsEnabled_(mode == Mode::Implicit || config.output.every > 0) {}

bool PicOutputCadence::fields(int step) const noexcept
{ return fieldsEnabled_ && step >= 0 && step % fieldsEvery_ == 0; }
bool PicOutputCadence::periodicSeries(int step) const noexcept
{ return step >= 0 && step % seriesEvery_ == 0; }
bool PicOutputCadence::series(int step) const noexcept
{ return step >= 0 && (periodicSeries(step) || step == finalStep_); }
bool PicOutputCadence::amounts(int step) const noexcept
{ return amountEvery_ > 0 && step >= 0 && (step % amountEvery_ == 0 || step == finalStep_); }
bool PicOutputCadence::explicitDetailFlush(int step) const noexcept
{ return step >= 0 && (step % 50 == 0 || step == finalStep_); }
bool PicOutputCadence::explicitFinalFields() const noexcept
{ return mode_ == Mode::Explicit && fieldsEnabled_ && finalStep_ > 0 && !fields(finalStep_); }

void PicOutputManifest::write(const Input& input, const PicOutputCadence& cadence)
{
    if (cadence.mode() == PicOutputCadence::Mode::Explicit)
        writeManifest(input, "ExplicitPIC", explicitFiles, std::size(explicitFiles),
            {{"fieldsEvery", cadence.fieldsEvery()}, {"seriesEvery", cadence.seriesEvery()},
             {"diagnosticsEvery", cadence.diagnosticsEvery()}});
    else
        writeManifest(input, "ImplicitPIC", implicitFiles, std::size(implicitFiles),
            {{"fieldsEvery", cadence.fieldsEvery()}, {"perfEvery", cadence.seriesEvery()},
             {"physicsEvery", cadence.seriesEvery()}, {"speciesEvery", cadence.seriesEvery()},
             {"amountEvery", cadence.amountEvery()}, {"boundaryEvery", cadence.seriesEvery()},
             {"flushEvery", cadence.seriesEvery()}});
}
