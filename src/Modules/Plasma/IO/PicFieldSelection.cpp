#include "Modules/Plasma/IO/PicFieldSelection.h"
#include <algorithm>
#include <cctype>
#include <limits>

namespace
{
using Aliases = std::initializer_list<const char*>;
const Aliases phi = {"phi"}, rho = {"rho"}, n = {"n"};
const Aliases muCell = {"mucell"}, muFace = {"muface"};
const Aliases ex = {"ex", "e", "efield", "electricfield", "electric"};
const Aliases ey = {"ey", "e", "efield", "electricfield", "electric"};
const Aliases psi = {"psi", "levelset", "level_set", "ls", "level-set"};
const Aliases psiComposite = {"psi_composite", "levelset_composite", "level_set_composite",
    "ls_composite", "composite_levelset", "composite_level_set"};
const Aliases psiAll = {"psi_all", "levelset_all", "level_set_all", "ls_all", "all_levelsets", "all_level_sets"};
const Aliases te = {"te", "te_e", "electron_temperature", "electron_temp", "temperature"};
const Aliases meanEnergy = {"electron_mean_energy", "electron_mean_energy_ev", "emean", "mean_energy"};
const Aliases weightMean = {"w", "wmean", "w_mean", "weight", "weights"};
const Aliases weightSum = {"w_sum", "wsum", "weight_sum", "weightsum", "wphys"};
const Aliases weightMax = {"w_max", "wmax", "weight_max"};
const Aliases macroCount = {"n_macro", "macro", "macrocount", "particle", "particlecount",
    "particle_count", "n_particles", "np"};

std::string lower(std::string s)
{
    for (auto& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string fieldKey(const std::string& raw, PicFieldMode mode)
{
    const std::string s = mode == PicFieldMode::Explicit ? lower(raw) : raw;
    std::string mapped;
    mapped.reserve(s.size() + 8);
    for (char ch : s)
    {
        if ((ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9'))
            mapped.push_back(ch);
        else if (ch >= 'A' && ch <= 'Z')
            mapped.push_back(static_cast<char>(ch - 'A' + 'a'));
        else if (mode == PicFieldMode::Explicit && ch == '+')
            mapped += "plus";
        else if (mode == PicFieldMode::Explicit && ch == '-')
            mapped += "minus";
        else if (mode == PicFieldMode::Explicit && ch == '*')
            mapped += "star";
        else
            mapped.push_back('_');
    }
    std::string out;
    out.reserve(mapped.size());
    bool prevUnderscore = false;
    for (char ch : mapped)
    {
        if (ch == '_')
        {
            if (!prevUnderscore) out.push_back(ch);
            prevUnderscore = true;
        }
        else
        {
            out.push_back(ch);
            prevUnderscore = false;
        }
    }
    while (!out.empty() && out.front() == '_') out.erase(out.begin());
    while (!out.empty() && out.back() == '_') out.pop_back();
    return out.empty() ? std::string("unnamed") : out;
}
}

PicFieldKeys::PicFieldKeys(PicFieldMode mode, std::size_t expected) : mode_(mode)
{
    counts_.reserve(expected);
}

std::string PicFieldKeys::next(const std::string& rawName)
{
    const std::string base = fieldKey(rawName, mode_);
    int& count = counts_[base];
    count += 1;
    if (count == 1) return base;
    return base + "_" + std::to_string(count);
}

PicFieldSelection::PicFieldSelection(PicFieldMode mode, const std::vector<std::string>& tokens) : mode_(mode)
{
    static const std::vector<std::string> defaults = {"phi", "rho", "psi", "E", "n", "Te_e", "n_macro"};
    const auto& selected = mode == PicFieldMode::Explicit && tokens.empty() ? defaults : tokens;
    tokens_.reserve(selected.size());
    orderedTokens_.reserve(selected.size());
    for (const auto& raw : selected)
    {
        const std::string token = lower(raw);
        tokens_.insert(token);
        orderedTokens_.push_back(token);
    }
    fields_.phi=has(phi); fields_.rho=has(rho); fields_.psi=has(psi);
    fields_.ex=has(ex); fields_.ey=has(ey); fields_.te=has(te); fields_.meanEnergy=has(meanEnergy);
    fields_.n=mode == PicFieldMode::Explicit || has(n);
    if (mode == PicFieldMode::Implicit)
    {
        fields_.muCell=has(muCell); fields_.muFace=has(muFace);
        fields_.psiComposite=fields_.psi || has(psiComposite); fields_.psiAll=has(psiAll);
    }
    else
    {
        fields_.weightMean=has(weightMean); fields_.weightSum=has(weightSum);
        fields_.weightMax=has(weightMax); fields_.macroCount=has(macroCount);
    }
}

bool PicFieldSelection::has(Aliases aliases) const
{
    for (const char* token : aliases)
        if (tokens_.count(token)) return true;
    return false;
}

int PicFieldSelection::first(Aliases aliases) const
{
    int result = std::numeric_limits<int>::max() / 4;
    for (const char* token : aliases)
        for (std::size_t i = 0; i < orderedTokens_.size(); ++i)
            if (orderedTokens_[i] == token)
            {
                result = std::min(result, static_cast<int>(i));
                break;
            }
    return result;
}

int PicFieldSelection::priority(const std::string& outputName) const
{
    const std::string k = lower(outputName);
    if (k == "phi_si_v") return first(phi);
    if (k == "rho_si_cpm3") return first(rho);
    if (k == "mucell_si_fpm") return first(muCell);
    if (k == "mufaceproxy_si_fpm") return first(muFace);
    if (k == "psi_si_m") return first(psi);
    if (k == "psi_composite_si_m") return first(psiComposite);
    if ((k.rfind("psi_", 0) == 0) && (k != "psi_si_m") && (k != "psi_composite_si_m")) return first(psiAll);
    if (k == "electron_mean_energy_ev") return first(meanEnergy);
    if (k == "te_e_ev") return first(te);
    if (k == "ex_si_vpm") return first(ex);
    if (k == "ey_si_vpm") return first(ey);
    if (k.rfind("n_", 0) == 0) return first(n);
    return std::numeric_limits<int>::max() / 2;
}

std::vector<std::size_t> PicFieldSelection::permutation(const std::vector<std::string>& names) const
{
    std::vector<std::size_t> perm(names.size());
    for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;
    std::stable_sort(perm.begin(), perm.end(), [&](std::size_t a, std::size_t b) {
        const int pa = priority(names[a]);
        const int pb = priority(names[b]);
        if (pa != pb) return pa < pb;
        return a < b;
    });
    return perm;
}

std::string PicFieldSelection::key(const std::string& rawName) const
{
    return fieldKey(rawName, mode_);
}
