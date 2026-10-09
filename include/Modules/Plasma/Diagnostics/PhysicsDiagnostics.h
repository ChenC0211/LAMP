#pragma once

#include <algorithm>
#include <cctype>
#include <string>
#include <unordered_map>
#include <vector>
#include <utility>


namespace PhysicsDiagnostics
{

inline std::string lowerCopy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

inline std::string sanitizeToken(const std::string& raw)
{
    std::string out;
    out.reserve(raw.size());
    bool prevUnderscore = false;
    for (unsigned char ch : raw)
    {
        if (std::isalnum(ch) != 0)
        {
            out.push_back(static_cast<char>(std::tolower(ch)));
            prevUnderscore = false;
        }
        else if (!prevUnderscore)
        {
            out.push_back('_');
            prevUnderscore = true;
        }
    }
    while (!out.empty() && out.front() == '_')
        out.erase(out.begin());
    while (!out.empty() && out.back() == '_')
        out.pop_back();
    if (out.empty())
        out = "diag";
    return out;
}

inline int resolveEvery(int sharedEvery, int localEvery)
{
    return (localEvery > 0) ? localEvery : sharedEvery;
}

struct ElectrodeCurrentGroup
{
    std::string name;
    std::vector<int> boundaryIds;
};

struct EnergyHistogramSample
{
    std::string name;
    std::string species = "all";
    std::string location = "bulk";
};

// Default EnergyHistogram behavior when the dictionary omits explicit `sample` entries:
// generate one bulk-only histogram per configured species, using `<species>_bulk` names.
inline std::vector<EnergyHistogramSample> defaultBulkSamples(const std::vector<std::string>& speciesNames)
{
    std::vector<EnergyHistogramSample> samples;
    samples.reserve(speciesNames.size());
    for (const auto& speciesName : speciesNames)
    {
        EnergyHistogramSample sample;
        sample.name = sanitizeToken(speciesName) + "_bulk";
        sample.species = speciesName;
        sample.location = "bulk";
        samples.push_back(std::move(sample));
    }
    return samples;
}

struct ElectrodeCurrentBucket
{
    double chargeC = 0.0;
    double electronChargeC = 0.0;
    double ionChargeC = 0.0;
    double neutralCountPhys = 0.0;

    void add(double chargeCIn, double neutralWeight)
    {
        chargeC += chargeCIn;
        if (chargeCIn < 0.0)
            electronChargeC += chargeCIn;
        else if (chargeCIn > 0.0)
            ionChargeC += chargeCIn;
        neutralCountPhys += neutralWeight;
    }
};

// Per-step grouped charge/neutral counts, independent of file availability.
// A repeated positive boundary ID belongs to its first configured group.
class ElectrodeCurrentGroups
{
public:
    explicit ElectrodeCurrentGroups(std::vector<ElectrodeCurrentGroup> groups = {})
        : groups_(std::move(groups)), buckets_(groups_.size())
    {
        for (std::size_t gi = 0; gi < groups_.size(); ++gi)
            for (int boundaryId : groups_[gi].boundaryIds)
                if (boundaryId > 0) boundaryToGroup_.emplace(boundaryId, gi);
    }
    const std::vector<ElectrodeCurrentGroup>& groups() const { return groups_; }
    const std::vector<ElectrodeCurrentBucket>& buckets() const { return buckets_; }
    std::size_t indexOf(int boundaryId) const
    {
        const auto it = boundaryToGroup_.find(boundaryId);
        return it == boundaryToGroup_.end() ? groups_.size() : it->second;
    }
    void beginStep() { std::fill(buckets_.begin(), buckets_.end(), ElectrodeCurrentBucket{}); }
    void record(int boundaryId, double chargeC, double neutralWeight)
    {
        if (boundaryId <= 0) return;
        const auto index = indexOf(boundaryId);
        if (index < buckets_.size()) buckets_[index].add(chargeC, neutralWeight);
    }
private:
    std::vector<ElectrodeCurrentGroup> groups_;
    std::vector<ElectrodeCurrentBucket> buckets_;
    std::unordered_map<int, std::size_t> boundaryToGroup_;
};

} // namespace PhysicsDiagnostics
