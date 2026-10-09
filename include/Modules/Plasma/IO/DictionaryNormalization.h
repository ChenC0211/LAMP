#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Core/Units/NormalizationContext.h"

class UnitsContext;
struct UnitsScales;

// Resolve legacy manual/plasma-frequency/thermal input reference scales.
// Keeps the original warning/fallback rules and derived-scale reset sentinels;
// the caller constructs UnitsContext and retains geometry closure/report policy.
// Does not read files, mutate inputs or retain references.
UnitsScales resolvePlasmaReferenceScales(const UnitsContext& baseUnits,
    const NormalizationConfig& config, std::string& warning);

// ----------------------------------------------------------------------------
// DictionaryNormalization
//
// Parse optional system/normalization.dict and populate NormalizationConfig.
//
// Supported syntax (block form recommended):
//
// Normalization
// {
//     strategy      manual            // manual | plasmaFrequency | thermal
//     dimension     2                 // 2 | 3
//     slabThickness 1 mm              // only for 2D
//     refSpecies    electron          // optional
//     Te0           3 eV              // optional (thermal)
// }
//
// Keys may also appear at top-level.
// Unknown keys are ignored for forward compatibility.
//
// If the file is missing, readDict() returns false.
// ----------------------------------------------------------------------------

class DictionaryNormalization : public Dictionary
{
public:
    DictionaryNormalization() = default;
    ~DictionaryNormalization() = default;

    // Return false if file does not exist / cannot be opened.
    bool readDict(const std::string& dictPath = "./system/normalization.dict");

    const NormalizationConfig& config() const { return cfg_; }
    bool hasDimension() const { return hasDimension_; }
    bool hasSlabThickness() const { return hasSlabThickness_; }
    bool hasStrategy() const { return hasStrategy_; }
    bool hasAny() const { return hasAny_; }

private:
    NormalizationConfig cfg_;
    bool hasDimension_ = false;
    bool hasSlabThickness_ = false;
    bool hasStrategy_ = false;
    bool hasAny_ = false;

    static std::string toLowerCopy(std::string s);
};
