#pragma once

#include "Core/Mesh/MeshRefinementIndicator.h"
#include <utility>

// Own the two PIC indicators supplied synchronously to the Core refiner.
// Callers retain their distinct averaging, gating and mesh-refresh policies.
// Empty names mark unavailable indicators; values are never registered globally.
class PicRefinementIndicators
{
public:
    void setDebye(std::vector<double> values, double scale)
    { fields_[0] = {"DebyeLength", std::move(values), scale, -1.0, -1.0}; }
    void setEField(std::vector<double> values, double reference, double percentile = -1.0)
    { fields_[1] = {"ElectricField", std::move(values), 1.0, reference, percentile}; }
    void clearDebye() { clear(fields_[0]); }
    void clearEField() { clear(fields_[1]); }

    const std::vector<double>& debyeValues() const { return fields_[0].cellValue; }
    const std::vector<double>& electricFieldValues() const { return fields_[1].cellValue; }
    const std::vector<RefinementIndicatorField>& fields() const { return fields_; }

private:
    static void clear(RefinementIndicatorField& field)
    { field.name.clear(); field.cellValue.clear(); }
    std::vector<RefinementIndicatorField> fields_{2};
};
