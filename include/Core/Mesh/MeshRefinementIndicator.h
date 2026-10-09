#pragma once

#include <string>
#include <vector>

struct RefinementIndicatorField
{
    std::string name;              // e.g. "DebyeLength"
    std::vector<double> cellValue; // value per leaf cell (code units)
    double scale = 1.0;            // optional scaling (e.g. C_D for DebyeLength)
    double reference = -1.0;       // optional reference value (e.g. E_ref for ElectricField)
    double percentile = -1.0;      // optional percentile (0-100) if reference not provided
};
