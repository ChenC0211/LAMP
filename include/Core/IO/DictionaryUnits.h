#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Core/Units/UnitsContext.h"

// ----------------------------------------------------------------------------
// DictionaryUnits
//
// Parse optional system/units.dict and build a UnitsContext.
//
// Supported syntax (block form recommended):
//
// Units
// {
//     mode  SI            // SI | CodeUnits
//     writeSI false        // optional: write code+SI in selected outputs
//     L0    3 cm          // reference length
//     T0    1 ns          // reference time
//     Phi0  500 V         // reference potential
//     n0    1e20 1/m^3     // reference number density (unit token optional)
//     v0    1e5 m/s        // optional (otherwise derived as L0/T0)
//     E0    1e6 V/m        // optional (otherwise derived as Phi0/L0)
//     rho0  1 C/m^3        // optional (defaults to 1 if omitted)
// }
//
// Also accepts single-line key-value pairs without the Units{} block.
//
// If the file is missing, readDict() returns false.
// ----------------------------------------------------------------------------

class DictionaryUnits : public Dictionary
{
public:
    DictionaryUnits() = default;
    ~DictionaryUnits() = default;

    // Return false if file does not exist / cannot be opened.
    bool readDict(const std::string& unitsDictPath = "./system/units.dict");

    const UnitsContext& context() const { return context_; }
    bool writeSIOutputs() const { return writeSIOutputs_; }

private:
    UnitsContext context_ = UnitsContext::MakeCodeUnits();
    bool writeSIOutputs_ = false;

    static std::string toLowerCopy(std::string s);
    static bool parseBoolToken(const std::string& s, bool& out);
};
