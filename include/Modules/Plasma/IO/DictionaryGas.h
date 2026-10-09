#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Modules/Plasma/Physics/BackgroundGas.h"

// ----------------------------------------------------------------------------
// DictionaryGas
//
// Parse optional system/gas.dict and build a GasState (SI).
//
// Supported syntax (block form recommended):
//
// Gas
// {
//     name        Ar
//     pressure    50 mTorr
//     temperature 300 K
// }
//
// Also accepts single-line key-value pairs without the Gas{} block.
//
// If the file is missing, readDict() returns false.
// ----------------------------------------------------------------------------

class DictionaryGas : public Dictionary
{
public:
    DictionaryGas() = default;
    ~DictionaryGas() = default;

    // Return false if file does not exist / cannot be opened.
    bool readDict(const std::string& gasDictPath = "./system/gas.dict");

    const GasState& state() const { return state_; }

private:
    GasState state_;

    static std::string toLowerCopy(std::string s);
};
