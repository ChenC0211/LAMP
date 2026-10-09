#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Modules/Plasma/Physics/CrossSectionDB.h"

// ----------------------------------------------------------------------------
// DictionaryCrossSections
//
// Parse optional system/crossSections.dict and build a CrossSectionDB.
//
// Supported syntax:
//
// CrossSections
// {
//   "e-Ar:ELASTIC"
//   {
//     type table
//     file "data/e_Ar_elastic.txt"
//     interp loglog
//     energyUnit eV
//     sigmaUnit m2
//   }
//
//   "LXCat_Ar"
//   {
//     type lxcatFile
//     file "data/lxcat_Ar_BOLSIG.txt"
//     projectile e
//     interp loglog
//     prefix "e-Ar:"
//   }
//
//   "builtin_Ar_v1"
//   {
//     type builtIn
//     library ArV1
//     prefix "builtin/Ar/v1:"
//     overwrite true
//   }
// }
//
// If the file is missing, readDict() returns false.
// ----------------------------------------------------------------------------

class DictionaryCrossSections : public Dictionary
{
public:
    DictionaryCrossSections() = default;
    ~DictionaryCrossSections() = default;

    bool readDict(const std::string& path = "./system/crossSections.dict");

    const CrossSectionDB& db() const { return db_; }

private:
    CrossSectionDB db_;

    static std::string toLowerCopy(std::string s);
    static std::string stripQuotesAndSemicolon(const std::string& s);
};
