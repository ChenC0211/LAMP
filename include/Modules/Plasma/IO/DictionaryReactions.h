#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Modules/Plasma/Physics/ReactionDatabase.h"

class CrossSectionDB;

// ----------------------------------------------------------------------------
// DictionaryReactions
//
// Parse optional system/reactions.dict and build a ReactionDatabase.
// The parser performs strict validation and (optionally) checks that all
// referenced cross-section keys exist in a provided CrossSectionDB.
//
// Supported syntax:
//
// Reactions
// {
//   "e-Ar"
//   {
//     projectile e
//     target     Ar
//     Channels
//     {
//       elastic
//       {
//         crossSection "e-Ar:ELASTIC"
//         scatter isotropic
//       }
//     }
//   }
// }
//
// If the file is missing, readDict() returns false.
// ----------------------------------------------------------------------------

class DictionaryReactions : public Dictionary
{
public:
    DictionaryReactions() = default;
    ~DictionaryReactions() = default;

    bool readDict(const std::string& path = "./system/reactions.dict",
                  const CrossSectionDB* xs_db = nullptr);

    const ReactionDatabase& db() const { return db_; }

private:
    ReactionDatabase db_;

    static std::string toLowerCopy(std::string s);
    static std::string stripQuotesAndSemicolon(const std::string& s);
};
