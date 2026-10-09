#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Modules/Plasma/Physics/CollisionsConfig.h"

// -----------------------------------------------------------------------------
// DictionaryCollisions
//
// Parse optional system/collisions.dict and build a CollisionsConfig.
//
// Supported syntax:
//
// Collisions
// {
//   enabled true;
//   method nullCollision;
//
//   // RNG control (for reproducibility, especially under threading):
//   rngMode global; // global|perParticle|counterBased
//   rngSeed 12345;
//
//   // Optional: enable fast-neutral tracking for ion-neutral CX channels when
//   // reactions.dict provides neutral products (e.g., "Ar_fast").
//   trackFastNeutrals true;
//
//   majorantEmax 1000 eV;
//   safetyFactor 1.2;
//   energySamples 2048;
//   clampEnergy true;
//
//   "e-Ar"
//   {
//     projectile e;
//     reactionSet "e-Ar";
//   }
// }
//
// If the file is missing, readDict() returns false.
// -----------------------------------------------------------------------------

class DictionaryCollisions : public Dictionary
{
public:
    DictionaryCollisions() = default;
    ~DictionaryCollisions() = default;

    bool readDict(const std::string& path = "./system/collisions.dict");

    const CollisionsConfig& config() const { return cfg_; }

private:
    CollisionsConfig cfg_;

    static std::string toLowerCopy(std::string s);
    static std::string stripQuotesAndSemicolon(const std::string& s);
};
