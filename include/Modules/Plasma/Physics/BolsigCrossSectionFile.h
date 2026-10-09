#pragma once

#include <string>
#include <vector>

#include "Modules/Plasma/Physics/CrossSectionDB.h"
#include "Modules/Plasma/Physics/TabulatedCrossSection.h"

// Parser/loader for the standard BOLSIG+ cross-section text format.
// This is also the most common export format for LXCat electron scattering data.
//
// File structure (per process block):
//   Line 1: collision type keyword in CAPITALS (ELASTIC, EXCITATION, IONIZATION, ...)
//   Line 2: target species name (optionally with an arrow and excited state)
//   Line 3: parameter line (type-dependent; may be omitted for ATTACHMENT)
//   Optional: comment lines
//   Table enclosed by dashed lines "-----":
//       E[eV]  sigma[m^2]
//     (or a single line "p  sigma1" denoting sigma(E)=sigma1*E^p)
//
// This loader splits each block into a TabulatedCrossSection and registers it into CrossSectionDB.
class BolsigCrossSectionFile
{
public:
    struct Options
    {
        // Projectile name used when forming DB keys. For LXCat/BOLSIG+ electron sets this is usually "e".
        std::string projectile = "e";

        // Optional prefix prepended to all registered DB keys.
        // This is useful for grouping datasets or avoiding key collisions.
        std::string key_prefix;

        // Options applied to each TabulatedCrossSection created.
        TabulatedCrossSection::Options xs_options;

        // Overwrite existing keys in DB.
        bool overwrite = false;

        // If true, unknown keywords or malformed blocks fail; otherwise they are skipped.
        bool strict = true;
    };

    struct BlockInfo
    {
        std::string key;      // registered DB key
        std::string type;     // keyword (CAPS)
        std::string target;   // target species
        std::string state;    // optional state name (for excitation)
        double parameter = 0.0;      // type-dependent (e.g. me/M for ELASTIC)
        double threshold_eV = 0.0;   // threshold energy for inelastic blocks
        std::vector<std::string> comments;
    };

    // Loads a BOLSIG+ cross-section file and registers all parsed blocks into db.
    // Returns false on error (or if strict=false, returns true but may skip bad blocks).
    static bool LoadAndRegister(const std::string& file_path,
                                CrossSectionDB* db,
                                const Options& options = Options(),
                                std::vector<BlockInfo>* out_blocks = nullptr,
                                std::string* err = nullptr);
};
