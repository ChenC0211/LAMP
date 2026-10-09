#pragma once
#include "Modules/Plasma/PlasmaCollisionInputs.h"
#include <string>

// Called from the case working directory. Keep the existing per-block source
// priority: top-level blocks in mcc.dict, then pic.dict, then the separate
// system dictionaries. Nested output/diagnostic blocks are not physics input.
// Writes the existing metadata reports. Null means failure; no partial bundle
// is published. Optional absent dictionaries produce a valid empty bundle.
std::shared_ptr<const PlasmaCollisionInputs> loadPlasmaCollisionInputs(
    const std::string& picPath = "./system/pic.dict");
