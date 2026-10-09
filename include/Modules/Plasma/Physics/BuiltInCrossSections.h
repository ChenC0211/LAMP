#pragma once

#include <string>

#include "Modules/Plasma/Physics/CrossSectionDB.h"

namespace BuiltInCrossSections
{

// Register a built-in Argon (Ar) cross-section set (version v1).
//
// Keys are registered as:
//   <prefix><process>
// where process is one of:
//   e-Ar:ELASTIC
//   e-Ar:IONIZATION
//   e-Ar:EXCITATION_1
//   e-Ar:EXCITATION_2
//   e-Ar:EXCITATION_3
//   Ar+-Ar:ELASTIC
//   Ar+-Ar:CHARGE_EXCHANGE
//   Ar+-Ar:IONIZATION
//   Ar+-Ar:EXCITATION_1
//   Ar+-Ar:EXCITATION_2
//   Ar+-Ar:EXCITATION_3
//   ArFast-Ar:ELASTIC
//   ArFast-Ar:IONIZATION
//   ArFast-Ar:EXCITATION_1
//   ArFast-Ar:EXCITATION_2
//   ArFast-Ar:EXCITATION_3
//
// Convention:
//   - Input energy: eV
//   - Output sigma: m^2
bool RegisterArV1(CrossSectionDB* db,
                  const std::string& prefix = "builtin/Ar/v1:",
                  bool overwrite = false,
                  std::string* err = nullptr);

}  // namespace BuiltInCrossSections
