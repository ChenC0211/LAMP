#pragma once

#include <string>

#include "Core/IO/Dictionary.h"
#include "Solvers/IO/KVBlockParser.h"
#include "Solvers/SolversConfig.h"

namespace Solvers
{

// ----------------------------------------------------------------------------
// DictionarySolvers
//
// Parses ./system/solvers.dict.
// The file is intended to carry solver-layer algorithm/linear solver knobs,
// not physics model parameters.
//
// Supported blocks (minimal for now):
//   PoissonCellCentered { ... }
//
// Keys in PoissonCellCentered:
//   useSuperLU  (bool)
//   useCG       (bool)
//   useIC       (bool)
//   maxIter     (int)
//   relTol      (double)
//   verbose     (bool)
//
// Unknown keys are ignored.
// ----------------------------------------------------------------------------
class DictionarySolvers : public Dictionary
{
public:
    DictionarySolvers() = default;
    ~DictionarySolvers() override = default;

    bool readDict(const std::string& path = "./system/solvers.dict");

    const SolversConfig& config() const { return config_; }

private:
    SolversConfig config_;
    IO::KVBlockParser parser_;
};

} // namespace Solvers
