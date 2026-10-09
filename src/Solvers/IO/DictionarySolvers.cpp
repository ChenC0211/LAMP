#include "Solvers/IO/DictionarySolvers.h"

#include <string>

using namespace std;

namespace Solvers
{

bool DictionarySolvers::readDict(const string& path)
{
    string text;
    if (!parser_.loadFile(path, text))
    {
        // solvers.dict is optional; defaults remain.
        return false;
    }

    // -----------------------------
    // PoissonCellCentered block
    // (Compatibility: also accept "Poisson" as the block name.)
    // -----------------------------
    string blockText;
    if (parser_.extractFirstBlockText(text, {"PoissonCellCentered", "Poisson"}, blockText))
    {
        const IO::KVBlock kv = parser_.parseKeyValues(blockText);
        auto& opt = config_.poissonCellCentered;

        opt.useSuperLU = parser_.getBool(kv, "useSuperLU", opt.useSuperLU);
        opt.useCG      = parser_.getBool(kv, "useCG",      opt.useCG);
        opt.useIC      = parser_.getBool(kv, "useIC",      opt.useIC);
        opt.maxIter    = parser_.getInt(kv,  "maxIter",    opt.maxIter);
        opt.relTol     = parser_.getDouble(kv, "relTol",   opt.relTol);
        // Compatibility: relTolerance
        opt.relTol     = parser_.getDouble(kv, "relTolerance", opt.relTol);
        opt.verbose    = parser_.getBool(kv, "verbose",    opt.verbose);

        // Factorization cache controls
        opt.cacheEnabled = parser_.getBool(kv, "cacheEnabled", opt.cacheEnabled);
        opt.cacheRebuildEvery = parser_.getInt(kv, "cacheRebuildEvery", opt.cacheRebuildEvery);
        opt.cacheLogEvery = parser_.getInt(kv, "cacheLogEvery", opt.cacheLogEvery);
        opt.cacheLogOnMiss = parser_.getBool(kv, "cacheLogOnMiss", opt.cacheLogOnMiss);

        // Embedded BC missing policy (error | skip)
        const std::string missPolicy = parser_.getString(kv, "embeddedBCMissingPolicy", "");
        if (!missPolicy.empty())
        {
            const std::string p = IO::KVBlockParser::toLowerCopy(missPolicy);
            if (p == "error" || p == "strict")
                opt.embeddedBCMissingPolicy = PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Error;
            else if (p == "skip" || p == "ignore")
                opt.embeddedBCMissingPolicy = PoissonFVSolverCellCentered::Options::EmbeddedBCMissingPolicy::Skip;
        }

        // Enforce mutually-exclusive solver mode selection for:
        //   useSuperLU / useCG / useIC
        // Priority when multiple are true: SuperLU > useIC(CG+IC) > useCG(CG no-IC).
        const bool wantSuper = opt.useSuperLU;
        const bool wantIC = opt.useIC;
        const bool wantCG = opt.useCG;
        if (wantSuper)
        {
            opt.useSuperLU = true;
            opt.useCG = false;
            opt.useIC = false;
        }
        else if (wantIC)
        {
            // "useIC" mode means iterative CG with IC preconditioner.
            opt.useSuperLU = false;
            opt.useCG = false;
            opt.useIC = true;
        }
        else if (wantCG)
        {
            // "useCG" mode means iterative CG without IC.
            opt.useSuperLU = false;
            opt.useCG = true;
            opt.useIC = false;
        }
        else
        {
            // Safe fallback: direct route.
            opt.useSuperLU = true;
            opt.useCG = false;
            opt.useIC = false;
        }
    }

    return true;
}

} // namespace Solvers
