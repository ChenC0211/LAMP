#pragma once
#include <memory>

class UnitsContext;
class NormalizationContext;
struct GasState;
class CrossSectionDB;
class ReactionDatabase;
struct CollisionsConfig;

// Physical services used by both PIC drivers. The caller supplies existing
// immutable owners; no Runtime/App selection, CLI or schema state crosses this
// boundary. Prepared case settings are supplied separately to withSession();
// the caller drives macro-steps and PIC retains each step's physical substeps.
struct PlasmaRunInputs
{
    std::shared_ptr<const UnitsContext> units;
    std::shared_ptr<const NormalizationContext> normalization;
    std::shared_ptr<const GasState> gas;
    std::shared_ptr<const CrossSectionDB> crossSections;
    std::shared_ptr<const ReactionDatabase> reactions;
    std::shared_ptr<const CollisionsConfig> collisions;
    bool writeSIOutputs = false;
};
