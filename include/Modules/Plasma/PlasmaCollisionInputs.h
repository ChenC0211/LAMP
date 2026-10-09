#pragma once
#include <memory>

struct GasState;
class CrossSectionDB;
class ReactionDatabase;
struct CollisionsConfig;

// Immutable service owners shared by PIC and its collision operator.
// A successful vacuum load has four empty owners.
struct PlasmaCollisionInputs
{
    std::shared_ptr<const GasState> gas;
    std::shared_ptr<const CrossSectionDB> crossSections;
    std::shared_ptr<const ReactionDatabase> reactions;
    std::shared_ptr<const CollisionsConfig> collisions;
};
