#pragma once

#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include <array>
#include <unordered_set>

class DictionaryGeometry;

// Boundary setup consumes prepared records without loading files or retaining
// input references. Returned snapshots own the mappings for the run lifetime.
namespace PicBoundarySetup
{
enum class Mode { Explicit, Implicit };
struct Fields
{
    PicDomainBC domain;
    PicEmbeddedBC embedded;
    PicEmbeddedBCMap embeddedById;
    std::vector<int> lineToBoundaryId;
    // First assigned ID on each side, west/east/south/north; -1 if unassigned.
    std::array<int, 4> domainIds{{-1, -1, -1, -1}};
};
struct Particles
{
    SEEConfig see;
    std::unordered_set<int> embeddedIds;
};
struct Explicit
{
    Fields fields;
    Particles particles;
    bool periodicX = false, periodicY = false;
};

// Both integrators use the same line ownership, side inference, duplicate
// checks and embedded fallback. Mode only selects legacy warning/error text.
Fields mapFields(const PicConfig& config, const DictionaryGeometry& geometry,
    bool haveGeometry, Mode mode);
// Explicit startup retains the additional geometry-pair periodic gate,
// domain SEE sides and original startup logging order/verbosity.
Explicit prepareExplicit(const PicConfig& config, const DictionaryGeometry& geometry, bool haveGeometry);
// Implicit SEE is selected later, after original particle initialization and
// restore. It enables embedded boundaries only and keeps its template policy.
Particles prepareImplicitParticles(const PicConfig& config, const DictionaryGeometry& geometry);
}
