#pragma once

#include "Modules/Plasma/PicCaseInput.h"
#include "Modules/Plasma/Physics/PicMeshRefinement.h"
#include "Modules/Plasma/Boundary/PicAmrGeometry.h"
#include <iosfwd>
#include <optional>

class UnitsContext;
class PicSpeciesState;
namespace Solvers { struct SolversConfig; }

// Compatibility IO invoked by the App (or a standalone caller) before withSession().
// Reads the selected PIC case in the current directory, using the supplied
// unit system and solver settings. A missing pic.dict returns no input.
std::optional<PicCaseInput> readPicCaseInput(PicIntegrator integrator,
    const UnitsContext* units, const Solvers::SolversConfig* solvers);

// Capture the legacy implicit PIC environment once, before starting a session.
// Programmatic callers may instead supply typed overrides directly.
ImplicitPicOverrides readImplicitPicOverrides();

// Explicit AMR input has two late gates: the initialized mesh first, then the
// restored particle population. Keep their warnings/output decisions at those
// stages. Own configuration/refiner input; never retain live mesh or particles.
class ExplicitPicAmrSetup
{
public:
    explicit ExplicitPicAmrSetup(const PicCaseInput& input);
    ExplicitPicAmrSetup(const ExplicitPicAmrSetup&) = delete;
    ExplicitPicAmrSetup& operator=(const ExplicitPicAmrSetup&) = delete;
    const AmrConfig& config() const { return config_; }
    const DictionaryMesh& refinerMesh() const { return refinerMesh_; }
    void activateForMesh(int maxLevel, bool geometryRefinementActive, std::ostream& errors);
    void selectElectron(const PicSpeciesState& particles, std::ostream& errors);
    bool enabled() const { return enabled_; }
    int electronIndex() const { return electronIndex_; }
    // Consume the original configured temperature after the particle/output
    // gate; preserve nonpositive fallback and the old nonfinite behavior.
    double temperatureEV(const std::vector<SpeciesConfig>& species) const;
private:
    AmrConfig config_;
    DictionaryMesh refinerMesh_;
    bool meshPhysicsEnabled_ = false, geometryRefinementActive_ = false, enabled_ = false;
    bool meshActivated_ = false, electronSelected_ = false;
    int electronIndex_ = -1;
};

// Resolved from prepared input only: no dictionary/environment reads or live
// mesh/particle references. The actual initialized geometry gate and selected
// electron index are supplied by their owners. Reports stay at the old stage.
struct ImplicitPicAmrSetup
{
    ImplicitPicRefinement::Options indicators;
    PicAmrGeometry::Options geometry;
    int adaptEvery = 0;
    bool fastRemap = false, resetTransients = false;
    void writeReport(std::ostream& output, std::ostream& errors, bool consoleInfo) const;
private:
    friend ImplicitPicAmrSetup resolveImplicitPicAmr(const PicCaseInput&, int, bool);
    bool meshRefreshPresent = false, picRefreshPresent = false, refreshConflict = false;
    // Legacy compatibility report only. Actual exact refresh is selected from
    // the current MovingEB geometry by PicAmrGeometry, independently of this token.
    bool exactMovingRequested = false, exactMovingAllowed = false;
};
ImplicitPicAmrSetup resolveImplicitPicAmr(const PicCaseInput& input,
    int electronIndex, bool geometryRefinementActive);
