#pragma once

struct ElectrosprayModuleInputs;

#include "Electrospray/Coupling/ElectrosprayInterfaceMotion.h"
#include "Electrospray/Setup/ElectrosprayMaterials.h"
#include "Modules/Hydrodynamics/HydroFittedMeshTransfer.h"
#include "Core/Mesh/MeshRefinementIndicator.h"
#include <functional>
#include <memory>

class QuadtreeMeshRefiner;
class DictionaryGeometry;
class UnitsContext;
struct ElectrosprayElectricState;
struct ElectrosprayConductingState;
namespace Quantum::Electricity { class FittedChargeState; }

// Synchronous read-only observation of one fully accepted coupled state.
// References expire when the observer returns; owning module snapshots may
// be retained through their existing APIs. This view cannot accept a candidate.
struct ElectrosprayAcceptedStateView
{
    const QuadtreeMesh& mesh;
    const HydroFittedDomain& domain;
    const ElectrosprayElectricState& binding;
    const ElectrosprayInterfaceState& interfaceState;
    const Hydrodynamics& hydro;
    const ElectrosprayConductingState& conducting;
    const ElectrosprayMaterialState& materials;
};

// Borrowed, already configured public refiner and geometry definition. No
// dictionary reader, private mesh engine, physical time or solver lives here.
struct ElectrosprayMeshAdaptation
{
    QuadtreeMeshRefiner& refiner;
    const DictionaryGeometry& definition;
    int everySteps = 1, maxSweeps = 2;
    bool currentPhases = false;
    // Configured name/reference only; populated from accepted Hydro each event.
    std::vector<RefinementIndicatorField> flowIndicators;
};

enum class ElectrosprayMeshStage { Grid, InterfaceAndCharge, Flow, ElectricAndForces, Volumes };
using ElectrosprayMeshObserver = std::function<void(ElectrosprayMeshStage, int)>;

struct ElectrosprayMeshReport
{
    int oldCells = 0, newCells = 0, refinedSourceCells = 0, coarsenedSourceCells = 0;
    int protectedCells = 0, sweeps = 0;
    double seconds = 0;
    long double volumeBefore = 0, volumeAfter = 0, volumeRoundoffBound = 0;
    double chargeBefore = 0, chargeAfter = 0;
    long double liquidFluxBefore = 0, liquidFluxAfter = 0;
    HydroFittedMeshTransfer::Report hydro;
};

// A fully prepared TARGET bundle. Physical owners cannot be replaced or
// extracted separately. The source must remain alive and unchanged until
// commit; validateSource also rejects changed owners, mesh and module records.
// Commit exchanges every owner once; the retired bundle lives until destruction.
class ElectrosprayMeshCandidate
{
public:
    ~ElectrosprayMeshCandidate();
    ElectrosprayMeshCandidate(const ElectrosprayMeshCandidate&) = delete;
    ElectrosprayMeshCandidate& operator=(const ElectrosprayMeshCandidate&) = delete;
    bool pending() const noexcept;
    const ElectrosprayMeshReport& report() const;
    const ElectrosprayInterfaceState& targetInterfaceState() const;
    const Hydrodynamics& targetHydro() const;
    const Quantum::Electricity::FittedChargeState& targetElectrical() const;
    void validateSource(const ElectrosprayAcceptedStateView&) const;
    void commit(QuadtreeMesh&, std::unique_ptr<HydroFittedDomain>&,
        std::unique_ptr<ElectrosprayElectricState>&, std::unique_ptr<ElectrosprayInterfaceState>&,
        std::unique_ptr<Hydrodynamics>&, std::unique_ptr<ElectrosprayConductingState>&,
        ElectrosprayMaterialState&);
private:
    struct Data;
    std::unique_ptr<Data> data_;
    explicit ElectrosprayMeshCandidate(const ElectrosprayAcceptedStateView&);
    friend std::unique_ptr<ElectrosprayMeshCandidate> prepareElectrosprayMeshAdaptation(
        ElectrosprayMeshAdaptation&, const ElectrosprayAcceptedStateView&,
        const ElectrosprayModuleInputs&, const UnitsContext&, const ElectrosprayMeshObserver&);
};

// Zero physical time: old adapt + Core Q4 transfer + module charge/history
// transfer + Hydro rebind + existing charged electric/load solvers. Does NOT
// update charge with current, move/reinitialize the surface, write output, or
// accept/replace state. Returns null if no actual topology change occurred.
// Include Core's native scalar-support closure BEFORE target preflight and
// cloning. With no requested change, no mesh/flow/history copy is constructed.
// Accepted Hydro/electric steps, their time levels and the retained scalar
// shared by Interface/electric/loads are checked even on the no-op path.
// Failure leaves source mesh and all accepted physical state untouched; the
// refiner's current-boundary polylines are merely scheduling data.
// Scalar transfer still uses the original acceptance mask and exact polynomial
// gates. A failed candidate is never accepted or retried with different gates.
std::unique_ptr<ElectrosprayMeshCandidate> prepareElectrosprayMeshAdaptation(
    ElectrosprayMeshAdaptation& policy, const ElectrosprayAcceptedStateView& source,
    const ElectrosprayModuleInputs& inputs,
    const UnitsContext& units,
    const ElectrosprayMeshObserver& observer = {});
