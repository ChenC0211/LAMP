#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceColumnMode.h"
#include "Modules/Interface/Compatibility/InterfaceMeshTransfer.h"
#include "Modules/Interface/Compatibility/InterfaceNeckCap.h"
#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <fstream>
#include <optional>

class Hydrodynamics;
class HydroFittedVelocitySampler;
struct HydrodynamicsConfig;
struct InterfaceNodalConfig;
struct HydroSolveReport;
namespace Quantum::Hydrodynamics { struct ConservativeCapillaryReport; }
namespace HydroGeometryInput { struct AxisymmetricColumn; }
template<typename T> class VertField;

// App adapter only: owns the one native conservative Interface state and
// supplies existing module services with the case's explicit clock/BCs.
class HydroConservativeInterface
{
public:
    HydroConservativeInterface(QuadtreeMesh&,const VertField<double>* initialScalar,
        const HydrodynamicsConfig&,const InterfaceNodalConfig&,const HydroGeometryInput::AxisymmetricColumn* column=nullptr,
        bool retainHistory=false);
    struct MeshCandidate {
        MeshCandidate();
        ~MeshCandidate();
        MeshCandidate(MeshCandidate&&) noexcept;
        MeshCandidate& operator=(MeshCandidate&&) noexcept=delete;
        // Destruction order retains the mesh until both consumers are gone.
        std::unique_ptr<QuadtreeMesh> mesh;
        std::unique_ptr<Hydrodynamics> hydro;
        std::unique_ptr<HydroConservativeInterface> adapter;
        InterfaceResolutionPlan::Result plan;
        InterfaceMeshTransfer::ConservativeReport currentGeometry,previousGeometry;
        Quantum::Hydrodynamics::HydroConservativeTransport::StateRemapReport flow;
    };
    MeshCandidate prepareAdaptation(Hydrodynamics&,const InterfaceAdaptationControls&) const;
    void inheritOutput(HydroConservativeInterface& old) noexcept;
    const InterfaceState& state() const { return *state_; }
    const NativeLeafDomain& domain() const { return domain_; }
    void materials(Hydrodynamics&,double time) const;
    void materials(Hydrodynamics&,const InterfaceState&,double time) const;
    struct MaterialMotion {
        std::shared_ptr<const InterfaceSurfaceTransfer> physical;
        std::optional<InterfaceNeckCap::MaterialTransfer> event;
    };
    std::unique_ptr<InterfaceState> prepare(const Hydrodynamics&,const InterfaceStepTime&,
        InterfaceConservativeGeometry::AdvanceReport* report=nullptr,
        MaterialMotion* material=nullptr) const;
    void accept(std::unique_ptr<InterfaceState> next) noexcept {
        if(retainHistory_)previous_=std::move(state_);
        state_=std::move(next);
    }
    std::optional<InterfaceColumnMode> columnMode(const InterfaceState&) const;
    void write(const Hydrodynamics&,int step,double time,const std::string& prefix,bool fields,double forceTime,
        const HydroSolveReport* report=nullptr,const Quantum::Hydrodynamics::ConservativeCapillaryReport* capillary=nullptr,
        const InterfaceColumnMode* columnMode=nullptr);
    // Preserve an accepted pre-event frame without appending history again.
    void writeEventSourceFrame(const Hydrodynamics&,const std::string& prefix) const;
private:
    HydroConservativeInterface(QuadtreeMesh&,const HydrodynamicsConfig&,NativeLeafDomain,
        std::unique_ptr<InterfaceState>,std::unique_ptr<InterfaceState>,double columnRadius,
        std::optional<InterfaceNodalConfig::NeckBreakup>,int halfWaves);
    HydroFittedVelocitySampler velocity(const Hydrodynamics&) const;
    Vector2D rootSpacing() const;
    void writeFields(const Hydrodynamics&,int step,double time,const std::string& prefix,double forceTime,
        const Quantum::Hydrodynamics::ConservativeCapillaryReport*) const;
    QuadtreeMesh& mesh_;
    const HydrodynamicsConfig& flow_;
    NativeLeafDomain domain_;
    std::unique_ptr<InterfaceState> state_;
    std::unique_ptr<InterfaceState> previous_;
    bool retainHistory_=false;
    std::vector<int> phases_;
    std::ofstream history_;
    std::ofstream componentHistory_;
    double columnRadius_=0;
    int columnHalfWaves_=1;
    std::optional<InterfaceNodalConfig::NeckBreakup> neckBreakup_;
    std::ofstream columnHistory_;
    std::ofstream eventHistory_;
    std::size_t writtenEvents_=0;
};
