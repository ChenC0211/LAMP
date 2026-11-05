#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "Modules/Electric/ElectricVofField.h"
#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

// The VOF-EHD reader owns seven strict documents plus a declared raw CSG asset through fluid.
// It is distinct from the retained HG EHD reader and owns no physical state.
struct VofElectrohydrodynamicsInput
{
    enum class Application { General,Electrospray };
    Application application=Application::General; // Explicit entry profile, never inferred from case paths.
    enum class Method { VofFiniteVolume,VofSharpSurface };
    Method method=Method::VofFiniteVolume;
    enum class InitialChargeType { CellDensity,SurfaceChargeProjection };
    InitialChargeType initialChargeType=InitialChargeType::CellDensity; // A-only initialization; B1 keeps its surface owner.
    HydrodynamicsInput fluid; // Sole trace/documents owner, including electric.dict.
    LAMP::Electric::ElectricVofMaterials materials;
    LAMP::Electric::ElectricVofControls controls;
    // A-only, explicit AMR protection of signed volume-charge layers. The
    // indicator has units code_charge/code_volume and uses the same End Q/V.
    std::optional<double> chargeDensityJumpTolerance;
    // W,S,E,N: static constant side data, mapped to actual native faces by Run.
    std::array<LAMP::Electric::ElectricVofBoundary,4> boundaries;
    struct CsgElectricBoundary {
        std::string name;
        int boundaryId=0; // Physical Core definition ID, never a native row.
        LAMP::Electric::ElectricVofBoundary condition;
    };
    // Explicit selected-solid boundary table, mutually exclusive with four sides.
    std::vector<CsgElectricBoundary> csgElectricBoundaries;
    // Same physical-patch order as fluid.vofBoundaryPatches, NOT native rows.
    // Potential incoming entries are steady prescribed inlets or pressure
    // outlets with explicit Hydro Backflow. Their independent signed rhoQ is
    // defaulted/recorded as zero when absent; all other entries are derived zero.
    // Only actual incoming gross volume consumes a configured reservoir value.
    std::vector<double> inletChargeDensityByPatch;
    double initialChargeDensity(Vector2D) const; // A CellDensity: signed code charge / true volume.
    double initialSurfaceChargeDensity(Vector2D) const; // B1 or A explicit projection: charge / actual interface measure.
    LAMP::Interface::VofSmoothGeometryConfig smoothGeometry;
    LAMP::Electric::ElectricSurfaceChargeTransferConfig surfaceTransport;
private:
    std::shared_ptr<const LAMP::ScalarExpression> chargeExpression_;
    friend VofElectrohydrodynamicsInput readVofElectrohydrodynamicsInput(const std::string&,std::optional<int>,Application);
};
// Canonical dispatch/provenance name for the explicit production profile.
const char* vofElectrohydrodynamicsApplicationName(VofElectrohydrodynamicsInput::Application);
VofElectrohydrodynamicsInput readVofElectrohydrodynamicsInput(const std::string& systemDir="system",
    std::optional<int> threadsOverride={},
    VofElectrohydrodynamicsInput::Application=VofElectrohydrodynamicsInput::Application::General);
