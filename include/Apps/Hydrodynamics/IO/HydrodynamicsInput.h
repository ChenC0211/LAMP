#pragma once
#include "IO/CaseDocuments.h"
#include "Core/Base/ScalarExpression.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Modules/Interface/GeometricVofPublication.h"
#include "Modules/Interface/GeometricVofTransport.h"
#include "RunControls.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include <memory>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

class DictionaryGeometry;
namespace FittedDomainGeometry { class InitialPhaseSampler; }
struct HydrodynamicsInput;
namespace AppIO {
struct TwoPhaseInputComponents;
TwoPhaseInputComponents readTwoPhaseInputComponents(HydrodynamicsInput&,const CaseDocuments&,
    ConstantsUsage&,std::optional<int>,int);
}
// Pure fluid typed input. Its standalone reader owns six strict documents.
struct HydrodynamicsInput
{
    // HermiteGals remains a shared-reader/library reference value. The
    // standalone production reader admits Vof explicitly before dispatch.
    enum class InterfaceMethod { HermiteGals, Vof };
    InterfaceMethod interfaceMethod=InterfaceMethod::Vof;
    LAMP::Hydrodynamics::HydroVofAcConfig vof;
    LAMP::Interface::GeometricVofInitializationConfig vofInitialization;
    // App boundary declarations are plain values, independent of a geometry
    // publication or a transported high-order interface state.
    struct VofSymmetryPlane {int normalAxis=0;double coordinate=0;};
    std::vector<VofSymmetryPlane> vofSymmetryPlanes;
    // Explicit straight-S geometric contact law: a physical tangent interval,
    // no mesh/publication/native IDs, pinned radius or flow-patch inference.
    std::optional<Vector2D> vofSouthEven90Interval;

    // Physical side intervals or original Core CSG boundary IDs, never native rows.
    // Empty is the old impermeable case.
    // Only VofElectrohydrodynamics currently parses physical flow patches.
    // Positive sigma requires its separately declared South Even90 geometry.
    struct VofBoundaryPatch {
        enum class Side { W, S, E, N };
        enum class Kind { PrescribedVelocity, PressureOutlet };
        std::string name;
        Side side=Side::W;
        std::optional<int> boundaryId; // Core physical ID; exclusive with side/interval.
        Kind kind=Kind::PrescribedVelocity;
        Vector2D interval{0,0}; // Ordered physical tangent endpoints on base-grid lines.
        LAMP::Interface::GeometricVofAdvectorSlab::Boundary::Mode reservoir=
            LAMP::Interface::GeometricVofAdvectorSlab::Boundary::Mode::PureLiquidReservoir;
        Vector2D velocity{0,0}; // Steady complete vector; normal g is derived.
        double pressure=0; // Physical pressure outlet only, never remap lambda.
    };
    std::vector<VofBoundaryPatch> vofBoundaryPatches;
    // Explicit CSG physical walls; actual positive-area native DOFs are rebound per mesh.
    std::vector<int> vofStationaryNoSlipBoundaryIds;
    // Explicit all-pressure reservoir capability; no fixed normal feed or endpoint vector.
    bool vofPressureDrivenReservoirs=false;
    // Explicit weak incoming momentum on ALL pressure outlets of this flow.
    // Absent keeps the original one-way model and an empty query.
    struct VofPressureOutletBackflow {
        LAMP::Hydrodynamics::HydroVofPressureOutletMomentumModel model=
            LAMP::Hydrodynamics::HydroVofPressureOutletMomentumModel::Ambient;
        std::optional<Vector2D> velocity; // Required only for steady P0 Ambient.
    };
    std::optional<VofPressureOutletBackflow> vofPressureOutletBackflow;
    std::optional<HydrodynamicsVofRegridPolicy> vofRegrid;
    enum class InitialBandScale { BaseCell, CurrentCell, TargetCell };
    InitialBandScale initialBandScale=InitialBandScale::BaseCell;
    bool initialInterfaceBand=false;
    int initialRefinementLevel=1;
    double initialInterfaceBandCells=2;
    bool radial=false;
    Vector2D domainMin{0,0},domainMax{1,1};
    int nx=0,ny=0,negativePhase=1,positivePhase=2,steps=0,characteristicSubsteps=1;
    RunControls run;
    int parallelThreads=2;
    Vector2D initialVelocity{0,0};
    std::string initialPhi;
    // Explicit immutable Core CSG source: no evolving state or substitute jet.
    std::shared_ptr<const DictionaryGeometry> csgGeometry;
    // Owns primitive callbacks/CSG expressions, never a mesh/state reference.
    std::shared_ptr<const FittedDomainGeometry::InitialPhaseSampler> csgSampler;
    std::string csgFile;
    double csgLengthScaleSI=0;
    std::vector<int> csgExcludedPhases;
    struct CsgPin {
        int pointId=0,wetBoundaryId=0,dryBoundaryId=0;
        LAMP::Interface::GeometricVofPublication::PinCertificate::WallTraceLaw wallTrace=
            LAMP::Interface::GeometricVofPublication::PinCertificate::WallTraceLaw::StrictWholeCell;
    };
    std::optional<CsgPin> csgPin;
    LAMP::ScalarExpression::Constants constants;
    LAMP::Hydrodynamics::HydroTwoPhaseStepConfig step;
    DictionaryInput::ConfigurationTrace trace;
    using Document=AppIO::CaseDocument;
    std::vector<Document> documents;
    LAMP::ScalarExpression::Jet initialJet(Vector2D) const;
private:
    std::shared_ptr<const LAMP::ScalarExpression> expression_;
    friend AppIO::TwoPhaseInputComponents AppIO::readTwoPhaseInputComponents(
        HydrodynamicsInput&,const AppIO::CaseDocuments&,AppIO::ConstantsUsage&,std::optional<int>,int);
};
HydrodynamicsInput readHydrodynamicsInput(const std::string& systemDir="system",
    std::optional<int> threadsOverride={});


