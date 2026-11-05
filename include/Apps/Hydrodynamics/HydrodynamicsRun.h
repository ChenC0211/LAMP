#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Base/TimeLayer.h"
#include <chrono>
#include <cstdint>
#include <functional>
#include <optional>
#include <stdexcept>
#include <utility>
#include "RunConfig.h"

namespace LAMP::Interface { class GeometricVofMeshTransfer; }
namespace LAMP::Hydrodynamics { class HydroVofMeshTransfer; }
struct HydrodynamicsVofRegridReport;
// Borrowed only until releaseStep(). Null receipt: regrid disabled; no-op:
// decision/physical pointers set, transfer pointers null; changed: all set.
// On changed commit sourceMesh retains the old storage for the actual step.
struct HydrodynamicsVofRemapReceipt {
    const QuadtreeMesh* sourceMesh=nullptr;
    const LAMP::Hydrodynamics::HydroTwoPhaseState* physicalEnd=nullptr;
    const LAMP::Hydrodynamics::HydroVofStep* physicalStep=nullptr;
    const LAMP::Interface::GeometricVofMeshTransfer* interfaceTransfer=nullptr;
    const LAMP::Hydrodynamics::HydroVofMeshTransfer* fluidTransfer=nullptr;
    const HydrodynamicsVofRegridReport* decision=nullptr;
    double inventoryTransferSeconds=0;
    double regridSeconds=0; // Preparation only; commit is also in step_total.
    // Changed-tree preparation only. Empty means no retained target stress
    // owner, no corresponding operation, or unsupported passive timing.
    std::optional<double> targetStressGeometrySeconds;
    std::optional<double> finalJointValidationSeconds;
};

// A transaction supplies an accepted fluid state and prepares/commits one step.
// This fluid workflow owns the single accepted clock, retries and scheduling.
// The calling App binds its own writer; no physical output extension is owned here.
class HydrodynamicsTransaction {
public:
    using Clock=std::chrono::steady_clock;
    virtual ~HydrodynamicsTransaction()=default;
    virtual const char* name() const=0;
    virtual const char* bannerSuffix() const=0;
    virtual const char* inputOkMessage() const=0;
    virtual const char* completionSuffix() const=0;
    virtual const char* prepareStage() const=0;
    // Actual Initial native state from the common geometry/seed. Default is
    // the original closed factory; an App may prepare declared open normals
    // with the original M/Saddle before publishing any accepted state.
    virtual LAMP::Hydrodynamics::HydroTwoPhaseState makeInitialVofState(const QuadtreeMesh& mesh,
        LAMP::Hydrodynamics::HydroPhaseGeometry phase,const LAMP::Hydrodynamics::HydroNativeValues& seed) const {
        return LAMP::Hydrodynamics::HydroTwoPhaseState::createVof(mesh,std::move(phase),seed);
    }
    // Null until initializeVof() has bound the first accepted state.
    virtual const LAMP::Hydrodynamics::HydroTwoPhaseState* accepted() const=0;
    // Complete accepted owner plus passive output baseline, captured while
    // the real previous physical step is still alive. No file IO in Modules.
    virtual bool supportsCheckpoint() const noexcept{return false;}
    virtual std::string checkpoint(const QuadtreeMesh&) const{
        throw std::invalid_argument("this transaction does not support accepted checkpoints");}
    virtual void restoreCheckpoint(const QuadtreeMesh&,const std::string&){
        throw std::invalid_argument("this transaction does not support checkpoint restore");}
    // Retained reference compatibility only; production run rejects non-VOF.
    virtual void initialize(const QuadtreeMesh&,const LAMP::Interface::HermiteGalsPublication&,
        LAMP::Hydrodynamics::HydroTwoPhaseState,Clock::time_point,const char*&){
        throw std::invalid_argument("this transaction does not support Hermite-GALS reference initialization");}
    virtual void initializeVof(const QuadtreeMesh&,LAMP::Hydrodynamics::HydroTwoPhaseState,
        Clock::time_point,const char*&){throw std::invalid_argument("this transaction does not support pure VOF");}
    // May throw HydroTwoPhaseStepRetry; stage names the phase reported on failure.
    virtual void prepare(const QuadtreeMesh&,const LAMP::TimeLayer& target,
        LAMP::Hydrodynamics::HydroTwoPhaseStepConfig,bool detailedLog,const char*& stage)=0;
    // A changed-tree commit swaps prevalidated mesh storage and state handles
    // once. Every throwing preparation/validation precedes the storage swap.
    virtual void commit(QuadtreeMesh&)=0;
    virtual const HydrodynamicsVofRemapReceipt* remapReceipt() const noexcept{return nullptr;}
    virtual void output(const QuadtreeMesh&,Clock::time_point start,double compute,bool final)=0;
    // Drops the prepared candidate once its accepted output is written.
    virtual void releaseStep()=0;
    // Destroys every accepted object before the mesh they borrow is destroyed.
    virtual void clear()=0;
};

struct HydrodynamicsTimeStepReport;
// File lifecycle callbacks are provided by the App that owns the writer.
struct HydrodynamicsRunIO {
    std::function<void()> begin;
    std::function<void(std::int64_t,std::int64_t,double,double,double,const std::string&)> retry;
    std::function<void(std::int64_t,std::int64_t,double,const char*,const std::string&)> failure;
    std::function<void()> complete;
    // Passive adaptive proposal provenance paired with the committed clock.
    std::function<void(const HydrodynamicsTimeStepReport&,const LAMP::TimeLayer&,int)> timeStep;
    // Called once only after successful accepted output, releaseStep and
    // recovery bookkeeping. Duration excludes this callback's own output.
    std::function<void(const LAMP::TimeLayer&,double)> stepCompleted;
};

int runHydrodynamics(const HydrodynamicsInput&,const RunConfig&,HydrodynamicsTransaction&,
    const HydrodynamicsRunIO&);
