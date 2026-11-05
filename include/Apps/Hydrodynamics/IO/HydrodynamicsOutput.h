#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "Core/Base/TimeLayer.h"
#include <chrono>
#include <cstdint>
#include <memory>
#include <string>

namespace LAMP::Hydrodynamics { class HydroVofStep; }
struct HydrodynamicsTimeStepReport;
struct HydrodynamicsVofRemapReceipt;
// Passive accepted-state output owned by this App; no numerical preparation or commit.
class HydrodynamicsOutput {
public:
    using Clock=std::chrono::steady_clock;
    explicit HydrodynamicsOutput(const HydrodynamicsInput&);
    ~HydrodynamicsOutput();
    void begin();
    // Passive accepted-output baseline, including original volume/work
    // references. Restore only into a fresh begun writer; no physical row
    // or initialization is emitted for the checkpoint clock itself.
    std::string checkpoint() const;
    void restoreCheckpoint(const QuadtreeMesh&,const std::string&,
        const LAMP::Hydrodynamics::HydroTwoPhaseState& restored);
    void accepted(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const LAMP::Hydrodynamics::HydroTwoPhaseStep*,Clock::time_point,double,bool final);
    void acceptedVof(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const LAMP::Hydrodynamics::HydroVofStep*,Clock::time_point,double,bool final,
        const HydrodynamicsVofRemapReceipt* receipt=nullptr);
    void timeStep(const HydrodynamicsTimeStepReport&,const LAMP::TimeLayer&,int retries);
    void retry(std::int64_t,std::int64_t,double,double,double,const std::string&);
    void failure(std::int64_t,std::int64_t,double,const std::string&,const std::string&);
    void complete();
private:
    void acceptedImpl(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const LAMP::Hydrodynamics::HydroTwoPhaseStep*,
        const LAMP::Hydrodynamics::HydroVofStep*,
        Clock::time_point,double,bool,const HydrodynamicsVofRemapReceipt* = nullptr);
    struct Data;std::unique_ptr<Data> data_;
};


