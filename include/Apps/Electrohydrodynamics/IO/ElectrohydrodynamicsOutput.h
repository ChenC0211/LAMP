#pragma once
#include "Electrohydrodynamics/IO/ElectrohydrodynamicsInput.h"
#include <chrono>
#include <cstdint>
#include <memory>
namespace LAMP::Electric {
class ElectricSurfaceChargeTimeState;class ElectricPhaseResult;class ElectricPanelCurrentOperator;class ElectricMaxwellTrace;struct ElectricChargeTimeReport;
}
// Required borrow of this App's actual jointly accepted receipts; no second field/clock.
struct ElectrohydrodynamicsOutputView {
    const LAMP::Electric::ElectricSurfaceChargeTimeState& accepted;
    const LAMP::Electric::ElectricPhaseResult& potential;
    const LAMP::Electric::ElectricPanelCurrentOperator& current;
    const LAMP::Electric::ElectricMaxwellTrace& maxwell;
    const LAMP::Electric::ElectricChargeTimeReport* report=nullptr;
};

// Passive accepted-state output owned by this App; no numerical preparation or commit.
class ElectrohydrodynamicsOutput {
public:
    using Clock=std::chrono::steady_clock;
    explicit ElectrohydrodynamicsOutput(const ElectrohydrodynamicsInput&);
    ~ElectrohydrodynamicsOutput();
    void begin();
    void accepted(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const LAMP::Hydrodynamics::HydroTwoPhaseStep*,Clock::time_point,double,bool final,const ElectrohydrodynamicsOutputView&);
    void retry(std::int64_t,std::int64_t,double,double,double,const std::string&);
    void failure(std::int64_t,std::int64_t,double,const std::string&,const std::string&);
    void complete();
private:
    struct Data;std::unique_ptr<Data> data_;
};
