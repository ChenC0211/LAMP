#pragma once
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>

struct HydrodynamicsTimeStepReport;
struct HydrodynamicsVofRemapReceipt;
struct VofEhdCsgBoundaryReport;
namespace LAMP::Electric { class ElectricVofChargeMeshTransfer; class ElectricVofChargeTransport; }
struct VofElectrohydrodynamicsOutputContext
{
    std::filesystem::path runDirectory=".";
};
// Passive App-owned serialization of actual accepted owners. No physical
// preparation, field solve, time advancement or accepted-state mutation.
class VofElectrohydrodynamicsOutput
{
public:
    using Context=VofElectrohydrodynamicsOutputContext;
    using Clock=std::chrono::steady_clock;
    explicit VofElectrohydrodynamicsOutput(const VofElectrohydrodynamicsInput&,Context={});
    ~VofElectrohydrodynamicsOutput();
    void begin();
    // A-only observer baseline. Pending completed-step timing is not physical
    // state and is cleared on restore; subsequent normal output starts at N+1.
    std::string checkpoint() const;
    void restoreCheckpoint(const QuadtreeMesh&,const std::string&,
        const LAMP::Hydrodynamics::HydroTwoPhaseState& restored,
        const LAMP::Electric::ElectricVofCharge& restoredCharge);
    void accepted(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const LAMP::Electric::ElectricVofCharge&,const LAMP::Electric::ElectricVofField&,
        const LAMP::Hydrodynamics::HydroVofStep*,Clock::time_point,double,bool final,
        const HydrodynamicsVofRemapReceipt* = nullptr,
        const LAMP::Electric::ElectricVofChargeMeshTransfer* = nullptr,
        // Paired actual physical-End receipts, retained on the old tree after
        // regrid. Null/null preserves old passive consumers; Initial is null.
        const LAMP::Electric::ElectricVofChargeTransport* = nullptr,
        const LAMP::Electric::ElectricVofField* physicalEndField = nullptr,
        const VofEhdCsgBoundaryReport* csgBoundaryPreparation = nullptr);
    void accepted(const QuadtreeMesh&,const LAMP::Hydrodynamics::HydroTwoPhaseState&,
        const LAMP::Electric::ElectricVofSurfaceCharge&,const LAMP::Electric::ElectricVofSharpField&,
        const LAMP::Hydrodynamics::HydroVofStep*,Clock::time_point,double,bool final);
    // Complete controller duration through release/recovery; the publication
    // itself is outside that sampled interval. Requires this exact accepted End.
    void stepCompleted(const LAMP::TimeLayer&,double wallSeconds);
    void timeStep(const HydrodynamicsTimeStepReport&,const LAMP::TimeLayer&,int retries);
    void retry(std::int64_t,std::int64_t,double,double,double,const std::string&);
    void failure(std::int64_t,std::int64_t,double,const std::string&,const std::string&);
    void complete();
private:
    struct Data;
    std::unique_ptr<Data> data_;
};
