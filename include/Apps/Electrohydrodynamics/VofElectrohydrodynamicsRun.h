#pragma once
#include "Hydrodynamics/HydrodynamicsRun.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"
#include "Modules/Electric/ElectricVofField.h"
#include "Modules/Electric/ElectricVofChargeTransport.h"
#include <memory>
#include <array>
#include <cstddef>
#include <cstdint>

class VofElectrohydrodynamicsOutput;
// Actual App CSG boundary-ID preparation only. Scope 0 is Initial/restore,
// 1 is all physical AC/End consumers, and 2 is target-AMR consumers. Separate
// target-domain selection and Core-internal graph construction are excluded.
struct VofEhdCsgBoundaryStageReport
{
    std::uint64_t requests=0,hits=0,rebuilds=0,retainedBuilds=0;
    std::uint64_t selectorCalls=0,completeGraphBuilds=0;
    std::size_t peakRetainedBytes=0;
    bool timingAvailable=true;
    double totalSeconds=0,rebuildSeconds=0;
};
struct VofEhdCsgBoundaryReport
{
    std::array<VofEhdCsgBoundaryStageReport,3> stages{};
    std::size_t byteLimit=32u*1024u*1024u;
};
// One accepted fluid/charge/field owner. The shared Hydro workflow remains
// the only clock/retry owner; module preparation never commits a subset.
class VofElectrohydrodynamicsTransaction final:public HydrodynamicsTransaction
{
public:
    VofElectrohydrodynamicsTransaction(const VofElectrohydrodynamicsInput&,VofElectrohydrodynamicsOutput&);
    ~VofElectrohydrodynamicsTransaction() override;
    const char* name() const override;
    const char* bannerSuffix() const override;
    const char* inputOkMessage() const override;
    const char* completionSuffix() const override;
    const char* prepareStage() const override;
    LAMP::Hydrodynamics::HydroTwoPhaseState makeInitialVofState(const QuadtreeMesh&,
        LAMP::Hydrodynamics::HydroPhaseGeometry,const LAMP::Hydrodynamics::HydroNativeValues&) const override;
    const LAMP::Hydrodynamics::HydroTwoPhaseState* accepted() const override;
    bool supportsCheckpoint() const noexcept override;
    std::string checkpoint(const QuadtreeMesh&) const override;
    void restoreCheckpoint(const QuadtreeMesh&,const std::string&) override;
    const LAMP::Electric::ElectricVofCharge* acceptedCharge() const noexcept;
    const LAMP::Electric::ElectricVofField* acceptedField() const noexcept;
    const LAMP::Electric::ElectricVofSurfaceCharge* acceptedSurfaceCharge() const noexcept;
    const LAMP::Electric::ElectricVofSharpField* acceptedSharpField() const noexcept;
    // Complete joint preparation only; borrowed through releaseStep().
    const LAMP::Hydrodynamics::HydroVofStep* preparedStep() const noexcept;
    // Actual final physical-T receipts, valid after successful prepare until
    // releaseStep, including after a changed commit. Not the target snapshot.
    const LAMP::Electric::ElectricVofChargeTransport* physicalChargeTransport() const noexcept;
    const LAMP::Electric::ElectricVofField* physicalEndField() const noexcept;
    // Borrowed until releaseStep/clear. Null when no CSG helper was requested.
    const VofEhdCsgBoundaryReport* csgBoundaryPreparation() const noexcept;
    void initializeVof(const QuadtreeMesh&,LAMP::Hydrodynamics::HydroTwoPhaseState,
        Clock::time_point,const char*&) override;
    void prepare(const QuadtreeMesh&,const LAMP::TimeLayer&,LAMP::Hydrodynamics::HydroTwoPhaseStepConfig,
        bool,const char*&) override;
    void commit(QuadtreeMesh&) override;
    const HydrodynamicsVofRemapReceipt* remapReceipt() const noexcept override;
    void output(const QuadtreeMesh&,Clock::time_point,double,bool) override;
    void releaseStep() override;
    void clear() override;
private:
    struct Data;
    std::unique_ptr<Data> data_;
};

