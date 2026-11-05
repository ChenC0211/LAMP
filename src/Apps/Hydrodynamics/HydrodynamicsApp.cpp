#include "Hydrodynamics/HydrodynamicsApp.h"
#include "Hydrodynamics/HydrodynamicsRun.h"
#include "Hydrodynamics/IO/HydrodynamicsOutput.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Modules/Hydrodynamics/HydroVofMeshTransfer.h"
#include "Core/IO/BinaryArchive.h"
#include <memory>
#include <stdexcept>
#include <chrono>
#include <optional>
#include <utility>

namespace {
using namespace LAMP::Hydrodynamics;

// Single-state fluid transaction: HYDRO A+C over one accepted HydroTwoPhaseState.
class FluidTransaction final:public HydrodynamicsTransaction {
public:
    explicit FluidTransaction(HydrodynamicsOutput& output,const HydrodynamicsInput& input):output_(output),input_(input),vofConfig_(input.vof){}
    const char* name() const override{return "Hydrodynamics";}
    const char* bannerSuffix() const override{return "";}
    const char* inputOkMessage() const override{return "HYDRO_INPUT_OK";}
    const char* completionSuffix() const override{return "";}
    const char* prepareStage() const override{return "HYDRO AC preparation";}
    const HydroTwoPhaseState* accepted() const override{return state_?&*state_:nullptr;}
    bool supportsCheckpoint() const noexcept override{return true;}
    std::string checkpoint(const QuadtreeMesh& mesh) const override{
        if(!state_)throw std::logic_error("fluid checkpoint has no accepted state");
        LAMP::IO::Checkpoint::Writer out("LAMP/Hydrodynamics/transaction/1");
        out.text(state_->checkpoint(mesh,vofStep_?&vofStep_->beginState():nullptr));
        out.text(output_.checkpoint());return out.finish();
    }
    void restoreCheckpoint(const QuadtreeMesh& mesh,const std::string& bytes) override{
        if(state_)throw std::logic_error("fluid checkpoint restore over accepted state");
        LAMP::IO::Checkpoint::Reader in(bytes,"LAMP/Hydrodynamics/transaction/1");
        const auto fluidBytes=in.text(),observer=in.text();in.finish();
        auto restored=std::make_shared<HydroTwoPhaseState>(HydroTwoPhaseState::restoreCheckpoint(mesh,fluidBytes));
        restored->validateFor(mesh);output_.restoreCheckpoint(mesh,observer,*restored);state_.swap(restored);
    }
    void initializeVof(const QuadtreeMesh& mesh,HydroTwoPhaseState initial,Clock::time_point start,const char*& stage) override{
        state_=std::make_shared<HydroTwoPhaseState>(std::move(initial));stage="accepted VOF output";
        output_.acceptedVof(mesh,*state_,nullptr,start,std::chrono::duration<double>(Clock::now()-start).count(),false);
    }
    void prepare(const QuadtreeMesh& mesh,const LAMP::TimeLayer& target,HydroTwoPhaseStepConfig,
        bool,const char*& stage) override{
        // A typed retry starts from the same accepted mesh/state. Drop every
        // previous attempt's fields/transfers before its candidate mesh.
        releaseStep();stage="HYDRO AC preparation";
        vofStep_=HydroVofStep::prepare(mesh,*state_,target,vofConfig_);
        if(!input_.vofRegrid||!hydrodynamicsVofRegridDue(*input_.vofRegrid,target))return;
        const auto& physicalEnd=vofStep_->candidateState();
        stage="VOF regrid tree preparation";const auto regridStart=Clock::now();
        regrid_=prepareHydrodynamicsVofRegrid(mesh,physicalEnd,*input_.vofRegrid);
        if(!regrid_->candidate){
            receipt_={&mesh,&physicalEnd,&*vofStep_,nullptr,nullptr,&regrid_->report};
            receipt_.regridSeconds=std::chrono::duration<double>(Clock::now()-regridStart).count();return;
        }
        auto& grid=*regrid_->candidate;
        stage="VOF regrid inventory transfer";const auto inventoryStart=Clock::now();
        auto transfer=LAMP::Interface::GeometricVofMeshTransfer::prepare(mesh,
            physicalEnd.phaseGeometry().vofPublication(),*grid.mesh,grid.domain,target,vofConfig_.transport.geometry);
        const auto inventorySeconds=std::chrono::duration<double>(Clock::now()-inventoryStart).count();
        stage="VOF regrid fluid transfer";
        HydroVofMeshTransferConfig config;config.momentumDensityRepresentation=vofConfig_.momentumDensityRepresentation;
        remap_=HydroVofMeshTransfer::prepare(mesh,physicalEnd,*grid.mesh,std::move(transfer),
            vofConfig_.liquidDensity,vofConfig_.gasDensity,config);
        remappedState_=std::make_shared<HydroTwoPhaseState>(remap_->candidateState());
        receipt_={grid.mesh.get(),&remap_->sourceState(),&*vofStep_,&remap_->interfaceTransfer(),&*remap_,&regrid_->report};
        receipt_.inventoryTransferSeconds=inventorySeconds;
        receipt_.regridSeconds=std::chrono::duration<double>(Clock::now()-regridStart).count();
    }
    void commit(QuadtreeMesh& mesh) override{
        if(!remap_){vofStep_->commit(mesh,*state_);receiptCommitted_=regrid_.has_value();return;}
        // Commit the physical step into a LOCAL handle first. A stale step or
        // invalid old/target state must fail before touching accepted storage.
        HydroTwoPhaseState physical=*state_;vofStep_->commit(mesh,physical);
        if(!physical.sameIdentity(remap_->sourceState())||!regrid_||!regrid_->candidate||!remappedState_)
            throw std::logic_error("VOF regrid transaction has no matching physical End/candidate");
        auto& grid=*regrid_->candidate;
        remappedState_->validateFor(*grid.mesh);
        if(!sameClock(physical.timeLayer(),remappedState_->timeLayer())
            ||!remappedState_->sameIdentity(remap_->candidateState()))
            throw std::logic_error("VOF regrid transaction changed accepted End/source");
        // swapStorage performs ALL compatibility checks before its nonallocating
        // swaps. Afterwards only noexcept shared_ptr swaps and flags remain.
        mesh.swapStorage(*grid.mesh);
        state_.swap(remappedState_);receiptCommitted_=true;
    }
    const HydrodynamicsVofRemapReceipt* remapReceipt() const noexcept override{
        return receiptCommitted_?&receipt_:nullptr;
    }
    void output(const QuadtreeMesh& mesh,Clock::time_point start,double compute,
        bool final) override{
        const auto* receipt=remapReceipt();
        if(receipt&&receipt->fluidTransfer){
            // D5_REMAP_OUTPUT_RECEIPT: old mesh/step remain alive below.
            output_.acceptedVof(mesh,*state_,&*vofStep_,start,compute,final,receipt);return;
        }
        // D5_NO_REMAP_OUTPUT_RECEIPT: optional null/no-op receipt goes last.
        output_.acceptedVof(mesh,*state_,&*vofStep_,start,compute,final,receipt);
    }
    void releaseStep() override{
        receiptCommitted_=false;receipt_={};
        // remappedState_ holds the OLD accepted state after storage exchange.
        // Destroy every old borrowed field/transfer before grid.mesh (old tree).
        remappedState_.reset();remap_.reset();vofStep_.reset();regrid_.reset();
    }
    void clear() override{releaseStep();state_.reset();}
private:
    static bool sameClock(const LAMP::TimeLayer& a,const LAMP::TimeLayer& b) noexcept{
        return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;
    }
    HydrodynamicsOutput& output_;
    HydrodynamicsInput input_;
    std::shared_ptr<HydroTwoPhaseState> state_,remappedState_;
    std::optional<HydroVofStep> vofStep_;
    std::optional<HydrodynamicsVofRegridDecision> regrid_;
    std::optional<HydroVofMeshTransfer> remap_;
    HydrodynamicsVofRemapReceipt receipt_;
    bool receiptCommitted_=false;
    HydroVofAcConfig vofConfig_;
};
}

int HydrodynamicsApp::run(){
    const auto input=readHydrodynamicsInput("system",config_.threadsOverride);
    HydrodynamicsOutput output(input);
    FluidTransaction transaction(output,input);
    const HydrodynamicsRunIO io{
        [&]{output.begin();},
        [&](auto attempted,auto accepted,double time,double oldDt,double newDt,const std::string& reason){
            output.retry(attempted,accepted,time,oldDt,newDt,reason);},
        [&](auto attempted,auto accepted,double time,const char* stage,const std::string& reason){
            output.failure(attempted,accepted,time,stage,reason);},
        [&]{output.complete();},
        [&](const auto& report,const auto& acceptedClock,int retries){
            output.timeStep(report,acceptedClock,retries);}};
    return runHydrodynamics(input,config_,transaction,io);
}


