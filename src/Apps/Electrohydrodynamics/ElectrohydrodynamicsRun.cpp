#include "Electrohydrodynamics/ElectrohydrodynamicsRun.h"
#include "Electrohydrodynamics/IO/ElectrohydrodynamicsOutput.h"
#include "Electrohydrodynamics/ElectrohydrodynamicStep.h"
#include "Hydrodynamics/HydrodynamicsRun.h"
#include "Modules/Electric/SurfaceChargeField.h"
#include "Modules/Electric/ElectricPhaseGeometry.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include <chrono>
#include <iostream>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace {
using namespace LAMP::Hydrodynamics;
using namespace LAMP::Interface;
namespace Electric=LAMP::Electric;
namespace Coupled=LAMP::Coupling;

// Joint transaction: one accepted Hydro state and one accepted charge state on the same publication and clock.
class ElectroTransaction final:public HydrodynamicsTransaction {
public:
    ElectroTransaction(const ElectrohydrodynamicsInput& input,ElectrohydrodynamicsOutput& output)
        :input_(input),output_(output){}
    const char* name() const override{return "Electrohydrodynamics";}
    const char* bannerSuffix() const override{return " / sharp-interface electrostatics coupled";}
    const char* inputOkMessage() const override{return "HYDRO_ELECTRIC_INPUT_OK";}
    const char* completionSuffix() const override{return " coupling=sharp-interface electrostatics";}
    const char* prepareStage() const override{return "joint electric/HYDRO AC preparation";}
    const HydroTwoPhaseState* accepted() const override{return coupled_?&coupled_->hydro():nullptr;}
    void initialize(const QuadtreeMesh& mesh,const HermiteGalsPublication& pub,HydroTwoPhaseState initial,
        Clock::time_point start,const char*& stage) override{
        const auto& electric=input_.electric;
        stage="initial electric phase geometry";
        std::cout<<"[INIT] preparing electric geometry on the same FirstJet publication"<<std::endl;
        const auto geometry=Electric::ElectricPhaseGeometry::build(mesh,pub,electric.geometry);
        stage="initial panel charge";
        std::vector<Electric::SurfaceChargeField*> mutableFields;
        std::vector<const Electric::SurfaceChargeField*> chargeViews;
        std::vector<double> initialCharge;
        for(const auto& curve:pub.topology().curves){
            auto field=std::make_unique<Electric::SurfaceChargeField>(mesh,curve,0.);
            std::vector<double> charge;
            for(const auto& panel:field->geometry().panels()){
                long double integral=0;
                for(const auto& node:panel.quadratureNodes)
                    integral+=node.physicalWeight*input_.initialSurfaceChargeDensity(node.frame.point);
                charge.push_back(static_cast<double>(integral));
            }
            field->assignIntegratedCharge(charge);
            initialCharge.insert(initialCharge.end(),charge.begin(),charge.end());
            mutableFields.push_back(field.get());chargeViews.push_back(field.get());
            initialFields_.push_back(std::move(field));
        }
        stage="initial electric/Q solve";
        std::cout<<"[INIT] solving initial electric/Q and binding the joint accepted state"<<std::endl;
        const auto moments=Electric::ElectricSurfaceChargeMomentOperator::prepare(
            geometry,chargeViews,geometry.timeLayer(),electric.moments);
        const auto potential=Electric::solveElectricPhaseElectrostatics(
            geometry,electric.electricProblem,moments.apply(initialCharge),electric.electric);
        const auto solver=Electric::ElectricFixedGeometryChargeSolver::prepare(
            potential,moments,electric.conductivityByPhase);
        const auto current=potential.preparePanelCurrent(moments,electric.conductivityByPhase);
        const auto maxwell=Electric::ElectricMaxwellTrace::prepare(potential);
        maxwell.validateSource(mesh,pub,potential.timeLayer());
        auto chargeState=Electric::ElectricSurfaceChargeTimeState::initialize(
            solver,mutableFields,geometry.timeLayer());
        coupled_.emplace(Coupled::ElectrohydrodynamicState::create(mesh,std::move(initial),std::move(chargeState)));
        const ElectrohydrodynamicsOutputView view{coupled_->electric(),potential,current,maxwell,nullptr};
        stage="accepted output";
        output_.accepted(mesh,coupled_->hydro(),nullptr,start,
            std::chrono::duration<double>(Clock::now()-start).count(),false,view);
    }
    void prepare(const QuadtreeMesh& mesh,const HermiteGalsTimeLayer& target,HydroTwoPhaseStepConfig config,
        bool detailedLog,const char*& stage) override{
        auto jointConfig=input_.electric;jointConfig.hydro=std::move(config);
        if(detailedLog)
            jointConfig.onElectricCandidate=[&](const Electric::ElectricMovingGeometryChargeStep& candidate){
                const auto clock=candidate.step().timeLayer();
                std::cout<<"[ELECTRIC] step="<<clock.step
                    <<" position="<<(clock.position==HermiteGalsTimeLayer::Position::Midpoint?"Midpoint":"End")
                    <<" t="<<clock.time()<<" charge="<<candidate.step().report().relativeTotalEquationResidual<<std::endl;
            };
        step_.emplace(Coupled::ElectrohydrodynamicStep::prepare(mesh,*coupled_,target,jointConfig));
        stage="End Maxwell preparation";
        endMaxwell_.emplace(Electric::ElectricMaxwellTrace::prepare(step_->electricEnd().step().potential()));
        endMaxwell_->validateSource(mesh,step_->hydroStep().candidateState().publication(),target);
    }
    void commit(const QuadtreeMesh& mesh) override{step_->commit(mesh,*coupled_);}
    void output(const QuadtreeMesh& mesh,Clock::time_point start,double compute,
        bool final) override{
        const auto& electric=step_->electricEnd().step();
        const ElectrohydrodynamicsOutputView view{
            coupled_->electric(),electric.potential(),electric.current(),*endMaxwell_,&electric.report()};
        output_.accepted(mesh,coupled_->hydro(),&step_->hydroStep(),start,compute,final,view);
    }
    void releaseStep() override{endMaxwell_.reset();step_.reset();}
    // The initial Fields are borrowed by the accepted charge owner and are destroyed after it.
    void clear() override{endMaxwell_.reset();step_.reset();coupled_.reset();initialFields_.clear();}
private:
    const ElectrohydrodynamicsInput& input_;
    ElectrohydrodynamicsOutput& output_;
    std::vector<std::unique_ptr<Electric::SurfaceChargeField>> initialFields_;
    std::optional<Coupled::ElectrohydrodynamicState> coupled_;
    std::optional<Coupled::ElectrohydrodynamicStep> step_;
    std::optional<Electric::ElectricMaxwellTrace> endMaxwell_;
};
}

int runElectrohydrodynamics(const ElectrohydrodynamicsInput& input,const RunConfig& config){
    ElectrohydrodynamicsOutput output(input);
    ElectroTransaction transaction(input,output);
    const HydrodynamicsRunIO io{
        [&]{output.begin();},
        [&](auto attempted,auto accepted,double time,double oldDt,double newDt,const std::string& reason){
            output.retry(attempted,accepted,time,oldDt,newDt,reason);},
        [&](auto attempted,auto accepted,double time,const char* stage,const std::string& reason){
            output.failure(attempted,accepted,time,stage,reason);},
        [&]{output.complete();}};
    return runHydrodynamics(input.fluid,config,transaction,io);
}
