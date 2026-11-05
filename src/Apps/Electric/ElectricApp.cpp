#include "Electric/ElectricApp.h"
#include "Electric/IO/ElectricCaseInput.h"
#include "Electric/IO/ElectricCaseOutput.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/TecplotWriter.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Modules/Electric/SurfaceChargeField.h"
#include "Modules/Electric/ElectricPotential.h"
#include "Modules/Electric/ElectricMaxwellTrace.h"
#include <algorithm>
#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

int ElectricApp::run()
{
    using namespace LAMP::Electric;
    using namespace LAMP::Interface;
    using Clock=std::chrono::steady_clock;
    const auto elapsed=[](Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now()-start).count();
    };
    const auto programStart=Clock::now();
    ElectricPerformanceReport initialization;
    auto stageStart=programStart;
    const char* stage="input";
    int attempted=0, acceptedStep=0;
    bool ownsOutput=false;
    double acceptedTime=0;
    // Fields outlive the borrowed P_Q and the sole accepted charge-time owner.
    std::unique_ptr<QuadtreeMesh> mesh;
    std::vector<std::unique_ptr<SurfaceChargeField>> fields;
    std::optional<ElectricSurfaceChargeTimeState> accepted;
    try {
        const auto input=readElectricCaseInput(config_.threadsOverride);
        const int steps=input.run.fixedStepCount().value();
        std::cout<<std::setprecision(10)<<"[Electric] sharp-interface electrostatics / HHO / BE+BDF2, grid="
            <<input.nx<<'x'<<input.ny<<" cells="<<static_cast<long long>(input.nx)*input.ny
            <<" coordinates="<<(input.radial?"rz":"xy")
            <<" dt="<<input.run.dt<<" endTime="<<input.run.endTime
#ifdef _OPENMP
            <<" threads="<<omp_get_max_threads()
#else
            <<" threads=1"
#endif
            <<std::endl;
        if(config_.checkCase) {
            std::cout<<"ELECTRIC_INPUT_OK (read-only schema/capability check; no mesh, solve or output)"<<std::endl;
            return 0;
        }
        initialization.inputSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        stage="output startup";
        ElectricCaseOutput output(input,steps);
        ownsOutput=true;
        initialization.outputStartupSeconds=elapsed(stageStart);
        stage="initial FirstJet publication";
        std::cout<<"[INIT] building native mesh and initial FirstJet publication"<<std::endl;
        DictionaryMesh meshInput;
        meshInput.minCoord=input.domainMin;meshInput.maxCoord=input.domainMax;
        meshInput.Nx=input.nx;meshInput.Ny=input.ny;
        meshInput.maxInitialIterations=0;
        meshInput.maxCellNum=std::max(16384,input.nx*input.ny);
        stageStart=Clock::now();
        mesh=std::make_unique<QuadtreeMesh>(meshInput);
        initialization.meshSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        const auto publication=initializeHermiteGalsFirstJet(*mesh,input.radial,[&](Vector2D point) {
            const auto jet=input.initialJet(point);
            return HermiteGalsBoundaryFirstJet{jet.value,{jet.derivative,jet.partialSecond}};
        },input.publication);
        initialization.publicationSeconds=elapsed(stageStart);
        stage="electric phase geometry";
        std::cout<<"[INIT] preparing electric phase quadrature on the same FirstJet geometry"<<std::endl;
        stageStart=Clock::now();
        const auto geometry=ElectricPhaseGeometry::build(*mesh,publication,input.geometry);
        initialization.geometrySeconds=elapsed(stageStart);
        stageStart=Clock::now();
        stage="initial panel charge";
        std::vector<SurfaceChargeField*> mutableFields;
        std::vector<const SurfaceChargeField*> chargeViews;
        std::vector<double> initialCharge;
        for(const auto& curve:publication.topology().curves) {
            auto field=std::make_unique<SurfaceChargeField>(*mesh,curve,0.);
            std::vector<double> charge;
            for(const auto& panel:field->geometry().panels()) {
                long double integral=0;
                for(const auto& node:panel.quadratureNodes)
                    integral+=node.physicalWeight*input.initialChargeDensity(node.frame.point);
                charge.push_back(static_cast<double>(integral));
            }
            field->assignIntegratedCharge(charge);
            initialCharge.insert(initialCharge.end(),charge.begin(),charge.end());
            mutableFields.push_back(field.get());chargeViews.push_back(field.get());
            fields.push_back(std::move(field));
        }
        initialization.panelChargeSeconds=elapsed(stageStart);
        stage="initial electric/Q solve";
        std::cout<<"[INIT] solving initial electric/Q and preparing the reusable electric factor"<<std::endl;
        stageStart=Clock::now();
        const auto moments=ElectricSurfaceChargeMomentOperator::prepare(geometry,chargeViews,geometry.timeLayer(),input.moments);
        initialization.momentPrepareSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        const auto initial=solveElectricPhaseElectrostatics(geometry,input.electric,moments.apply(initialCharge),input.linear);
        initialization.initialSolveSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        const auto solver=ElectricFixedGeometryChargeSolver::prepare(initial,moments,input.conductivity);
        initialization.chargePrepareSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        const auto& initialCurrent=solver.initialCurrent();
        initialization.initialCurrentSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        const auto initialMaxwell=ElectricMaxwellTrace::prepare(initial,input.maxwell);
        initialMaxwell.validateSource(*mesh,publication,initial.timeLayer());
        initialization.maxwellSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        accepted.emplace(ElectricSurfaceChargeTimeState::initialize(solver,mutableFields,geometry.timeLayer()));
        initialization.commitSeconds=elapsed(stageStart);
        stageStart=Clock::now();
        stage="initial accepted output";
        output.write(*mesh,initial,initialCurrent,initialMaxwell,0);
        initialization.outputSeconds=elapsed(stageStart);
        initialization.totalSeconds=elapsed(programStart);
        output.writePerformance(0,0,initialization,initial);
        std::cout<<"[TIMING INIT] geometry="<<initialization.geometrySeconds
            <<"s electric="<<initialization.initialSolveSeconds
            <<"s output="<<initialization.outputSeconds
            <<"s total="<<initialization.totalSeconds<<"s"<<std::endl;
        std::cout<<"[INIT] accepted step=0 t=0 panels="<<moments.panelCount()
            <<" localGauss="<<initial.report().maxRelativeLocalGaussResidual
            <<" globalGauss="<<initial.report().relativeGlobalGaussResidual<<std::endl;
        for(int step=1;step<=steps;++step) {
            attempted=step;
            const double end=step==steps?input.run.endTime:step*input.run.dt;
            const HermiteGalsTimeLayer target{step,accepted->timeLayer().time(),end,HermiteGalsTimeLayer::Position::End};
            const auto start=Clock::now();
            ElectricPerformanceReport timing;
            stage="BDF2-Schur preparation";
            if(input.run.detailedLog)
                std::cout<<"[ADVANCE] step="<<step<<" accepted_t="<<acceptedTime
                    <<" proposed_dt="<<end-target.begin<<std::endl;
            stageStart=Clock::now();
            const auto candidate=accepted->solver().advance(*accepted,target,input.charge);
            timing.chargeAdvanceSeconds=elapsed(stageStart);
            stageStart=Clock::now();
            stage="coeval Maxwell preparation";
            const auto maxwell=ElectricMaxwellTrace::prepare(candidate.potential(),input.maxwell);
            maxwell.validateSource(*mesh,publication,target);
            timing.maxwellSeconds=elapsed(stageStart);
            stageStart=Clock::now();
            stage="accepted charge commit";
            accepted->commit(candidate);
            acceptedStep=step;acceptedTime=accepted->timeLayer().time();
            timing.commitSeconds=elapsed(stageStart);
            stageStart=Clock::now();
            stage="accepted output";
            output.write(*mesh,candidate.potential(),candidate.current(),maxwell,step,&candidate.report());
            timing.outputSeconds=elapsed(stageStart);
            timing.totalSeconds=elapsed(start);
            output.writePerformance(step,acceptedTime,timing,candidate.potential(),&candidate.report());
            if(step==1 || step==steps || step%input.run.logEverySteps==0) {
                const auto& report=candidate.report();
                const auto& electric=candidate.potential().report();
                long double total=0;for(double value:accepted->integratedCharge())total+=value;
                std::cout<<"[Step "<<step<<"] t="<<acceptedTime<<" dt="<<report.dt
                    <<" order="<<report.temporalOrder<<" Schur="<<report.relativeSchurResidual
                    <<" charge="<<report.relativeTotalEquationResidual
                    <<" localGauss="<<electric.maxRelativeLocalGaussResidual
                    <<" globalGauss="<<electric.relativeGlobalGaussResidual
                    <<" current="<<candidate.current().report().relativeGlobalCurrentBalance
                    <<" Q="<<total
                    <<" iterations="<<report.iterations<<" electricSolves="<<report.electricSolves
                    <<" wall="<<timing.totalSeconds<<"s"
                    <<" (charge="<<timing.chargeAdvanceSeconds<<"s Maxwell="<<timing.maxwellSeconds
                    <<"s output="<<timing.outputSeconds<<"s)"
                    <<std::endl;
            }
        }
        stage="series Tecplot export";
        TecplotWriter::exportCsvHistories("output/series");
        std::cout<<"[COMPLETED] Electric steps="<<acceptedStep<<" time="<<acceptedTime
            <<" method=HHO charge=BE+BDF2 interfaceMotion=fixed"<<std::endl;
        return 0;
    } catch(const std::exception& error) {
        if(ownsOutput){
            try{TecplotWriter::exportCsvHistories("output/series");}
            catch(const std::exception& e){std::cerr<<"[Electric] series Tecplot export failed: "<<e.what()<<'\n';}
            catch(...){std::cerr<<"[Electric] series Tecplot export failed with an unknown error\n";}
        }
        std::cerr<<"[FAILED] Electric stage="<<stage<<" attempted_step="<<attempted
            <<" accepted_step="<<acceptedStep<<" accepted_t="<<acceptedTime
            <<" reason="<<error.what()<<std::endl;
        return std::string(stage)=="input"?2:3;
    }
}
