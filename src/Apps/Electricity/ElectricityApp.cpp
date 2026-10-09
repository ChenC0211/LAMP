#include "Electricity/ElectricityApp.h"
#include "Electricity/IO/ElectricityCaseInput.h"
#include "Electricity/IO/ElectricityCaseOutput.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Geometry/InterfaceGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <iostream>

int ElectricityApp::run()
{
    const char* stage="input";
    try
    {
        auto input=readElectricityCaseInput();
        if (config_.checkCase) {
            std::cout<<"[Electricity] Read-only input/capability check passed; no mesh, solve or output.\n";
            return 0;
        }
        stage="fixed-interface";
        QuadtreeMesh mesh(input.domain.mesh);
        ImmersedGeometry<QuadtreeMesh> geometry(input.domain.geometry,mesh);
        geometry.initializeLevelSets(mesh,0.,false);
        const auto& scalar=geometry.levelSets.front();
        std::unique_ptr<InterfaceGeometry> scalarGeometry;
        if(input.domain.electric.interfaceDiscretization=="sharpHermiteP3")
            scalarGeometry=std::make_unique<InterfaceGeometry>(mesh,scalar,input.domain.electric.axisymmetric);
        Quantum::Electricity::DielectricEvolution service(mesh,input.domain.electric,input.domain.diffuseThickness);
        stage="initial-electric-solve";
        auto evaluation=service.evaluate(mesh,scalar,{1,0.,input.run.dt},scalarGeometry.get());
        stage="output";
        const int steps=input.run.fixedStepCount().value();
        ElectricityCaseOutput output(input,steps);
        output.write(mesh,*evaluation,0);
        for(int completed=0;completed<steps;++completed) {
            const int step=completed+1;
            stage="charge-step";
            std::optional<Quantum::Electricity::DielectricChargeUpdate> update;
            if(input.domain.electric.model=="leakyDielectric")
                update=service.completeCharge({step,(step-1)*input.run.dt,input.run.dt});
            stage="coeval-electric-solve";
            // Evaluate accepted Q^(n+1) before publication. EHD keeps the same
            // service's old field until its next force evaluation instead.
            evaluation=service.evaluate(mesh,scalar,{step+1,step*input.run.dt,input.run.dt},scalarGeometry.get());
            stage="accepted-output";
            output.write(mesh,*evaluation,step,update?&*update:nullptr);
            const auto& report=evaluation->electric.report;
            if(step==1 || step==steps || step%input.run.logEverySteps==0)
                std::cout<<"[Electricity] accepted step "<<step<<'/'<<steps
                    <<", time="<<step*input.run.dt<<", residual="<<report.relativeResidual
                    <<", SuperLU="<<report.usedSuperLU<<'\n';
            if(input.run.detailedLog && report.conductionCurrent)
                std::cout<<"[Electricity] coeval boundary current="<<report.conductionCurrent->netBoundaryCurrent<<'\n';
        }
        std::cout<<"[Electricity] Completed "<<steps<<" steps on a fixed interface; output/fields and output/series.\n";
        return 0;
    }
    catch(const std::exception& error) {
        std::cerr<<"[Electricity] failed at "<<stage<<": "<<error.what()<<'\n';
        return std::string(stage)=="input"?2:3;
    }
}
