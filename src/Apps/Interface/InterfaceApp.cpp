#include "Interface/InterfaceApp.h"
#include "Interface/IO/InterfaceTransportInput.h"
#include "Interface/IO/InterfaceFirstJetOutput.h"
#include "Interface/Setup/InterfaceRunSetup.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Modules/Interface/HermiteGalsMeshTransfer.h"
#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>

namespace
{
using Config=InterfaceTransportConfig;
using namespace LAMP::Interface;
using Position=HermiteGalsTimeLayer::Position;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[Interface FirstJet App] ")+why);}
struct Accepted
{
    std::unique_ptr<QuadtreeMesh> mesh;
    HermiteGalsPublication publication;
    bool amr;
};

int runInterfaceTransport(const InterfaceTransportInput& input,bool checkOnly)
{
    const auto& c=input.config;c.validate();const InterfacePrescribedData data(c);
    std::cout<<std::scientific<<std::setprecision(16)
        <<"[Interface FirstJet] method=hermiteGals transport=FirstJet coordinates="<<(c.radial?"rz":"xy")
        <<" resolution="<<c.resolution<<" dt="<<c.dt<<" steps="<<c.steps<<" initialMesh="<<(c.initialAmr?"amr":"uniform")
        <<" velocity="<<(c.velocity==Config::Velocity::Constant?"constant":"polynomial")
        <<" boundary="<<(c.prescribedBoundary?"prescribed":"none")
        <<" domain="<<c.domainMin.x<<','<<c.domainMin.y<<':'<<c.domainMax.x<<','<<c.domainMax.y
        <<" writeFields="<<(c.writeFields?"true":"false")<<" units=code configurationEntries="<<input.trace.entries.size()<<'\n';
    if(checkOnly){std::cout<<"FIRST_JET_APP_INPUT_OK (no mesh or output)\n";return 0;}
    const char* stage="output startup";int attempted=0;std::unique_ptr<Accepted> accepted;
    InterfaceFirstJetOutput output(input);using Clock=InterfaceFirstJetOutput::Clock;
    try {
        output.begin();stage="initialization";const auto initialStart=Clock::now();
        auto mesh=makeInterfaceTransportMesh(c,c.initialAmr);auto initial=initializeHermiteGalsFirstJet(*mesh,c.radial,[&](Vector2D p){const auto j=data.initial(p);return HermiteGalsBoundaryFirstJet{j.value,{j.derivative,j.partialSecond}};});
        accepted=std::make_unique<Accepted>(Accepted{std::move(mesh),std::move(initial),c.initialAmr});
        const auto initialVolume=accepted->publication.volumes().report().totalRawLiquid;
        require(initialVolume>0 && std::isfinite(initialVolume),"invalid initial geometric volume");
        stage="accepted output";
        output.accepted(*accepted->mesh,accepted->publication,initialStart,
            std::chrono::duration<double>(Clock::now()-initialStart).count(),0);
        long double maximumDrift=0;std::size_t regridCount=0;
        for(int step=1;step<=c.steps;++step){attempted=step;stage="first-jet step";
            const auto stepStart=Clock::now();
            const auto& old=accepted->publication.state();const auto begin=old.timeLayer().time();
            const HermiteGalsTimeLayer time{old.timeLayer().step+1,begin,begin+c.dt,Position::End};
            HermiteGalsCandidate velocityGeometry(*accepted->mesh,time,old.primary());
            const auto velocity=c.velocity==Config::Velocity::Constant
                ?HermiteGalsTemporalVelocitySource::constant(old,velocityGeometry.identity(),time,c.constantVelocity)
                :HermiteGalsTemporalVelocitySource::polynomial(old,velocityGeometry.identity(),time,data.polynomial(begin,time.time()));
            const HermiteGalsTemporalCharacteristicSource characteristic(*accepted->mesh,interfaceTransportTraceDomain(*accepted->mesh,c.radial));
            HermiteGalsBoundaryFirstJetQuery boundary;
            if(c.prescribedBoundary)boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){const auto j=data.boundary(q,t.time());return HermiteGalsBoundaryFirstJet{j.value,{j.derivative,j.partialSecond}};};
            const auto prepared=prepareHermiteGalsFirstJetStep(*accepted->mesh,accepted->publication,characteristic.domain(),velocity.previous(),velocity.arrival(),boundary);
            auto next=prepared.publication();
            const auto stepSeconds=std::chrono::duration<double>(Clock::now()-stepStart).count();
            double transferSeconds=0;
            std::unique_ptr<Accepted> staged;std::optional<HermiteGalsMeshTransferReport> regrid;
            if(std::binary_search(c.regridSteps.begin(),c.regridSteps.end(),step)){
                const auto transferStart=Clock::now();
                stage="mesh transfer";auto target=makeInterfaceTransportMesh(c,!accepted->amr);
                const auto kind=accepted->amr?HermiteGalsMeshTransferKind::FirstJetProjection:HermiteGalsMeshTransferKind::ExactRefinement;
                const auto result=transferHermiteGalsFirstJetMesh(*accepted->mesh,next,*target,kind);regrid=result.report();
                staged=std::make_unique<Accepted>(Accepted{std::move(target),result.publication(),!accepted->amr});
                transferSeconds=std::chrono::duration<double>(Clock::now()-transferStart).count();
            }
            // A scheduled mesh change and its physical step are one App
            // transaction. No clock or mesh is committed on preparation failure.
            if(staged)accepted.swap(staged);else accepted->publication=std::move(next);
            if(regrid){++regridCount;const auto& r=*regrid;
                std::cout<<"FIRST_JET_APP_REGRID,"<<step<<','<<accepted->publication.state().timeLayer().time()<<','
                    <<(r.kind==HermiteGalsMeshTransferKind::ExactRefinement?"exact":"project")<<','<<accepted->amr<<','
                    <<r.refinedCells<<','<<r.coarsenedCells<<','<<r.unchangedCells<<','<<r.relativeDomainVolumeChange<<'\n';}
            const auto volume=accepted->publication.volumes().report().totalRawLiquid,drift=(volume-initialVolume)/initialVolume;
            maximumDrift=std::max(maximumDrift,std::abs(drift));
            std::cout<<"[Step "<<step<<"] t="<<accepted->publication.state().timeLayer().time()<<" dt="<<c.dt
                <<" cells="<<accepted->mesh->cellIndexNum()<<" geometry="<<accepted->publication.state().identity().value()
                <<" V0="<<initialVolume<<" V="<<volume<<" dV/V="<<drift
                <<" kappaMax="<<accepted->publication.report().maxAbsCurvature<<" inventory=NA flux_ledger=NA\n";
            stage="accepted output";
            output.accepted(*accepted->mesh,accepted->publication,stepStart,stepSeconds,transferSeconds,regrid?&*regrid:nullptr,&prepared.report());
        }
        stage="final output";output.complete(maximumDrift);
        std::cout<<"FIRST_JET_APP_RESULT,"<<(c.radial?"rz":"xy")<<','<<c.resolution<<','<<c.dt<<','<<c.steps<<','
            <<accepted->publication.state().timeLayer().time()<<",NA,NA,NA,NA,"<<maximumDrift<<','<<regridCount<<",COMPLETED\n";
        return 0;
    } catch(const std::exception& e) {
        try{output.rejected(attempted,accepted?accepted->publication.state().timeLayer().step:0,
            accepted?accepted->publication.state().timeLayer().time():0,stage,e.what());}
        catch(const std::exception& io){std::cerr<<"[Interface FirstJet] failure diagnostic IO failed: "<<io.what()<<'\n';}
        std::cerr<<"[REJECT] step="<<attempted<<" accepted_step="<<(accepted?accepted->publication.state().timeLayer().step:0)
            <<" accepted_time="<<(accepted?accepted->publication.state().timeLayer().time():0)<<" stage="<<stage<<" reason="<<e.what()<<'\n';
        throw;
    }
}
} // App-local transport workflow

int InterfaceApp::run()
{
    const char* stage="input";
    try {
        const auto input=readInterfaceTransportInput();
        if(!config_.checkCase) stage="hermite-gals-first-jet";
        return runInterfaceTransport(input,config_.checkCase);
    } catch(const std::exception& error) {
        std::cerr<<"[LAMP Interface] stage="<<stage<<"; "<<error.what()<<'\n';
        return std::string(stage)=="input"?2:3;
    }
}
