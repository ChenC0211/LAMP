#include "Modules/Hydrodynamics/HydroVofRefinement.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Hydrodynamics/HydrodynamicsVofRegrid.h"
#include "Modules/Hydrodynamics/HydroTwoPhaseStep.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
using namespace LAMP::Hydrodynamics;using namespace LAMP::Interface;
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(F&& f,const char* why)
{bool rejected=false;try{f();}catch(const std::invalid_argument&){rejected=true;}check(rejected,why);}
std::unique_ptr<QuadtreeMesh> mesh(int n,int maxLevel=1,double width=1){
    DictionaryMesh d;d.minCoord={0,0};d.maxCoord={width,1};d.Nx=d.Ny=n;
    int fullLeaves=n*n;for(int level=0;level<maxLevel;++level)fullLeaves*=4;
    d.maxCellNum=std::max(16384,2*fullLeaves+1024);d.maxInitialIterations=0;
    d.maxLevel_local={maxLevel};return std::make_unique<QuadtreeMesh>(d);
}
long double measure(long double x0,long double x1,long double y0,long double y1,bool rz){
    if(x1<=x0||y1<=y0)return 0;
    return (x1-x0)*(y1-y0)*(rz?std::acos(-1.L)*(x0+x1):1.L);
}
HydroTwoPhaseState state(const QuadtreeMesh& mesh,double center,bool pure=false,bool vortex=false,bool rz=false){
    FaceSegmentGraph graph;graph.build(mesh);auto domain=NativeLeafDomain::build(mesh,graph,rz);
    std::vector<std::array<long double,2>> stocks(mesh.cellIndexNum());
    for(int ci=0;ci<mesh.cellIndexNum();++ci){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);const auto p=cell.centerCoord(),h=cell.size();
        const long double x0=p.x-h.x/2,x1=p.x+h.x/2,y0=p.y-h.y/2,y1=p.y+h.y/2;
        if(pure){stocks[ci]={measure(x0,x1,y0,y1,rz),0};continue;}
        const long double bx0=center-.125L,bx1=center+.125L,by0=.375L,by1=.625L;
        const auto a=std::max(x0,bx0),b=std::min(x1,bx1);
        const auto liquid=measure(a,b,std::max(y0,by0),std::min(y1,by1),rz);
        // Four disjoint gas rectangles are measured independently.
        const auto gas=measure(x0,std::min(x1,bx0),y0,y1,rz)+measure(std::max(x0,bx1),x1,y0,y1,rz)
            +measure(a,b,y0,std::min(y1,by0),rz)+measure(a,b,std::max(y0,by1),y1,rz);
        stocks[ci]={liquid,gas};
    }
    const auto publication=GeometricVofPublication::prepareLocalPlicInventory(mesh,domain,
        {1,0,.125,LAMP::TimeLayer::Position::End},stocks,GeometricVofConfig{});
    const auto layout=HydroDofLayout::build(mesh,publication.domain());
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    std::vector<double> velocity(layout.velocityCount()),pressure(layout.pressureCount());
    const auto psi=[&](Vector2D p)->long double{
        const long double x=p.x,y=p.y;return (rz?x*x*(1-x)*(1-x):x*(1-x))*y*(1-y);
    };
    if(vortex)for(int u=0;u<layout.velocityCount();++u){
        const auto f=domain.face(mesh,layout.nativeFace(u));
        if(f.negativeCell<0||f.positiveCell<0)continue;
        Vector2D a=f.center,b=f.center;
        if(f.horizontal){a.x-=f.tangentSpacing/2;b.x+=f.tangentSpacing/2;}
        else{a.y-=f.tangentSpacing/2;b.y+=f.tangentSpacing/2;}
        velocity[u]=static_cast<double>((f.horizontal?-1:1)*(psi(b)-psi(a))*(rz?2*std::acos(-1.L):1.L)/f.area);
    }
    return HydroTwoPhaseState::createVof(mesh,HydroPhaseGeometry::build(mesh,metrics,publication),
        HydroNativeValues(layout,velocity,pressure));
}
HydroTwoPhaseState zeroState(const QuadtreeMesh& mesh,const GeometricVofPublication& publication){
    const auto layout=HydroDofLayout::build(mesh,publication.domain());
    const auto metrics=HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,layout));
    return HydroTwoPhaseState::createVof(mesh,HydroPhaseGeometry::build(mesh,metrics,publication),
        HydroNativeValues(layout,std::vector<double>(layout.velocityCount()),std::vector<double>(layout.pressureCount())));
}
HydroTwoPhaseState microTails(const QuadtreeMesh& mesh,const HydroTwoPhaseState& reference){
    const auto& old=reference.phaseGeometry().vofPublication();
    std::vector<std::array<long double,2>> stocks;stocks.reserve(old.plicCells().size());
    for(const auto& cell:old.plicCells()){
        auto q=cell.inventory;
        // Independently retained positive phase, never C - majority or alpha.
        if(q[0]==0)q[0]=1e-54L;if(q[1]==0)q[1]=1e-54L;stocks.push_back(q);
    }
    const auto publication=GeometricVofPublication::prepareLocalPlicInventory(mesh,old.domain(),old.timeLayer(),stocks,GeometricVofConfig{});
    return zeroState(mesh,publication);
}
void microInventoryIndicators(){
    HydrodynamicsVofRegridPolicy policy;policy.interfaceBandLayers=1;policy.vorticityTolerance=1000;
    for(bool rz:{false,true})for(bool liquid:{false,true}){
        auto source=mesh(4);const NativeMeshIdentity identity(*source);
        FaceSegmentGraph graph;graph.build(*source);const auto domain=NativeLeafDomain::build(*source,graph,rz);
        std::vector<std::array<long double,2>> inventory(source->cellIndexNum());
        for(int ci=0;ci<source->cellIndexNum();++ci){
            const auto& cell=source->cells(source->cellIndexToLabel[ci]);const auto c=cell.centerCoord(),h=cell.size();
            const auto capacity=measure(c.x-h.x/2,c.x+h.x/2,c.y-h.y/2,c.y+h.y/2,rz);
            inventory[ci]=liquid?std::array<long double,2>{capacity,1e-54L}:std::array<long double,2>{1e-54L,capacity};
        }
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*source,domain,
            {1,0,.125,LAMP::TimeLayer::Position::End},inventory,GeometricVofConfig{});
        for(const auto& cell:publication.plicCells())
            check(!cell.resolved&&cell.measured[liquid?0:1]>0&&cell.measured[liquid?1:0]==0,
                "micro fixture must be unresolved and actual geometric pure");
        const auto accepted=zeroState(*source,publication);const auto savedId=publication.identity();
        const auto savedClock=accepted.timeLayer();const auto savedU=accepted.values().normalVelocity();
        const auto decision=prepareHydrodynamicsVofRegrid(*source,accepted,policy);
        check(decision.report.mixedCells==inventory.size()&&decision.report.unresolvedMixedCells==inventory.size()
            &&decision.report.geometricMixedCells==0,"micro inventory/geometry diagnostic meanings");
        check(!decision.candidate&&!decision.report.requested&&!decision.report.changed&&decision.report.clones==0
            &&decision.report.sweeps==0&&decision.report.adaptationSeconds==0&&decision.report.interfaceSeedCells==0
            &&decision.report.interfaceBandCells==0&&decision.report.pureJumpFaces==0,
            "geometrically pure micro stocks must not trigger band, clone or adaptation");
        identity.validate(*source);accepted.validateFor(*source);
        const auto& after=accepted.phaseGeometry().vofPublication();const auto clock=accepted.timeLayer();
        check(after.identity()==savedId&&clock.step==savedClock.step&&clock.begin==savedClock.begin
            &&clock.end==savedClock.end&&clock.position==savedClock.position&&savedU==accepted.values().normalVelocity(),
            "indicator changed source owner, complete clock or velocity");
        for(std::size_t ci=0;ci<inventory.size();++ci)
            check(after.plicCells()[ci].inventory==inventory[ci],"AMR planning changed independent 1e-54 phase inventory");
    }
    // Pure L/G jumps remain geometric evidence when BOTH sides carry tails.
    auto source=mesh(8);const auto aligned=state(*source,.25);const auto tailJump=microTails(*source,aligned);
    const auto jump=prepareHydrodynamicsVofRegrid(*source,aligned,policy);
    const auto microJump=prepareHydrodynamicsVofRegrid(*source,tailJump,policy);
    check(microJump.candidate&&microJump.report.changed&&microJump.report.pureJumpFaces==jump.report.pureJumpFaces
        &&microJump.report.pureJumpFaces>0&&microJump.report.interfaceSeedCells==jump.report.interfaceSeedCells
        &&microJump.report.interfaceBandCells==jump.report.interfaceBandCells,
        "unresolved majority planes lost the real pure L/G jump or enlarged its band");
    // A resolved curved shape still seeds its measured mixed cells. The
    // reference and tailed variants differ only in independently retained
    // far-field stocks, not in their actual published interface geometry.
    const auto initial=GeometricVofInitialization::prepare(*source,false,[](Vector2D p){
        return std::hypot(p.x-.51,p.y-.49)-.219;
    });
    std::vector<std::array<long double,2>> curveStocks;
    for(const auto& cell:initial.plicCells())curveStocks.push_back(cell.inventory);
    const auto curvePub=GeometricVofPublication::prepareLocalPlicInventory(*source,initial.domain(),
        {1,0,.125,LAMP::TimeLayer::Position::End},curveStocks,GeometricVofConfig{});
    const auto curve=zeroState(*source,curvePub);const auto tailedCurve=microTails(*source,curve);
    const auto clean=prepareHydrodynamicsVofRegrid(*source,curve,policy);
    const auto tailed=prepareHydrodynamicsVofRegrid(*source,tailedCurve,policy);
    check(clean.candidate&&tailed.candidate&&clean.report.geometricMixedCells>0
        &&tailed.report.geometricMixedCells==clean.report.geometricMixedCells
        &&tailed.report.interfaceSeedCells==clean.report.interfaceSeedCells
        &&tailed.report.interfaceBandCells==clean.report.interfaceBandCells,
        "resolved curved interface seed/band changed because of far-field micro stocks");
    check(tailed.report.mixedCells==tailedCurve.phaseGeometry().vofPublication().plicCells().size()
        &&tailed.report.unresolvedMixedCells>clean.report.unresolvedMixedCells,
        "curved fixture must retain separate unresolved stock diagnostics");
}
void balanced(const QuadtreeMesh& mesh,const NativeLeafDomain& d){
    d.validateMeshIdentity(mesh);
    for(const auto f:d.faces(mesh))if(f.negativeCell>=0&&f.positiveCell>=0)
        check(std::abs(mesh.cells(mesh.cellIndexToLabel[f.negativeCell]).level()
            -mesh.cells(mesh.cellIndexToLabel[f.positiveCell]).level())<=1,"native 2:1");
}
void cellMinimumRequirements(){
    using Requirements=HydrodynamicsVofRegridRequirements;
    for(bool rz:{false,true}){
        auto source=mesh(4,3);const NativeMeshIdentity identity(*source);
        const auto accepted=state(*source,.5,true,false,rz);
        const auto& publication=accepted.phaseGeometry().vofPublication();
        const auto savedU=accepted.values().normalVelocity();const auto savedId=publication.identity();
        std::vector<std::array<long double,2>> stocks;
        for(const auto& cell:publication.plicCells())stocks.push_back(cell.inventory);
        HydrodynamicsVofRegridPolicy p;p.maxLevel=3;p.interfaceBandLayers=0;p.vorticityTolerance=1000;
        const auto legacy=prepareHydrodynamicsVofRegrid(*source,accepted,p);
        const auto explicitNull=prepareHydrodynamicsVofRegrid(*source,accepted,p,nullptr);
        std::vector<int> levels(source->cellIndexNum());
        const auto none=Requirements::fromCellMinimumLevels(*source,accepted,levels);
        const auto neutral=prepareHydrodynamicsVofRegrid(*source,accepted,p,&none);
        for(const auto* decision:{&legacy,&explicitNull,&neutral})
            check(!decision->candidate&&!decision->report.requested&&!decision->report.changed
                &&decision->report.clones==0&&decision->report.sweeps==0
                &&decision->report.interfaceBandCells==0&&decision->report.flowRefineCells==0
                &&decision->report.additionalRefineCells==0&&decision->report.additionalHoldCells==0,
                "absent or neutral requirements changed the legacy no-op path");
        int protectedCell=-1;
        for(int ci=0;ci<source->cellIndexNum();++ci){const auto c=source->cells(source->cellIndexToLabel[ci]).centerCoord();
            if(c.x<.25&&c.y<.25)protectedCell=ci;
        }
        check(protectedCell>=0,"local requirement fixture lacks its corner donor");levels[protectedCell]=3;
        const auto local=Requirements::fromCellMinimumLevels(*source,accepted,levels);
        std::fill(levels.begin(),levels.end(),0);
        check(local.minimumLevels()[protectedCell]==3,"requirements borrowed a mutable caller level array");
        const auto clock=local.sourceClock(),actual=accepted.timeLayer();
        check(clock.step==actual.step&&clock.begin==actual.begin&&clock.end==actual.end&&clock.position==actual.position,
            "requirements changed the actual complete source clock");
        const auto decision=prepareHydrodynamicsVofRegrid(*source,accepted,p,&local);
        check(decision.candidate&&decision.report.requested&&decision.report.changed&&decision.report.clones==1
            &&decision.report.interfaceBandCells==0&&decision.report.flowRefineCells==0
            &&decision.report.additionalRefineCells==1&&decision.report.additionalHoldCells==0,
            "a local requirement did not use the real native adaptation path");
        const auto& target=*decision.candidate->mesh;balanced(target,decision.candidate->domain);
        const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(*source,publication.domain(),target,decision.candidate->domain);
        check(!overlap.bySource()[protectedCell].empty(),"protected donor lost target coverage");
        for(const auto index:overlap.bySource()[protectedCell])
            check(target.cells(target.cellIndexToLabel[overlap.pieces()[index].targetCell]).level()==3,
                "a protected donor descendant stayed below its minimum level");
        bool quietCoarse=false;
        for(int ci=0;ci<target.cellIndexNum();++ci){const auto& cell=target.cells(target.cellIndexToLabel[ci]);
            const auto c=cell.centerCoord();if(c.x>.75&&c.y>.75&&cell.level()==0)quietCoarse=true;
        }
        check(quietCoarse&&target.cellIndexNum()<16*64,"a local requirement spuriously refined the whole native domain");
        // A requirement equal to current level is the hysteresis hold request.
        // The unprotected planner must genuinely coarsen this same fine stock.
        const auto refined=state(target,.5,true,false,rz);
        const auto freeCoarsen=prepareHydrodynamicsVofRegrid(target,refined,p);
        check(freeCoarsen.candidate&&freeCoarsen.report.changed,"fine pure fixture did not allow unprotected coarsening");
        const auto noProtection=Requirements::fromCellMinimumLevels(target,refined,std::vector<int>(target.cellIndexNum()));
        const auto neutralCoarsen=prepareHydrodynamicsVofRegrid(target,refined,p,&noProtection);
        check(neutralCoarsen.candidate&&neutralCoarsen.report.changed
            &&neutralCoarsen.candidate->mesh->cellIndexNum()==freeCoarsen.candidate->mesh->cellIndexNum()
            &&neutralCoarsen.report.additionalRefineCells==0&&neutralCoarsen.report.additionalHoldCells==0,
            "neutral requirements changed real unprotected coarsening");
        const NativeMeshIdentity unprotectedIdentity(*freeCoarsen.candidate->mesh);
        unprotectedIdentity.validate(*neutralCoarsen.candidate->mesh);
        std::vector<int> holds(target.cellIndexNum());
        for(int ci=0;ci<target.cellIndexNum();++ci){const auto& cell=target.cells(target.cellIndexToLabel[ci]);const auto c=cell.centerCoord();
            if(c.x<.25&&c.y<.25)holds[ci]=cell.level();
        }
        const auto hold=Requirements::fromCellMinimumLevels(target,refined,holds);
        const auto held=prepareHydrodynamicsVofRegrid(target,refined,p,&hold);
        check(held.report.additionalRefineCells==0&&held.report.additionalHoldCells>0,
            "current-level requirements did not report actual additional coarsen protection");
        const auto& heldMesh=held.candidate?*held.candidate->mesh:target;
        const auto& heldDomain=held.candidate?held.candidate->domain:refined.phaseGeometry().vofPublication().domain();
        balanced(heldMesh,heldDomain);
        const auto heldOverlap=LAMP::Mesh::NativeLeafOverlap::build(target,refined.phaseGeometry().vofPublication().domain(),heldMesh,heldDomain);
        for(const auto& piece:heldOverlap.pieces())
            check(heldMesh.cells(heldMesh.cellIndexToLabel[piece.targetCell]).level()>=holds[piece.sourceCell],
                "a current-level hysteresis request allowed actual coarsening");
        // Same mesh, geometry and clock are insufficient: the actual source
        // Hydro handle is frozen, rather than inferred from array size or IDs.
        const auto otherOwner=zeroState(*source,publication);
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,otherOwner,p,&local);},"foreign actual source owner accepted");
        for(const auto changedClock:{LAMP::TimeLayer{2,0,.125,LAMP::TimeLayer::Position::End},
            LAMP::TimeLayer{1,.01,.125,LAMP::TimeLayer::Position::End},
            LAMP::TimeLayer{1,0,.25,LAMP::TimeLayer::Position::End}}){
            const auto changedPublication=GeometricVofPublication::prepareLocalPlicInventory(*source,publication.domain(),changedClock,stocks,GeometricVofConfig{});
            const auto changed=zeroState(*source,changedPublication);
            rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,changed,p,&local);},"different full source clock accepted");
        }
        auto foreignMesh=mesh(4,3,2);
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*foreignMesh,accepted,p,&local);},"equal-size foreign native geometry accepted");
        auto bad=levels;bad.pop_back();
        rejects([&]{(void)Requirements::fromCellMinimumLevels(*source,accepted,bad);},"incomplete native requirements accepted");
        bad=levels;bad[0]=-1;
        rejects([&]{(void)Requirements::fromCellMinimumLevels(*source,accepted,bad);},"negative native requirement accepted");
        bad[0]=4;
        rejects([&]{(void)Requirements::fromCellMinimumLevels(*source,accepted,bad);},"unallocated native requirement accepted");
        auto narrow=p;narrow.maxLevel=2;
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,accepted,narrow,&local);},"requirement above active maxLevel accepted");
        narrow=p;narrow.minLevel=1;
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,accepted,narrow,&none);},"requirement below active minLevel accepted");
        const auto initialPublication=GeometricVofPublication::prepareLocalPlicInventory(*source,publication.domain(),{},stocks,GeometricVofConfig{});
        const auto initial=zeroState(*source,initialPublication);
        rejects([&]{(void)Requirements::fromCellMinimumLevels(*source,initial,levels);},"Initial requirements masqueraded as accepted End");
        FaceSegmentGraph selectedGraph;selectedGraph.build(*source,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x<.5;});
        const auto selectedDomain=NativeLeafDomain::build(*source,selectedGraph,rz);
        std::vector<std::array<long double,2>> selectedStocks(stocks.size());
        for(int ci:selectedDomain.cellIndexByRow())selectedStocks[ci]=stocks[ci];
        const auto selectedPublication=GeometricVofPublication::prepareLocalPlicInventory(*source,selectedDomain,actual,selectedStocks,GeometricVofConfig{});
        const auto selected=zeroState(*source,selectedPublication);
        const auto selectedRequirements=Requirements::fromCellMinimumLevels(*source,selected,levels);
        selectedRequirements.validateFor(*source,selected,p);
        check(selectedRequirements.minimumLevels()==levels
            &&selectedRequirements.sourceClock().step==actual.step
            &&selectedRequirements.sourceClock().begin==actual.begin
            &&selectedRequirements.sourceClock().end==actual.end
            &&selectedRequirements.sourceClock().position==actual.position,
            "selected requirements changed actual native-indexed payload or complete End clock");
        std::vector<int> excludedRequest=levels;int excluded=-1;
        for(int ci=0;ci<source->cellIndexNum();++ci)if(selectedDomain.cellRowByIndex()[ci]<0){excluded=ci;break;}
        check(excluded>=0,"selected requirements fixture has no actual excluded leaf");excludedRequest[excluded]=1;
        const auto nonminimumExcluded=Requirements::fromCellMinimumLevels(*source,selected,excludedRequest);
        rejects([&]{nonminimumExcluded.validateFor(*source,selected,p);},"excluded cell requirement above policy minLevel accepted");
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,selected,p,&selectedRequirements);},
            "selected planning accepted no explicit target domain provider");
        const HydrodynamicsVofTargetDomainProvider wrongProvider=[rz](const QuadtreeMesh& current){
            FaceSegmentGraph graph;graph.build(current,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x>=.5;});
            return NativeLeafDomain::build(current,graph,rz);
        };
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,selected,p,&selectedRequirements,wrongProvider);},
            "target provider with equal-size foreign source selection accepted");
        const HydrodynamicsVofTargetDomainProvider selectedProvider=[rz](const QuadtreeMesh& current){
            FaceSegmentGraph graph;graph.build(current,[](const QuadtreeMesh&,const Cell& c){return c.centerCoord().x<.5;});
            return NativeLeafDomain::build(current,graph,rz);
        };
        const auto selectedNoop=prepareHydrodynamicsVofRegrid(*source,selected,p,&selectedRequirements,selectedProvider);
        check(!selectedNoop.candidate&&!selectedNoop.report.requested&&selectedNoop.report.clones==0,
            "neutral selected-domain requirements spuriously adapted");
        auto selectedLevels=levels;selectedLevels[protectedCell]=1;
        const auto selectedLocal=Requirements::fromCellMinimumLevels(*source,selected,selectedLevels);
        const auto selectedRefine=prepareHydrodynamicsVofRegrid(*source,selected,p,&selectedLocal,selectedProvider);
        check(selectedRefine.candidate&&selectedRefine.report.changed&&selectedRefine.report.additionalRefineCells==1,
            "selected requirement did not refine the real active donor");
        balanced(*selectedRefine.candidate->mesh,selectedRefine.candidate->domain);
        const auto selectedOverlap=LAMP::Mesh::NativeLeafOverlap::build(*source,selectedDomain,
            *selectedRefine.candidate->mesh,selectedRefine.candidate->domain);
        for(int ci=0;ci<source->cellIndexNum();++ci)
            check((selectedDomain.cellRowByIndex()[ci]>=0)==!selectedOverlap.bySource()[ci].empty(),
                "selected refinement gave excluded solid physical donor overlap");
        for(const auto index:selectedOverlap.bySource()[protectedCell])
            check(selectedRefine.candidate->mesh->cells(selectedRefine.candidate->mesh->cellIndexToLabel[
                selectedOverlap.pieces()[index].targetCell]).level()>=1,"selected protected donor target stayed coarse");
        const int sourceLeaves=source->cellIndexNum();
        const HydrodynamicsVofTargetDomainProvider changedCoverage=[rz,sourceLeaves](const QuadtreeMesh& current){
            FaceSegmentGraph graph;graph.build(current,[&](const QuadtreeMesh&,const Cell& c){
                return current.cellIndexNum()==sourceLeaves?c.centerCoord().x<.5:c.centerCoord().x>=.5;});
            return NativeLeafDomain::build(current,graph,rz);
        };
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,selected,p,&selectedLocal,changedCoverage);},
            "candidate target provider changed selected physical coverage after actual adaptation");
        identity.validate(*source);accepted.validateFor(*source);
        check(publication.identity()==savedId&&accepted.values().normalVelocity()==savedU,"requirements changed the accepted source owner or velocity");
        for(std::size_t ci=0;ci<stocks.size();++ci)
            check(publication.plicCells()[ci].inventory==stocks[ci],"requirements changed independent source phase stocks");
    }
    std::cout<<"HYDRO_VOF_CELL_REQUIREMENTS_PASS scenes=XY,RZ frozen=1 local=level3 hysteresis=hold owner=fullclock selected_provider=1 2:1=1\n";
}
void phaseSpecificFloors(){
    // Floors follow the actual measured majority geometry. Independently
    // retained 1e-54 minority stocks must not promote the whole gas to liquid.
    for(bool rz:{false,true})for(bool liquid:{false,true}){
        auto source=mesh(4,3);const NativeMeshIdentity identity(*source);
        FaceSegmentGraph graph;graph.build(*source);const auto domain=NativeLeafDomain::build(*source,graph,rz);
        std::vector<std::array<long double,2>> inventory(source->cellIndexNum());
        for(int ci=0;ci<source->cellIndexNum();++ci){
            const auto& cell=source->cells(source->cellIndexToLabel[ci]);const auto c=cell.centerCoord(),h=cell.size();
            const auto capacity=measure(c.x-h.x/2,c.x+h.x/2,c.y-h.y/2,c.y+h.y/2,rz);
            inventory[ci]=liquid?std::array<long double,2>{capacity,1e-54L}:std::array<long double,2>{1e-54L,capacity};
        }
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*source,domain,
            {1,0,.125,LAMP::TimeLayer::Position::End},inventory,GeometricVofConfig{});
        const auto accepted=zeroState(*source,publication);const auto savedId=publication.identity();
        HydrodynamicsVofRegridPolicy p;p.maxLevel=3;p.interfaceBandLayers=0;p.vorticityTolerance=1000;
        p.liquidMinLevel=2;p.gasMinLevel=1;
        const auto decision=prepareHydrodynamicsVofRegrid(*source,accepted,p);
        check(decision.candidate&&decision.report.changed&&decision.report.interfaceSeedCells==0
            &&decision.report.interfaceBandCells==0&&decision.report.geometricMixedCells==0
            &&decision.report.mixedCells==inventory.size()&&decision.report.unresolvedMixedCells==inventory.size(),
            "pure majority with minority stocks lost its geometric phase classification");
        const auto& target=*decision.candidate->mesh;const int floor=liquid?2:1;
        check(target.cellIndexNum()==16*(1<<(2*floor)),"minor inventory broadened the phase floor to the other phase");
        for(int ci=0;ci<target.cellIndexNum();++ci)
            check(target.cells(target.cellIndexToLabel[ci]).level()==floor,"actual pure-phase target is not at its own floor");
        balanced(target,decision.candidate->domain);identity.validate(*source);accepted.validateFor(*source);
        check(accepted.phaseGeometry().vofPublication().identity()==savedId,"phase floor planning changed accepted owner");
        for(std::size_t ci=0;ci<inventory.size();++ci)
            check(accepted.phaseGeometry().vofPublication().plicCells()[ci].inventory==inventory[ci],
                "phase floor planning changed independently retained minority inventory");
    }
    // Independent planar geometry oracle: liquid x<.51, gas x>.51. Actual
    // mixed source cells must reach the interface level, far pure phases their
    // own floor; face balance may retain a finer transition but never a deficit.
    for(bool rz:{false,true}){
        auto source=mesh(8,3);const NativeMeshIdentity identity(*source);
        FaceSegmentGraph graph;graph.build(*source);const auto domain=NativeLeafDomain::build(*source,graph,rz);
        std::vector<std::array<long double,2>> inventory(source->cellIndexNum());
        for(int ci=0;ci<source->cellIndexNum();++ci){
            const auto& cell=source->cells(source->cellIndexToLabel[ci]);const auto c=cell.centerCoord(),h=cell.size();
            const long double x0=c.x-h.x/2,x1=c.x+h.x/2,y0=c.y-h.y/2,y1=c.y+h.y/2;
            inventory[ci]={measure(x0,std::min(x1,.51L),y0,y1,rz),measure(std::max(x0,.51L),x1,y0,y1,rz)};
        }
        const auto publication=GeometricVofPublication::prepareLocalPlicInventory(*source,domain,
            {1,0,.125,LAMP::TimeLayer::Position::End},inventory,GeometricVofConfig{});
        const auto accepted=zeroState(*source,publication);const auto savedId=publication.identity();
        HydrodynamicsVofRegridPolicy p;p.maxLevel=3;p.interfaceBandLayers=0;p.vorticityTolerance=1000;
        p.liquidMinLevel=2;p.gasMinLevel=1;
        const auto decision=prepareHydrodynamicsVofRegrid(*source,accepted,p);
        check(decision.candidate&&decision.report.changed&&decision.report.geometricMixedCells==8,
            "planar two-phase floor fixture lost actual mixed interface cells");
        const auto& target=*decision.candidate->mesh;
        const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(*source,domain,target,decision.candidate->domain);
        bool quietLiquid=false,quietGas=false;
        for(const auto& piece:overlap.pieces()){
            const auto& donor=source->cells(source->cellIndexToLabel[piece.sourceCell]);
            const auto c=donor.centerCoord(),h=donor.size();const double x0=c.x-h.x/2,x1=c.x+h.x/2;
            const int required=x1<=.51?2:x0>=.51?1:3;
            check(target.cells(target.cellIndexToLabel[piece.targetCell]).level()>=required,
                "actual overlap descendant is below analytic phase/interface floor");
        }
        for(int ci=0;ci<target.cellIndexNum();++ci){const auto& cell=target.cells(target.cellIndexToLabel[ci]);
            if(cell.centerCoord().x<.125&&cell.level()==2)quietLiquid=true;
            if(cell.centerCoord().x>.875&&cell.level()==1)quietGas=true;
        }
        check(quietLiquid&&quietGas,"phase floors did not retain distinct quiet liquid/gas regions");
        balanced(target,decision.candidate->domain);identity.validate(*source);accepted.validateFor(*source);
        check(accepted.phaseGeometry().vofPublication().identity()==savedId,"two-phase floor planning changed source publication");
        for(std::size_t ci=0;ci<inventory.size();++ci)
            check(accepted.phaseGeometry().vofPublication().plicCells()[ci].inventory==inventory[ci],
                "two-phase floor planning changed source phase stocks");
        auto bad=p;bad.liquidMinLevel=-1;
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,accepted,bad);},"negative liquid floor accepted");
        bad=p;bad.gasMinLevel=4;
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,accepted,bad);},"gas floor above allocation accepted");
        bad=p;bad.minLevel=2;bad.gasMinLevel=1;
        rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,accepted,bad);},"phase floor below common floor accepted");
    }
    std::cout<<"HYDRO_VOF_PHASE_FLOOR_PASS scenes=XY,RZ pure_micro=4 planar=2 liquid=2 gas=1 interface=3\n";
}
void regridInterval(){
    using Time=LAMP::TimeLayer;using Position=Time::Position;
    HydrodynamicsVofRegridPolicy policy;const auto end=[](std::int64_t step){return Time{step,.25,.5,Position::End};};
    check(policy.interval==1&&hydrodynamicsVofRegridDue(policy,end(1)),"default policy lost every accepted End");
    policy.interval=50;
    check(!hydrodynamicsVofRegridDue(policy,end(1))&&!hydrodynamicsVofRegridDue(policy,end(49))
        &&hydrodynamicsVofRegridDue(policy,end(50))&&hydrodynamicsVofRegridDue(policy,end(100)),
        "accepted global-step scheduling lost first/multiple boundaries");
    const auto repeated=end(50);
    check(hydrodynamicsVofRegridDue(policy,repeated)&&hydrodynamicsVofRegridDue(policy,repeated),
        "same candidate clock consumed scheduling state on retry");
    const auto last=std::numeric_limits<std::int64_t>::max();
    check(!hydrodynamicsVofRegridDue(policy,end(last))&&hydrodynamicsVofRegridDue(policy,end(last-7)),
        "large int64 accepted-step boundary was narrowed or converted through floating time");
    policy.interval=7;check(hydrodynamicsVofRegridDue(policy,end(last)),"INT64_MAX lost its exact divisible-by-seven boundary");
    policy.interval=50;
    rejects([&]{(void)hydrodynamicsVofRegridDue(policy,Time{50,.5,.5,Position::Initial});},"Initial clock accepted as an End schedule");
    rejects([&]{(void)hydrodynamicsVofRegridDue(policy,Time{50,.25,.5,Position::Midpoint});},"Midpoint clock accepted as an End schedule");
    for(const auto& bad:std::vector<Time>{{0,.25,.5,Position::End},{50,.5,.25,Position::End},
        {50,.25,std::numeric_limits<double>::infinity(),Position::End},{50,.25,.5,static_cast<Position>(17)}})
        rejects([&]{(void)hydrodynamicsVofRegridDue(policy,bad);},"invalid full candidate clock accepted by schedule");
    for(int interval:{0,-1}){auto bad=policy;bad.interval=interval;
        rejects([&]{(void)hydrodynamicsVofRegridDue(bad,end(50));},"nonpositive typed scheduling interval accepted");}
    // Scheduling belongs to the transaction. The explicit geometry factory
    // still performs its real interface decision even on a non-due End.
    auto source=mesh(8);const NativeMeshIdentity identity(*source);const auto accepted=state(*source,.25);
    const auto savedClock=accepted.timeLayer();const auto savedU=accepted.values().normalVelocity();
    const auto savedId=accepted.phaseGeometry().vofIdentity();
    std::vector<std::array<long double,2>> stocks;
    for(const auto& cell:accepted.phaseGeometry().vofPublication().plicCells())stocks.push_back(cell.inventory);
    policy.interfaceBandLayers=1;policy.vorticityTolerance=1000;
    check(!hydrodynamicsVofRegridDue(policy,savedClock),"step1 fixture unexpectedly due under interval50");
    const auto decision=prepareHydrodynamicsVofRegrid(*source,accepted,policy);
    check(decision.candidate&&decision.report.changed&&decision.report.clones==1
        &&decision.report.pureJumpFaces>0&&decision.report.targetLeaves>decision.report.sourceLeaves,
        "explicit regrid factory was incorrectly suppressed by transaction interval");
    balanced(*decision.candidate->mesh,decision.candidate->domain);identity.validate(*source);accepted.validateFor(*source);
    const auto clock=accepted.timeLayer();
    check(savedId==accepted.phaseGeometry().vofIdentity()&&savedU==accepted.values().normalVelocity()
        &&clock.step==savedClock.step&&clock.begin==savedClock.begin&&clock.end==savedClock.end&&clock.position==savedClock.position
        &&decision.report.sourceClock.step==savedClock.step&&decision.report.sourceClock.begin==savedClock.begin
        &&decision.report.sourceClock.end==savedClock.end&&decision.report.sourceClock.position==savedClock.position,
        "interval decision changed actual source identity/velocity/full clock");
    const auto& publication=accepted.phaseGeometry().vofPublication();
    for(std::size_t i=0;i<stocks.size();++i)check(publication.plicCells()[i].inventory==stocks[i],"interval decision changed independent phase stocks");
    auto bad=policy;bad.interval=0;
    rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,accepted,bad);},"explicit factory accepted invalid interval policy");
    auto foreign=mesh(9); // Different native topology/coordinates, not an identical-content clone.
    rejects([&]{(void)prepareHydrodynamicsVofRegrid(*foreign,accepted,policy);},"interval policy bypassed actual foreign mesh rejection");
    std::cout<<"HYDRO_VOF_REGRID_INTERVAL_PASS global_End_steps=1,49,50,100 retry_same_clock int64_boundary explicit_factory_source_and_balance_gates\n";
}
void movingInterface(){
    auto source=mesh(8);const NativeMeshIdentity identity(*source);
    const auto west=state(*source,.25);const auto savedU=west.values().normalVelocity();
    HydrodynamicsVofRegridPolicy p;p.interfaceBandLayers=1;p.vorticityTolerance=1000;
    auto first=prepareHydrodynamicsVofRegrid(*source,west,p);
    check(first.candidate&&first.report.changed&&first.report.clones==1,"aligned interface did not refine");
    check(first.report.mixedCells==0&&first.report.pureJumpFaces>0,"pure L/G jump not protected");
    check(first.report.targetLeaves>first.report.sourceLeaves,"no real refinement");
    identity.validate(*source);check(savedU==west.values().normalVelocity(),"source U changed");
    auto& refined=*first.candidate->mesh;const NativeMeshIdentity refinedIdentity(refined);
    balanced(refined,first.candidate->domain);
    // An independent accepted shape fixture, not a physical transport claim.
    const auto east=state(refined,.75);auto moved=prepareHydrodynamicsVofRegrid(refined,east,p);
    check(moved.candidate&&moved.report.changed,"moved shape did not adapt");
    const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(refined,east.phaseGeometry().vofPublication().domain(),
        *moved.candidate->mesh,moved.candidate->domain);
    bool refine=false,coarsen=false;
    for(const auto& q:overlap.pieces()){
        refine|=q.kind==LAMP::Mesh::NativeLeafOverlapPiece::Kind::Refined;
        coarsen|=q.kind==LAMP::Mesh::NativeLeafOverlapPiece::Kind::Coarsened;
    }
    check(refine&&coarsen,"moving shape needs actual refine/coarsen");
    refinedIdentity.validate(refined);identity.validate(*source);balanced(*moved.candidate->mesh,moved.candidate->domain);
    auto bad=p;bad.coarsenRatio=.5;
    rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,west,bad);},"invalid ratio accepted");
    bad=p;bad.maxLevel=2;
    rejects([&]{(void)prepareHydrodynamicsVofRegrid(*source,west,bad);},"unallocated level accepted");
}
void closedVortex(bool rz){
    auto source=mesh(4);const NativeMeshIdentity identity(*source);
    const auto accepted=state(*source,.5,true,true,rz);const auto savedU=accepted.values().normalVelocity();
    const auto indicator=HydroVofRefinement::sample(*source,accepted);
    check(indicator.vorticity.size()==static_cast<std::size_t>(source->cellIndexNum())
        &&indicator.maximumHAbsVorticity>1e-5,"closed vortex mislabeled curl-free by RT0");
    for(int ci=0;ci<source->cellIndexNum();++ci){
        const auto size=source->cells(source->cellIndexToLabel[ci]).size();
        const auto expected=std::min(size.x,size.y)*std::abs(indicator.vorticity[ci]);
        check(std::abs(indicator.hAbsVorticity[ci]-expected)<=1e-13*std::max(1.,expected),"native h omega units/order");
    }
    HydrodynamicsVofRegridPolicy p;p.vorticityTolerance=1e-5;
    auto refined=prepareHydrodynamicsVofRegrid(*source,accepted,p);
    check(refined.candidate&&refined.report.flowRefineCells>0&&refined.report.interfaceBandCells==0,"pure vortex flow refinement");
    balanced(*refined.candidate->mesh,refined.candidate->domain);
    p.vorticityTolerance=1e6;auto none=prepareHydrodynamicsVofRegrid(*source,accepted,p);
    check(!none.candidate&&none.report.clones==0&&none.report.sweeps==0&&!none.report.requested
        &&!none.report.changed&&none.report.adaptationSeconds==0,"no-op must have zero clone/adapt/transfer/solve");
    identity.validate(*source);check(savedU==accepted.values().normalVelocity(),"source vortex U changed");
}
// Deep ROI targets must survive the real coarse-neighbor prerequisites. This
// checks the App's actual grid planner; it does not simulate transfer or Run init.
void refinementPrerequisites(){
    std::size_t totalRequiredTargets=0;
    for(bool rz:{false,true}){
        const int n=rz?8:4,maxLevel=4;
        auto source=mesh(n,maxLevel);const NativeMeshIdentity sourceIdentity(*source);
        const auto accepted=state(*source,.5,false,false,rz);
        const auto& publication=accepted.phaseGeometry().vofPublication();
        const auto& domain=publication.domain();const auto savedId=publication.identity();
        const auto savedClock=accepted.timeLayer();const auto savedU=accepted.values().normalVelocity();
        std::vector<std::array<long double,2>> savedStocks;
        std::vector<bool> seed(publication.plicCells().size(),false);
        std::vector<int> pure(seed.size(),-1);
        for(std::size_t ci=0;ci<seed.size();++ci){
            const auto& cell=publication.plicCells()[ci];savedStocks.push_back(cell.inventory);
            seed[ci]=cell.resolved&&cell.measured[0]>0&&cell.measured[1]>0;
            if(cell.measured[0]>0&&cell.measured[1]==0)pure[ci]=0;
            else if(cell.measured[1]>0&&cell.measured[0]==0)pure[ci]=1;
        }
        // Independently enumerate genuine shared-face phase jumps and exactly
        // one face-neighbor halo, without calling the planner's mask helper.
        for(const auto face:domain.faces(*source))
            if(face.negativeCell>=0&&face.positiveCell>=0&&pure[face.negativeCell]>=0
                &&pure[face.positiveCell]>=0&&pure[face.negativeCell]!=pure[face.positiveCell])
                seed[face.negativeCell]=seed[face.positiveCell]=true;
        auto required=seed;
        for(const auto face:domain.faces(*source))if(face.negativeCell>=0&&face.positiveCell>=0){
            if(seed[face.negativeCell])required[face.positiveCell]=true;
            if(seed[face.positiveCell])required[face.negativeCell]=true;
        }
        const auto seeds=std::count(seed.begin(),seed.end(),true);
        const auto band=std::count(required.begin(),required.end(),true);
        check(seeds>0&&band>0&&band<source->cellIndexNum(),"deep ROI fixture needs a genuine local geometric band and coarse exterior");
        HydrodynamicsVofRegridPolicy policy;policy.minLevel=0;policy.maxLevel=maxLevel;
        policy.interfaceBandLayers=1;policy.vorticityTolerance=1000;
        const auto decision=prepareHydrodynamicsVofRegrid(*source,accepted,policy);
        check(decision.candidate&&decision.report.requested&&decision.report.changed&&decision.report.clones==1,
            "real deep ROI planner must return an actual changed candidate");
        check(decision.report.interfaceSeedCells==static_cast<std::size_t>(seeds)
            &&decision.report.interfaceBandCells==static_cast<std::size_t>(band),
            "independent geometric source seed/halo disagrees with planned band");
        const auto& target=*decision.candidate->mesh;
        const auto overlap=LAMP::Mesh::NativeLeafOverlap::build(*source,domain,target,decision.candidate->domain);
        overlap.validateFor(*source,target);balanced(target,decision.candidate->domain);
        std::vector<bool> checked(target.cellIndexNum(),false);std::size_t requiredTargets=0;
        const double expectedSpacing=std::ldexp(1./n,-maxLevel);
        for(std::size_t ci=0;ci<required.size();++ci)if(required[ci]){
            check(!overlap.bySource()[ci].empty(),"required source ROI lost actual target coverage");
            for(const auto pieceIndex:overlap.bySource()[ci]){
                check(pieceIndex<overlap.pieces().size(),"required source overlap index outside actual table");
                const auto& piece=overlap.pieces()[pieceIndex];
                check(piece.sourceCell==static_cast<int>(ci)&&piece.targetCell>=0
                    &&piece.targetCell<target.cellIndexNum(),"required source overlap has foreign actual donor/target");
                const auto& cell=target.cells(target.cellIndexToLabel[piece.targetCell]);
                check(cell.level()==maxLevel,"coarse-neighbor prerequisite left a required target below requested level four");
                const auto h=cell.size();check(h.x==expectedSpacing&&h.y==expectedSpacing,
                    "required ROI target physical spacing differs from requested level four");
                if(!checked[piece.targetCell]){checked[piece.targetCell]=true;++requiredTargets;}
            }
        }
        check(requiredTargets>static_cast<std::size_t>(band),"deep ROI must have real descendant coverage, not a target tag only");
        sourceIdentity.validate(*source);accepted.validateFor(*source);
        const auto clock=accepted.timeLayer();const auto reportClock=decision.report.sourceClock;
        check(publication.identity()==savedId&&clock.step==savedClock.step&&clock.begin==savedClock.begin
            &&clock.end==savedClock.end&&clock.position==savedClock.position
            &&reportClock.step==savedClock.step&&reportClock.begin==savedClock.begin
            &&reportClock.end==savedClock.end&&reportClock.position==savedClock.position
            &&accepted.values().normalVelocity()==savedU,"deep grid planning changed source owner/fullclock/velocity");
        for(std::size_t ci=0;ci<savedStocks.size();++ci){
            check(publication.plicCells()[ci].inventory==savedStocks[ci],"deep grid planning changed independent L/G stocks");
            check(source->cells(source->cellIndexToLabel[ci]).level()==0,"deep grid planning changed source topology");
        }
        totalRequiredTargets+=requiredTargets;
        std::cout<<"HYDRO_VOF_REFINEMENT_PREREQUISITES_SCENE shape="<<(rz?"RZ":"XY")<<" base="<<n
            <<" source="<<source->cellIndexNum()<<" seeds="<<seeds<<" halo="<<band
            <<" target="<<target.cellIndexNum()<<" required_level4="<<requiredTargets
            <<" h="<<expectedSpacing<<" sweeps="<<decision.report.sweeps<<'\n';
    }
    std::cout<<"HYDRO_VOF_REFINEMENT_PREREQUISITES_PASS scenes=2 level=4 requiredTargets="<<totalRequiredTargets<<'\n';
}
}
int main(int argc,char** argv){
    if(argc==2&&std::string(argv[1])=="--regrid-interval-only"){
        try{regridInterval();return 0;}
        catch(const std::exception& e){std::cerr<<"HYDRO_VOF_REGRID_INTERVAL_FAIL "<<e.what()<<'\n';return 1;}
    }
    if(argc==2&&std::string(argv[1])=="--cell-requirements-only"){
        try{cellMinimumRequirements();return 0;}
        catch(const std::exception& e){std::cerr<<"HYDRO_VOF_CELL_REQUIREMENTS_FAIL "<<e.what()<<'\n';return 1;}
    }
    if(argc==2&&std::string(argv[1])=="--refinement-prerequisites-only"){
        try{refinementPrerequisites();return 0;}
        catch(const std::exception& e){std::cerr<<"HYDRO_VOF_REFINEMENT_PREREQUISITES_FAIL "<<e.what()<<'\n';return 1;}
    }
    try{movingInterface();closedVortex(false);closedVortex(true);microInventoryIndicators();phaseSpecificFloors();cellMinimumRequirements();
        std::cout<<"HYDRO_VOF_REFINEMENT_PASS interface=2 vortex=XY,RZ noop=0clone\n";return 0;}
    catch(const std::exception& e){std::cerr<<"HYDRO_VOF_REFINEMENT_FAIL "<<e.what()<<'\n';return 1;}
}
