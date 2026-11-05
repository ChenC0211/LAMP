#include "Hydrodynamics/HydrodynamicsRun.h"
#include "Hydrodynamics/HydrodynamicsTimeStep.h"
#include "Hydrodynamics/IO/VofCheckpoint.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Geometry/FittedDomainGeometry.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include <map>
#include <functional>
#include "Core/Mesh/FaceSegmentGraph.h"
#include <set>
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/NativeTargetLevelField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Core/Base/TimeLayer.h"
#include "Modules/Interface/GeometricVofPublication.h"
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <cfenv>
#include <cerrno>
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace LAMP::Hydrodynamics;
using namespace LAMP::Interface;

namespace {
// The App selects declared physical L/G curves; Core owns every primitive
// distance and the original level-set sign/closure. This query is consumed
// only while constructing the initial mesh, never as transported geometry.
std::function<double(Vector2D)> initialCsgInterfaceDistance(const HydrodynamicsInput& input,QuadtreeMesh& mesh)
{
    if(!input.csgGeometry||!input.csgSampler)
        throw std::invalid_argument("CSG targetCell band requires its immutable source geometry and sampler");
    const auto& definition=*input.csgGeometry;
    const auto liquid=input.negativePhase-1,gas=input.positivePhase-1;
    std::map<int,std::set<int>> linesByLevelSet;std::map<int,int> selectedLines;
    const auto fluidPair=[&](const auto& boundary){return (boundary.phase==liquid&&boundary.opposedPhase==gas)
        ||(boundary.phase==gas&&boundary.opposedPhase==liquid);};
    for(const auto& boundary:definition.boundaries){
        if(boundary.ignored||boundary.kind==GeometricBoundary::domain||!fluidPair(boundary))continue;
        if(boundary.levelSetIndex<0||boundary.levelSetIndex>=static_cast<int>(definition.levelSetName.size())||boundary.lineIds.empty())
            throw std::invalid_argument("CSG targetCell band has no valid original L/G level-set/physical lines");
        for(const int line:boundary.lineIds){
            if(line<0||line>=static_cast<int>(definition.lines.size()))
                throw std::invalid_argument("CSG targetCell band physical line is outside the original definition");
            if(definition.lines[line].type()!=LineType::straight&&definition.lines[line].type()!=LineType::arc)
                throw std::invalid_argument("CSG targetCell band requires original straight/arc physical interface lines");
            const auto inserted=selectedLines.emplace(line,boundary.levelSetIndex);
            if(!inserted.second&&inserted.first->second!=boundary.levelSetIndex)
                throw std::invalid_argument("CSG targetCell band physical line has ambiguous level-set ownership");
            linesByLevelSet[boundary.levelSetIndex].insert(line);
        }
    }
    if(linesByLevelSet.empty())throw std::invalid_argument("CSG targetCell band requires declared physical liquid/gas boundaries");
    for(const auto& boundary:definition.boundaries){
        if(boundary.ignored||boundary.kind==GeometricBoundary::domain||fluidPair(boundary))continue;
        for(const int line:boundary.lineIds)if(selectedLines.count(line))
            throw std::invalid_argument("CSG targetCell band physical line has ambiguous adjacent phases");
    }
    DictionaryGeometry initialDefinition=definition;
    ImmersedGeometry<QuadtreeMesh> seed(initialDefinition,mesh);
    using Query=std::function<ImmersedGeometry<QuadtreeMesh>::InitialCurveSample(Vector2D)>;
    std::vector<Query> queries;queries.reserve(linesByLevelSet.size());
    for(const auto& entry:linesByLevelSet){std::vector<SignInt> lines;lines.reserve(entry.second.size());
        for(const int line:entry.second)lines.push_back({line,false});
        queries.push_back(seed.initialCurveSampler(mesh,lines,entry.first));
    }
    return [queries=std::move(queries)](Vector2D point){
        double distance=std::numeric_limits<double>::infinity();
        for(const auto& query:queries){const auto value=query(point).unsignedDistance;
            if(!std::isfinite(value)||value<0)throw std::invalid_argument("invalid original CSG physical-interface distance");
            distance=std::min(distance,value);
        }
        return distance;
    };
}

// Passive performance publication preserves the numerical caller's complete
// native/raw floating state, including sticky exception flags and errno.
struct CompletedStepPerformanceEnvironment {
    std::fenv_t native{};int error=errno;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr=_mm_getcsr();
#endif
    CompletedStepPerformanceEnvironment(){if(std::fegetenv(&native)!=0){errno=error;throw std::runtime_error("completed-step performance environment capture failed");}}
    ~CompletedStepPerformanceEnvironment() noexcept {
        if(std::fesetenv(&native)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr);
#endif
        errno=error;
    }
};
}

int runHydrodynamics(const HydrodynamicsInput& input,const RunConfig& config,HydrodynamicsTransaction& tx,const HydrodynamicsRunIO& io){
    if(input.interfaceMethod!=HydrodynamicsInput::InterfaceMethod::Vof)
        throw std::invalid_argument("production Hydrodynamics workflow requires Interface.method vof");
    std::cout<<std::setprecision(8)<<'['<<tx.name()<<"] HYDRO / VOF/PLIC / AC"<<tx.bannerSuffix()<<", grid="
        <<input.nx<<'x'<<input.ny<<" cells="<<static_cast<long long>(input.nx)*input.ny
        <<" force=heightFunction";
    std::cout<<" dt="<<input.run.dt<<" endTime="<<input.run.endTime
#ifdef _OPENMP
        <<" threads="<<omp_get_max_threads()
#else
        <<" threads=1"
#endif
        <<std::endl;
    if((input.run.checkpointEverySteps>0||!input.run.restartDirectory.empty())&&!tx.supportsCheckpoint())
        throw std::invalid_argument("accepted checkpoint/restart is unsupported by this transaction");
    VofCheckpointStore checkpoints(input,tx.name());
    const auto restart=checkpoints.read();
    if(config.checkCase){
        std::cout<<tx.inputOkMessage()<<" (no mesh or output)"<<std::endl;return 0;
    }
    std::unique_ptr<QuadtreeMesh> mesh;
    struct Clear{HydrodynamicsTransaction& tx;~Clear(){tx.clear();}} clearBeforeMesh{tx};
    const auto state=[&]() -> const HydroTwoPhaseState& {return *tx.accepted();};
    const char* stage="output startup";std::int64_t attempted=0;
    try{
        io.begin();stage="initialization";
        std::cout<<(restart?"[RESTART] restoring accepted mesh and physical owners":"[INIT] building mesh and initial VOF/PLIC geometry")<<std::endl;
        { // Initialization aliases end before physical time advancement.
        const auto start=HydrodynamicsTransaction::Clock::now();
        DictionaryMesh meshInput;meshInput.minCoord=input.domainMin;meshInput.maxCoord=input.domainMax;
        meshInput.Nx=input.nx;meshInput.Ny=input.ny;meshInput.maxInitialIterations=0;
        // This production workflow uses VOF/PLIC and target-only AMR, not a physical LS field.
        meshInput.levelSetCount=0;
        const int initialLevel=input.initialInterfaceBand?input.initialRefinementLevel:0;
        const int allocationLevel=std::max(initialLevel,input.vofRegrid?input.vofRegrid->maxLevel:0);
        long long maximumLeaves=static_cast<long long>(input.nx)*input.ny;
        const auto allocationLimit=static_cast<long long>(std::numeric_limits<int>::max());
        for(int level=0;level<allocationLevel;++level){
            if(maximumLeaves>allocationLimit/4)throw std::invalid_argument("requested native refinement exceeds mesh allocation range");
            maximumLeaves*=4;
        }
        // Each entity kind uses this label capacity; parent cells/edges also
        // survive refinement. Reserve headroom beyond the maximum leaf count.
        if(allocationLevel>0&&maximumLeaves>(allocationLimit-16384)/2)
            throw std::invalid_argument("requested VOF native storage exceeds mesh allocation range");
        const auto capacity=allocationLevel>0?2*maximumLeaves+16384:maximumLeaves;
        if(capacity>allocationLimit)throw std::invalid_argument("requested native mesh exceeds allocation range");
        meshInput.maxCellNum=std::max(16384,static_cast<int>(capacity));
        if(allocationLevel>0)meshInput.maxLevel_local={allocationLevel};
        mesh=std::make_unique<QuadtreeMesh>(meshInput);
        if(restart){
            stage="accepted checkpoint mesh restore";checkpoints.restoreMesh(*mesh,*restart);
            stage="accepted checkpoint physical restore";tx.restoreCheckpoint(*mesh,restart->transaction);
            const auto clock=state().timeLayer();
            if(clock.step!=restart->clock.step||clock.begin!=restart->clock.begin||clock.end!=restart->clock.end
                ||clock.position!=restart->clock.position)
                throw std::runtime_error("checkpoint physical and controller clocks differ");
            std::cout<<"[RESTART] step="<<clock.step<<" time="<<std::setprecision(17)<<clock.time()
                <<" next_dt="<<restart->nextDt<<" source="<<restart->directory.u8string()<<std::setprecision(8)<<std::endl;
        }else{
        const int globalMinimum=input.vofRegrid?input.vofRegrid->minLevel:0;
        const int liquidMinimum=input.vofRegrid?input.vofRegrid->liquidMinLevel.value_or(globalMinimum):globalMinimum;
        const int gasMinimum=input.vofRegrid?input.vofRegrid->gasMinLevel.value_or(globalMinimum):globalMinimum;
        const bool phaseSpecificFloors=liquidMinimum>globalMinimum||gasMinimum>globalMinimum;
        const bool targetInitialBand=input.initialInterfaceBand&&input.initialBandScale==HydrodynamicsInput::InitialBandScale::TargetCell;
        const bool localInitialBand=input.initialInterfaceBand&&(input.initialBandScale==HydrodynamicsInput::InitialBandScale::CurrentCell||targetInitialBand);
        if(targetInitialBand&&(!input.csgGeometry||liquidMinimum!=globalMinimum||gasMinimum!=globalMinimum))
            throw std::invalid_argument("initial targetCell band requires CSG geometry with the common phase floor");
        if(input.vofRegrid&&input.vofRegrid->minLevel>0){
            // This is a uniform floor, independent of the initial phi query.
            DictionaryMesh::RefineBox floor;floor.level=input.vofRegrid->minLevel;
            floor.min=input.domainMin;floor.max=input.domainMax;mesh->applyRefineBoxes({floor});
        }
        if(input.initialInterfaceBand&&!localInitialBand){
            std::vector<DictionaryMesh::RefineBox> boxes;
            for(int ci=0;ci<mesh->cellIndexNum();++ci){
                const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);const auto center=cell.centerCoord();
                const auto h=std::max(cell.size().x,cell.size().y);const auto j=input.initialJet(center);
                const auto gradient=std::hypot(j.derivative,j.partialSecond);
                bool negative=j.value<0,positive=j.value>0;
                for(auto v:cell.cellToVert){const auto value=input.initialJet(mesh->verts(v).coord()).value;negative|=value<0;positive|=value>0;}
                const bool band=gradient>0&&std::abs(j.value)<=input.initialInterfaceBandCells*h*gradient;
                if(!(band||(negative&&positive)||j.value==0))continue;
                DictionaryMesh::RefineBox box;box.level=input.initialRefinementLevel;
                box.min=center-cell.size()*.5;box.max=center+cell.size()*.5;boxes.push_back(box);
            }
            // Preserve the original base/floor phi-selected boxes. Only the
            // target-only native construction changes; no child phi resampling
            // or Interface state/phase transfer occurs before initialization.
            const int floorLevel=input.vofRegrid?input.vofRegrid->minLevel:0;
            const int maximumSweeps=allocationLevel+2;bool settled=false;
            for(int sweep=0;sweep<=maximumSweeps;++sweep){
                LAMP::Mesh::NativeTargetLevelField targets(*mesh,static_cast<double>(floorLevel)),gradient(*mesh,0.);
                for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto center=mesh->cells(mesh->cellIndexToLabel[ci]).centerCoord();
                    int level=floorLevel;
                    for(const auto& box:boxes)if(center.x>=box.min.x&&center.x<=box.max.x&&center.y>=box.min.y&&center.y<=box.max.y)
                        level=std::max(level,box.level);
                    targets.value[ci]=level;
                }
                // A requested fine leaf may first need its coarser face
                // neighbors refined. Complete this temporary target wave
                // before interpreting a false native preflight as settled.
                mesh->completeRefinementPrerequisites(targets);
                if(!mesh->mayAdaptToTargets(targets)){settled=true;break;}
                if(sweep>=maximumSweeps)throw std::invalid_argument("initial VOF interface-band targets did not settle within allocated-level bound");
                std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
                std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
                std::vector<LAMP::Mesh::NativeTargetLevelField*> cells{&targets};
                mesh->adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&targets);
            }
            if(!settled||mesh->cellIndexNum()>maximumLeaves)
                throw std::invalid_argument("initial VOF interface-band construction exceeds settled native leaf budget");
            for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
                const int level=cell.level();const auto center=cell.centerCoord();int required=floorLevel;
                for(const auto& box:boxes)if(center.x>=box.min.x&&center.x<=box.max.x&&center.y>=box.min.y&&center.y<=box.max.y)
                    required=std::max(required,box.level);
                if(level<required)throw std::invalid_argument("initial VOF interface-band did not reach the frozen requested leaf level");
                if(level<floorLevel||level>allocationLevel)throw std::invalid_argument("initial VOF interface-band leaf outside allocated floor/level");}
        }
        if(localInitialBand||phaseSpecificFloors){
            std::function<double(Vector2D)> csgDistance;double targetBandWidth=0;
            std::optional<FittedDomainGeometry::FixedSolidSelection> initialSolids;
            if(targetInitialBand){
                csgDistance=initialCsgInterfaceDistance(input,*mesh);
                const LAMP::Mesh::DyadicGridCoordinates coordinates(input.domainMin,input.domainMax,input.nx,input.ny);
                const auto h=coordinates.cellSize(static_cast<unsigned>(input.initialRefinementLevel));
                targetBandWidth=input.initialInterfaceBandCells*std::max(h.x,h.y);
                if(!(std::isfinite(targetBandWidth)&&targetBandWidth>0))
                    throw std::invalid_argument("initial targetCell band has an unrepresentable physical width");
            }
            const auto refreshInitialSolids=[&]{if(targetInitialBand)
                initialSolids=FittedDomainGeometry::selectFixedSolids(*mesh,*input.csgGeometry,*input.csgSampler,input.csgExcludedPhases);};
            // Geometry queries are allowed only here, before the one inventory
            // initialization. Re-evaluate child scale after every native wave.
            const auto requiredLevel=[&](const auto& cell){
                const auto center=cell.centerCoord();
                if(targetInitialBand){
                    if(!initialSolids)throw std::logic_error("initial CSG band lacks current-mesh solid certification");
                    if(initialSolids->excludedPhaseByCell.at(static_cast<std::size_t>(cell.index()))!=0)return globalMinimum;
                    double radius=0;for(const int v:cell.cellToVert){const auto point=mesh->verts(v).coord();
                        radius=std::max(radius,std::hypot(point.x-center.x,point.y-center.y));}
                    const auto distance=csgDistance(center);
                    return static_cast<long double>(distance)<=static_cast<long double>(targetBandWidth)+radius
                        ?std::max(globalMinimum,input.initialRefinementLevel):globalMinimum;
                }
                const auto j=input.initialJet(center);
                bool negative=j.value<0,positive=j.value>0;
                for(auto v:cell.cellToVert){const auto value=input.initialJet(mesh->verts(v).coord()).value;
                    negative|=value<0;positive|=value>0;}
                int required=globalMinimum;
                if(negative)required=std::max(required,liquidMinimum);
                if(positive)required=std::max(required,gasMinimum);
                if(!negative&&!positive)required=std::max(liquidMinimum,gasMinimum);
                const auto h=std::max(cell.size().x,cell.size().y);
                const auto gradient=std::hypot(j.derivative,j.partialSecond);
                const bool band=gradient>0&&std::abs(j.value)<=input.initialInterfaceBandCells*h*gradient;
                if(localInitialBand&&(band||(negative&&positive)||j.value==0))
                    required=std::max(required,input.initialRefinementLevel);
                return required;
            };
            const int maximumSweeps=allocationLevel+2;bool settled=false;
            for(int sweep=0;sweep<=maximumSweeps;++sweep){
                refreshInitialSolids();
                LAMP::Mesh::NativeTargetLevelField targets(*mesh,0.),gradient(*mesh,0.);
                for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
                    // adapt visits levels in one call. A target of L5 here
                    // would propagate a coarse box to L5 before child queries.
                    targets.value[ci]=std::max(cell.level(),std::min(cell.level()+1,requiredLevel(cell)));
                }
                mesh->completeRefinementPrerequisites(targets);
                if(!mesh->mayAdaptToTargets(targets)){settled=true;break;}
                if(sweep>=maximumSweeps)throw std::invalid_argument("initial VOF local-band/phase-floor targets did not settle");
                std::vector<VertField<double>> levelSets;std::vector<VertField<double>*> vertices;
                std::vector<HEdgeField<double>*> horizontal;std::vector<VEdgeField<double>*> vertical;
                std::vector<LAMP::Mesh::NativeTargetLevelField*> cells{&targets};
                mesh->adapt(levelSets,gradient,vertices,horizontal,vertical,cells,&targets);
            }
            if(!settled||mesh->cellIndexNum()>maximumLeaves)
                throw std::invalid_argument("initial VOF local-band/phase-floor exceeds settled native leaf budget");
            refreshInitialSolids();
            for(int ci=0;ci<mesh->cellIndexNum();++ci){const auto& cell=mesh->cells(mesh->cellIndexToLabel[ci]);
                if(cell.level()<requiredLevel(cell)||cell.level()>allocationLevel)
                    throw std::invalid_argument("initial VOF local-band/phase-floor did not reach its required level");}
        }
        const auto vofPublication=[&]{
            if(!input.csgGeometry)return GeometricVofInitialization::prepare(*mesh,input.radial,[&](Vector2D p){
                return input.initialJet(p).value;
            },input.vofInitialization);
            if((input.initialInterfaceBand&&!targetInitialBand)||input.run.checkpointEverySteps>0||!input.run.restartDirectory.empty())
                throw std::invalid_argument("CSG initialization requires uniform or explicit targetCell initial band without restart/checkpoint");
            const auto& definition=*input.csgGeometry;
            if(!input.csgSampler)throw std::invalid_argument("CSG initial source has no validated owned primitive sampler");
            const auto& sampler=*input.csgSampler;
            const auto solids=FittedDomainGeometry::selectFixedSolids(*mesh,definition,sampler,input.csgExcludedPhases);
            FaceSegmentGraph selected;selected.build(*mesh,[&](const QuadtreeMesh&,const Cell& cell){
                return solids.excludedPhaseByCell.at(cell.index())==0;});
            const auto domain=NativeLeafDomain::build(*mesh,selected,input.radial);
            auto publication=GeometricVofInitialization::prepare(*mesh,domain,[&](Vector2D p){
                return sampler.phaseValue(input.negativePhase,p);
            },input.vofInitialization);
            if(!input.csgPin)return publication;
            const auto& pin=*input.csgPin;const auto point=definition.points.at(static_cast<std::size_t>(pin.pointId-1));
            FaceSegmentGraph complete;complete.build(*mesh);
            if(solids.boundaryIdsBySegment.size()!=complete.segments().size())throw std::logic_error("different actual CSG boundary support");
            std::set<int> wet,dry;
            for(std::size_t i=0;i<complete.segments().size();++i){const int id=solids.boundaryIdsBySegment[i];
                if(id!=pin.wetBoundaryId&&id!=pin.dryBoundaryId)continue;
                const auto& segment=complete.segments()[i];const auto& edge=segment.isHorizontal?mesh->hEdges(segment.edgeLabel):mesh->vEdges(segment.edgeLabel);
                bool endpoint=false;for(int v:edge.edgeToVert){const auto p=mesh->verts(v).coord();endpoint=endpoint||(p.x==point.x&&p.y==point.y);}
                if(!endpoint)continue;
                const int row=edge.index()+(segment.isHorizontal?mesh->vEdgeIndexNum():0);
                (id==pin.wetBoundaryId?wet:dry).insert(row);
            }
            if(wet.size()!=1||dry.size()!=1)throw std::invalid_argument("CSG Pin requires one actual wet/dry boundary-ID wall leg meeting the exact declared point");
            return GeometricVofPublication::withPinnedCorner(*mesh,publication,point,*wet.begin(),*dry.begin(),pin.wallTrace);
        }();
        const auto domain=vofPublication.domain();
        if(input.initialInterfaceBand||phaseSpecificFloors)for(const auto face:domain.faces(*mesh)){
            if(face.negativeCell<0||face.positiveCell<0)continue;
            const auto& a=mesh->cells(mesh->cellIndexToLabel[face.negativeCell]);
            const auto& b=mesh->cells(mesh->cellIndexToLabel[face.positiveCell]);
            if(std::abs(a.level()-b.level())>1)
                throw std::invalid_argument("initial VOF interface-band violates actual native-face 2:1 levels");
        }
        const auto layout=HydroDofLayout::build(*mesh,domain);
        const auto metrics=HydroPhysicalDivergence::build(*mesh,HydroIntegralDivergence::build(*mesh,layout));
        std::vector<double> velocity(layout.velocityCount()),pressure(layout.pressureCount());
        for(int u=0;u<layout.velocityCount();++u){
            const auto f=domain.face(*mesh,layout.nativeFace(u));
            velocity[u]=f.negativeCell<0||f.positiveCell<0?0:(layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?input.initialVelocity.x:input.initialVelocity.y);
        }
        tx.initializeVof(*mesh,tx.makeInitialVofState(*mesh,
            HydroPhaseGeometry::build(*mesh,metrics,vofPublication),HydroNativeValues(layout,velocity,pressure)),start,stage);
        }
        }
        double nextDt=restart?restart->nextDt:input.run.dt;int withoutRetry=restart?restart->withoutRetry:0;
        if(input.run.checkpointEverySteps>0){
            stage="accepted checkpoint publication";CompletedStepPerformanceEnvironment environment;
            checkpoints.write(*mesh,{state().timeLayer(),nextDt,withoutRetry,tx.checkpoint(*mesh),{}});
        }
        while(state().timeLayer().time()<input.run.endTime){
            // Include active rate queries, all failed attempts and real AMR
            // planning/transfer/commit costs in this accepted step's compute.
            const auto stepStart=HydrodynamicsTransaction::Clock::now();
            const auto begin=state().timeLayer().time();attempted=state().timeLayer().step+1;
            if(attempted>1000000)throw std::runtime_error("accepted step budget exceeded");
            std::optional<HydrodynamicsTimeStepReport> timeStepReport;
            double end;
            if(input.run.adaptive){
                stage="VOF adaptive timestep";
                const auto rates=state().vofPrescribedFlow()
                    ?HydroVofStep::stabilityRatesPrescribedFlow(*mesh,state(),input.vof)
                    :HydroVofStep::stabilityRates(*mesh,state(),input.vof);
                timeStepReport=selectHydrodynamicsTimeStep(input.run,begin,nextDt,
                    static_cast<int>(state().timeLayer().step),rates.flowRate,rates.capillaryAvailable,
                    rates.capillaryRate,rates.interfaceMinimumSpacing,input.vof.maximumCapillaryCourant);
                end=timeStepReport->end;
            }else end=input.run.stepEnd(begin,nextDt,static_cast<int>(state().timeLayer().step));
            if(input.run.detailedLog)std::cout<<"[ADVANCE] step="<<attempted<<" accepted_t="<<begin<<" proposed_dt="<<end-begin<<std::endl;
            int retries=0;
            for(;;){
                stage=tx.prepareStage();
                try{
                    const LAMP::TimeLayer target{attempted,begin,end,LAMP::TimeLayer::Position::End};
                    tx.prepare(*mesh,target,input.step,input.run.detailedLog,stage);
                    break;
                }
                catch(const HydroTwoPhaseStepRetry& e){
                    if(!input.run.stepRetry)throw;
                    const auto& policy=*input.run.stepRetry;const auto oldDt=end-begin,newDt=oldDt*.5;
                    if(retries>=policy.maxRetries||newDt<policy.minDt||!(begin+newDt>begin))throw;
                    ++retries;io.retry(attempted,state().timeLayer().step,begin,oldDt,newDt,e.what());end=begin+newDt;
                }
            }
            stage="accepted commit";
            tx.commit(*mesh);
            const auto compute=std::chrono::duration<double>(HydrodynamicsTransaction::Clock::now()-stepStart).count();
            stage="accepted output";
            const bool final=state().timeLayer().time()==input.run.endTime;
            if(timeStepReport&&io.timeStep)io.timeStep(*timeStepReport,state().timeLayer(),retries);
            tx.output(*mesh,stepStart,compute,final);
            const bool save=input.run.checkpointEverySteps>0
                &&(final||state().timeLayer().step%input.run.checkpointEverySteps==0);
            std::string checkpointPayload;
            if(save){
                stage="accepted checkpoint capture";CompletedStepPerformanceEnvironment environment;
                checkpointPayload=tx.checkpoint(*mesh);
            }
            tx.releaseStep();
            nextDt=state().timeLayer().end-state().timeLayer().begin;if(retries)withoutRetry=0;else ++withoutRetry;
            if(input.run.stepRetry&&input.run.stepRetry->recoverAfter>0&&withoutRetry>=input.run.stepRetry->recoverAfter){
                nextDt=std::min(input.run.dt,nextDt*1.25);withoutRetry=0;
            }
            if(io.stepCompleted){
                stage="accepted full-step performance output";
                CompletedStepPerformanceEnvironment environment;
                const auto completed=HydrodynamicsTransaction::Clock::now();
                io.stepCompleted(state().timeLayer(),std::chrono::duration<double>(completed-stepStart).count());
            }
            if(save){
                stage="accepted checkpoint publication";CompletedStepPerformanceEnvironment environment;
                checkpoints.write(*mesh,{state().timeLayer(),nextDt,withoutRetry,std::move(checkpointPayload),{}});
            }
        }
        stage="completion output";io.complete();
        std::cout<<"[COMPLETED] steps="<<state().timeLayer().step<<" time="<<state().timeLayer().time()
            <<tx.completionSuffix()<<std::endl;return 0;
    }catch(const std::exception& e){
        const HydroTwoPhaseState* accepted=tx.accepted();
        const auto step=accepted?accepted->timeLayer().step:0;const double time=accepted?accepted->timeLayer().time():0;
        try{io.failure(attempted,step,time,stage,e.what());}catch(const std::exception& io){std::cerr<<"[FAILURE IO] "<<io.what()<<std::endl;}
        std::cerr<<"[FAILED] stage="<<stage<<" attempted_step="<<attempted<<" accepted_step="<<step<<" accepted_t="<<time<<" reason="<<e.what()<<std::endl;return 1;
    }
}


