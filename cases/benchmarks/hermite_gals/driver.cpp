#include "interface_first_jet_fixture.h"
namespace {
std::optional<Errors> runFirstJetCase(Case s,int n,double dt,int steps,bool dynamicMesh=false)
{
    check(s.pureFirstJet,"explicit FirstJet case required");
    const auto label=s.name+"_"+(s.radial?"RZ":"XY")+"_"+(s.amr?"AMR":"uniform");
    std::optional<HermiteGalsPublication> accepted;std::string stage="initialize";
    try{
        auto mesh=makeMesh(n,s.amr?1:0);accepted=initial(*mesh,s);const auto initialState=accepted->state();bool currentAmr=s.amr;
        std::vector<Vector2D> initialPoints;if(dynamicMesh)for(const auto& v:initialState.primary().vertices)initialPoints.push_back(point(*mesh,v.vertex));
        const auto initialVolume=accepted->volumes().report().totalRawLiquid;
        check(initialVolume>0,"pure initial volume must be positive");
        if(s.amr){std::size_t coarse=0,fine=0;for(const auto& c:accepted->volumes().cells())if(c.certifiedPhase==0){
            if(mesh->cells(mesh->cellIndexToLabel[c.cell]).level()==0)++coarse;else ++fine;}
            check(coarse && fine,"pure AMR interface must cut both levels");
            std::cout<<"FIRST_JET_AMR_CUT,"<<label<<','<<coarse<<','<<fine<<'\n';}
        long double maximumDrift=0;
        for(int k=0;k<steps;++k){const auto& old=accepted->state();const auto id=old.identity();
            const double begin=old.timeLayer().time();
            try{
                check(old.primary().transport==HermiteGalsTransport::FirstJet && old.primary().phaseVolumes.empty(),"pure advance received inventory");
                auto next=advanceFirstJet(*mesh,*accepted,s,{old.timeLayer().step+1,begin,begin+dt,Position::End},stage);
                accepted=std::move(next);
            }catch(...){check(accepted->state().identity()==id,"pure rejection changed accepted identity");throw;}
            if(dynamicMesh && (k+1==steps/4 || k+1==steps/2 || k+1==3*steps/4)){
                stage="first_jet_mesh_transfer";auto target=makeMesh(n,currentAmr?0:1);const auto sourceId=accepted->state().identity();
                const auto before=errors(*accepted,s);const auto time=accepted->state().timeLayer();
                const auto kind=currentAmr?HermiteGalsMeshTransferKind::FirstJetProjection:HermiteGalsMeshTransferKind::ExactRefinement;
                const auto transfer=transferHermiteGalsFirstJetMesh(*mesh,*accepted,*target,kind);const auto& next=transfer.publication();
                const auto after=errors(next,s);const auto t=next.state().timeLayer();
                check(accepted->state().identity()==sourceId && next.state().identity()!=sourceId && t.step==time.step
                    && t.begin==time.begin && t.end==time.end && t.position==time.position,"dynamic regrid changed source/time");
                const auto& r=transfer.report();currentAmr=!currentAmr;
                std::cout<<"FIRST_JET_REGRID,"<<label<<','<<n<<','<<k+1<<','<<time.time()<<','<<(kind==HermiteGalsMeshTransferKind::ExactRefinement?"exact":"project")
                    <<','<<currentAmr<<','<<r.refinedCells<<','<<r.coarsenedCells<<','<<r.unchangedCells<<','<<r.relativeDomainVolumeChange
                    <<','<<before.curvature<<','<<after.curvature<<'\n';
                accepted=next;mesh=std::move(target);
            }
            const auto volume=accepted->volumes().report().totalRawLiquid,drift=(volume-initialVolume)/initialVolume;
            maximumDrift=std::max(maximumDrift,std::abs(drift));
            std::cout<<"FIRST_JET_ACCEPTED,"<<label<<','<<n<<','<<k+1<<','<<accepted->state().timeLayer().time()
                <<','<<initialVolume<<','<<volume<<','<<drift<<",inventory=NA,flux_ledger=NA\n";
        }
        stage="error_reference";auto e=errors(*accepted,s);
        if(s.motion==0){const auto& a=initialState.primary();const auto& b=accepted->state().primary();
            if(dynamicMesh){for(std::size_t k=0;k<a.vertices.size();++k){const auto j=accepted->state().geometry().snapshot().scalarJetNative(initialPoints[k]);
                e.drift=std::max({e.drift,std::abs(a.vertices[k].phi-j.value),std::abs(a.vertices[k].dx-j.dx),std::abs(a.vertices[k].dy-j.dy)});}}
            else for(std::size_t k=0;k<a.vertices.size();++k)e.drift=std::max({e.drift,std::abs(a.vertices[k].phi-b.vertices[k].phi),
                    std::abs(a.vertices[k].dx-b.vertices[k].dx),std::abs(a.vertices[k].dy-b.vertices[k].dy)});
            check(e.drift<=1e-11,"pure stationary drift gate");}
        if(s.motion<=1 && !s.nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10 && e.cellVolume<=5e-12,
            "pure Q2 static/translation accuracy gate");
        std::cout<<"FIRST_JET_BENCHMARK,"<<label<<','<<n<<','<<dt<<','<<steps<<','<<accepted->state().timeLayer().time()
            <<','<<e.position<<','<<e.gradient<<','<<e.curvature<<','<<e.volume<<','<<e.cellVolume<<','<<maximumDrift<<",PASS\n";
        return e;
    }catch(const std::exception& ex){++failures;std::cout<<"FIRST_JET_FAILURE,"<<label<<','<<n<<','<<dt<<','<<steps
        <<','<<(accepted?accepted->state().timeLayer().time():0)<<','<<stage<<','<<ex.what()<<'\n';}
    return std::nullopt;
}
void firstJetBenchmark(bool smoke,bool dynamicMesh=false)
{
    int completed=0;
    for(bool radial:{false,true})for(bool amr:{false,true}){
        Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;if(radial)s.radius=.27;
        s.name=dynamicMesh?"regrid_static":"pure_static";if(runFirstJetCase(s,8,.005,4,dynamicMesh))++completed;
        s.name=dynamicMesh?"regrid_translation":amr?"pure_amr_crossing":"pure_translation";s.motion=1;
        if(runFirstJetCase(s,8,.01,smoke?(dynamicMesh?4:1):40,dynamicMesh))++completed;
    }
    for(int motion:{2,3})for(bool amr:{false,true}){
        Case s;s.pureFirstJet=true;s.name=dynamicMesh?(motion==2?"regrid_rotation":"regrid_vortex"):(motion==2?"pure_rotation":"pure_vortex");
        s.cx=.65;s.cy=.53;s.radius=.13;s.motion=motion;s.speed=motion==2?2:40;s.period=.4;s.amr=amr;
        if(runFirstJetCase(s,8,motion==2?.01:.001,smoke?(dynamicMesh?4:1):motion==2?30:400,dynamicMesh))++completed;
    }
    std::cout<<(dynamicMesh?"I12p2b FirstJet regrid ":"I12p1 FirstJet ")<<(smoke?"smoke":"full")<<' '<<(failures?"FAIL":"PASS")
        <<" completed="<<completed<<" planned=12 failures="<<failures<<'\n';
}
void firstJetCoarseningStudy()
{
    int completed=0;
    for(bool radial:{false,true})for(bool sourceAmr:{false,true}){
        Case s;s.pureFirstJet=true;s.radial=radial;s.nonpolynomial=true;s.motion=1;if(radial)s.radius=.27;
        std::optional<Errors> previous;bool complete=true,decreasing=true;double lastG=0,lastK=0;
        for(int n:{8,16,32}){std::string stage="initialize";
            try{
                auto fine=makeMesh(sourceAmr?n:2*n,sourceAmr?1:0),coarse=makeMesh(n);
                auto accepted=initial(*fine,s);
                const auto advanceOn=[&](const QuadtreeMesh& m){const auto& old=accepted.state();const auto t=old.timeLayer();
                    accepted=advanceFirstJet(m,accepted,s,{t.step+1,t.time(),t.time()+.01,Position::End},stage);};
                for(int k=0;k<2;++k)advanceOn(*fine);
                stage="advected_source_projection";const auto sourceId=accepted.state().identity();
                const auto transfer=transferHermiteGalsFirstJetMesh(*fine,accepted,*coarse,HermiteGalsMeshTransferKind::FirstJetProjection);
                check(accepted.state().identity()==sourceId && transfer.report().coarsenedCells>0
                    && transfer.publication().state().timeLayer().time()==.02,"coarsening study did not transfer advected source");
                const auto volumeChange=transfer.report().relativeDomainVolumeChange;
                accepted=transfer.publication();for(int k=0;k<2;++k)advanceOn(*coarse);
                check(accepted.state().timeLayer().time()==.04 && accepted.state().primary().phaseVolumes.empty(),"coarsening continuation lost time/mode");
                const auto e=errors(accepted,s);
                const double gp=previous && previous->gradient>0 && e.gradient>0?std::log2(previous->gradient/e.gradient):0,
                    kp=previous && previous->curvature>0 && e.curvature>0?std::log2(previous->curvature/e.curvature):0;
                if(previous){decreasing=decreasing && e.gradient<previous->gradient && e.curvature<previous->curvature;lastG=gp;lastK=kp;}
                std::cout<<"FIRST_JET_COARSENING_SPACE,"<<radial<<','<<sourceAmr<<','<<n<<','<<e.position<<','<<e.gradient<<','<<e.curvature
                    <<','<<e.volume<<','<<e.cellVolume<<','<<volumeChange<<','<<gp<<','<<kp<<'\n';
                previous=e;++completed;
            }catch(const std::exception& ex){complete=false;++failures;previous.reset();
                std::cout<<"FIRST_JET_COARSENING_FAILURE,"<<radial<<','<<sourceAmr<<','<<n<<','<<stage<<','<<ex.what()<<'\n';}
        }
        if(!complete || !decreasing || lastG<2 || lastK<1.5){++failures;
            std::cout<<"FIRST_JET_COARSENING_ORDER_FAILURE,"<<radial<<','<<sourceAmr<<','<<lastG<<','<<lastK<<'\n';}
    }
    std::cout<<"I12p2b FirstJet coarsening "<<(failures?"FAIL":"PASS")<<" completed="<<completed<<" planned=12 failures="<<failures<<'\n';
}
void firstJetSpaceStudy(bool stationary)
{
    int completed=0;
    for(bool radial:{false,true})for(bool amr:{false,true}){
        Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=true;s.motion=stationary?0:1;
        s.name=stationary?"pure_curvature":"pure_nonpolynomial_translation";if(radial)s.radius=.27;
        std::optional<Errors> previous;double lastG=0,lastK=0;bool decreasing=true,complete=true;
        for(int n:{8,16,32}){
            std::optional<Errors> current;
            if(stationary){try{auto mesh=makeMesh(n,amr?1:0);current=errors(initial(*mesh,s),s);}
                catch(const std::exception& e){++failures;std::cout<<"FIRST_JET_SPACE_FAILURE,"<<stationary<<','<<radial<<','<<amr<<','<<n<<','<<e.what()<<'\n';}}
            else current=runFirstJetCase(s,n,.01,4);
            if(current){++completed;const auto& e=*current;
                const double gp=previous && previous->gradient>0 && e.gradient>0?std::log2(previous->gradient/e.gradient):0,
                    kp=previous && previous->curvature>0 && e.curvature>0?std::log2(previous->curvature/e.curvature):0;
                if(previous){decreasing=decreasing && e.gradient<previous->gradient && e.curvature<previous->curvature;lastG=gp;lastK=kp;}
                std::cout<<"FIRST_JET_SPACE,"<<stationary<<','<<radial<<','<<amr<<','<<n<<','<<e.position<<','<<e.gradient<<','<<e.curvature
                    <<','<<e.volume<<','<<e.cellVolume<<','<<gp<<','<<kp<<'\n';
            }else complete=false;
            previous=current;
        }
        if(!complete || !decreasing || lastG<2 || lastK<1.5){++failures;
            std::cout<<"FIRST_JET_SPACE_ORDER_FAILURE,"<<stationary<<','<<radial<<','<<amr<<','<<lastG<<','<<lastK<<'\n';}
    }
    std::cout<<"I12p2a FirstJet "<<(stationary?"curvature":"translation space")<<' '<<(failures?"FAIL":"PASS")
        <<" completed="<<completed<<" planned=12 failures="<<failures<<'\n';
}
void firstJetTimeStudy()
{
    int completed=0;
    for(bool amr:{false,true}){Case s;s.pureFirstJet=true;s.amr=amr;s.name="pure_rotation_time";
        s.cx=.65;s.cy=.53;s.radius=.13;s.motion=2;s.speed=2;
        std::optional<Errors> previous;bool complete=true,decreasing=true;double lastOrder=0;
        for(int refinement:{1,2,4}){
            const double dt=.01/refinement;const int steps=30*refinement;const auto current=runFirstJetCase(s,8,dt,steps);
            if(current){++completed;const double order=previous && previous->position>0 && current->position>0
                    ?std::log2(previous->position/current->position):0;
                if(previous){decreasing=decreasing && current->position<previous->position;lastOrder=order;}
                std::cout<<"FIRST_JET_TIME,"<<amr<<','<<dt<<','<<steps<<','<<current->position<<','<<current->gradient<<','<<current->curvature
                    <<','<<current->volume<<','<<order<<'\n';
            }else complete=false;
            previous=current;
        }
        if(!complete || !decreasing || lastOrder<1.8){++failures;std::cout<<"FIRST_JET_TIME_ORDER_FAILURE,"<<amr<<','<<lastOrder<<'\n';}
    }
    std::cout<<"I12p2a FirstJet time "<<(failures?"FAIL":"PASS")<<" completed="<<completed<<" planned=6 failures="<<failures<<'\n';
}
}
int main(int argc,char** argv) try {
    std::cout<<std::scientific<<std::setprecision(16);
    if(argc!=2) throw std::invalid_argument("full | regrid | space | curvature | time | coarsening | smoke");
    const std::string mode=argv[1];
    if(mode=="full") firstJetBenchmark(false);
    else if(mode=="regrid") firstJetBenchmark(false,true);
    else if(mode=="smoke") firstJetBenchmark(true,true);
    else if(mode=="space") firstJetSpaceStudy(false);
    else if(mode=="curvature") firstJetSpaceStudy(true);
    else if(mode=="time") firstJetTimeStudy();
    else if(mode=="coarsening") firstJetCoarseningStudy();
    else throw std::invalid_argument("unknown study");
    return failures?1:0;
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}