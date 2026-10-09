#include "interface_first_jet_fixture.h"
namespace {
void firstJetPublicationContract()
{
    int cases=0,guards=0;
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool nonpolynomial:{false,true}){
        auto mesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=nonpolynomial;
        if(radial)s.radius=.27;
        const auto initialGeometry=initial(*mesh,s);auto p=initialGeometry.state().primary();
        check(p.transport==HermiteGalsTransport::FirstJet && p.phaseVolumes.empty(),"pure initial state contains M");
        HermiteGalsCandidate candidate(*mesh,{},p);const auto id=candidate.identity();
        const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate,s.symmetry(*mesh));
        HermiteGalsCenterPrediction predicted{id,centers(*mesh,s.profile(0))};
        const auto reject=[&](auto f,const std::string& message){rejects(f,message);++guards;};
        auto bad=p;bad.phaseVolumes.push_back({p.geometryCells.front(),0});
        reject([&]{HermiteGalsCandidate invalid(*mesh,{},bad);},"FirstJet state cannot carry conserved inventory");
        bad=p;bad.transport=static_cast<HermiteGalsTransport>(-1);
        reject([&]{HermiteGalsCandidate invalid(*mesh,{},bad);},"unknown transport mode");
        bad=p;bad.transport=HermiteGalsTransport::ConservativeVolume;
        reject([&]{HermiteGalsCandidate invalid(*mesh,{},bad);},"geometry cell lacks conserved inventory");
        HermiteGalsCandidate other(*mesh,{},p);const auto foreign=closeHermiteGalsMixedAmr(*mesh,other,s.symmetry(*mesh));
        auto wrong=predicted;wrong.primaryIdentity=other.identity();
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,wrong);},"FirstJet primary revision mismatch");
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),foreign,predicted);},"FirstJet primary revision mismatch");
        wrong=predicted;wrong.values.pop_back();
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,wrong);},"FirstJet closure layout mismatch");
        wrong=predicted;wrong.values.front()=std::numeric_limits<double>::quiet_NaN();
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,wrong);},"FirstJet nonfinite center prediction");
        HermiteGalsPublicationConfig invalidConfig;invalidConfig.unitNormalTolerance=0;
        reject([&]{publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,predicted,invalidConfig);},"invalid publication tolerance");
        reject([&]{closeHermiteGalsBubble(*mesh,candidate,mixed,{id,predicted.values});},"inventory requires complete scalar support");
        check(candidate.identity()==id && !candidate.hasGeometry() && candidate.primary().phaseVolumes.empty(),"pure publication failure mutated candidate");
        const auto pub=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,predicted);
        reject([&]{candidate.identity();},"moved-from or consumed handle");
        reject([&]{pub.volumes().initializeCandidate(*mesh);},"FirstJet state cannot carry conserved inventory");
        check(pub.state().identity()==id && pub.state().timeLayer().time()==0 && pub.state().primary().phaseVolumes.empty(),"pure publication authority");
        check(pub.state().closure().cellBubbleCoefficients==predicted.values
            && pub.state().closure().vertexMixedDerivatives==mixed.vertexMixedDerivatives(),"pure publication adjusted closure");
        const auto& scalar=pub.state().geometry().snapshot();
        check(&pub.volumes().geometry().snapshot()==&scalar && &pub.facePartitions().source()==&scalar,"pure geometry has multiple owners");
        check(std::isnan(pub.report().maxRelativeInventoryResidual) && std::isnan(pub.volumes().report().totalSourceInventory)
            && std::isnan(pub.volumes().report().geometricInventoryDifference),"pure mode claims inventory residual");
        check(pub.volumes().cells().size()==p.geometryCells.size() && pub.report().quadratureNodes>0,"pure geometric diagnostics missing");
        for(std::size_t k=0;k<p.geometryCells.size();++k){const int cell=p.geometryCells[k];
            close(scalar.scalarJetInNativeCell(cell,mesh->cells(mesh->cellIndexToLabel[cell]).centerCoord()).value,
                predicted.values[k],2e-15,"pure center interpolation");}
        for(std::size_t k=0;k<p.vertices.size();++k){const auto& a=p.vertices[k];const auto& b=pub.state().primary().vertices[k];
            check(a.phi==b.phi && a.dx==b.dx && a.dy==b.dy,"pure publication changed master first jet");}
        const auto e=errors(pub,s);if(!nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10,"pure Q2 geometry reproduction");
        std::cout<<"FIRST_JET_PUBLICATION,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<e.position<<','<<e.gradient<<','<<e.curvature<<','<<e.volume<<",PASS\n";
        ++cases;
    }
    std::cout<<"I12p1 FirstJet publication PASS cases="<<cases<<" guards="<<guards<<'\n';
}
void firstJetStepContract()
{
    int cases=0,guards=0;
    const auto reject=[&](auto f,const std::string& why){rejects(f,why);++guards;};
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool nonpolynomial:{false,true}){
        auto mesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=nonpolynomial;
        s.motion=1;if(radial)s.radius=.27;
        const auto source=initial(*mesh,s);const auto id=source.state().identity();const auto& scalar=source.state().geometry().snapshot();
        const auto domain=traceDomain(*mesh,radial);HermiteGalsPublicationConfig config;config.volume.symmetryPlanes=s.symmetry(*mesh);
        const HermiteGalsTimeLayer time{1,0,.01,Position::End};HermiteGalsCandidate authority(*mesh,time,source.state().primary());
        const auto velocity=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),time,s.velocity({},0).value);
        std::size_t queried=0;const HermiteGalsBoundaryFirstJetQuery boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){
            check(q.x==mesh->minCoord.x || q.x==mesh->maxCoord.x || q.y==mesh->minCoord.y || q.y==mesh->maxCoord.y,"boundary query used for interior reset");
            check(t.step==1 && t.begin==0 && t.end==.01,"boundary query lost target clock");++queried;
            const auto j=s.exact(q,t.time());return HermiteGalsBoundaryFirstJet{j.v,{j.x,j.y}};};
        for(Position position:{Position::Midpoint,Position::End}){
            const HermiteGalsTimeLayer target{1,0,.01,position};const auto arrival=velocity.constantArrivalView(authority.identity(),target);
            queried=0;auto result=prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),arrival,boundary,config);
            const auto& pub=result.publication();const auto& r=result.report();const auto& p=pub.state().primary();
            check(result.sourceIdentity()==id && result.arrivalVelocityIdentity()==authority.identity() && pub.state().identity()!=id
                && pub.state().identity()!=authority.identity() && pub.state().timeLayer().time()==target.time(),"step provenance/time authority");
            check(p.transport==HermiteGalsTransport::FirstJet && p.phaseVolumes.empty() && r.tracedCenters==p.geometryCells.size(),"step introduced inventory or lost centers");
            check(r.transportedVertices+r.prescribedVertices==p.vertices.size() && r.transportedEdges+r.prescribedEdges==p.edges.size()
                && r.prescribedVertices && r.prescribedEdges && queried==r.prescribedVertices+r.prescribedEdges,"step boundary/transport coverage");
            check(&pub.volumes().geometry().snapshot()==&pub.state().geometry().snapshot()
                && &pub.facePartitions().source()==&pub.state().geometry().snapshot(),"step has multiple scalar owners");
            const auto e=errors(pub,s);if(!nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10 && e.cellVolume<=5e-12,"step Q2 accuracy gate");
            std::cout<<"FIRST_JET_STEP,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<(position==Position::Midpoint?"midpoint":"end")
                <<','<<r.transportedVertices<<','<<r.transportedEdges<<','<<r.tracedCenters<<','<<r.prescribedVertices<<','<<r.prescribedEdges<<",PASS\n";
            if(position==Position::End){const auto retained=std::move(result);reject([&]{result.publication();},"moved-from result");
                check(retained.sourceIdentity()==id,"moved step lost provenance");}
            ++cases;
        }
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),{},config);},"requires prescribed first jet");
        const HermiteGalsBoundaryFirstJetQuery nan=[](Vector2D,HermiteGalsTimeLayer,Boundary){return HermiteGalsBoundaryFirstJet{std::numeric_limits<double>::quiet_NaN(),{0,0}};};
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),nan,config);},"nonfinite prescribed boundary");
        const HermiteGalsBoundaryFirstJetQuery throwing=[](Vector2D,HermiteGalsTimeLayer,Boundary)->HermiteGalsBoundaryFirstJet{throw std::runtime_error("fixture boundary failure");};
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),throwing,config);},"fixture boundary failure");
        const HermiteGalsVelocityView wrongPrevious(authority.identity(),source.state().timeLayer(),[&](Vector2D q){return s.velocity(q,0);});
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,wrongPrevious,velocity.arrival(),boundary,config);},"previous velocity geometry differs");
        auto wrongMesh=makeMesh(8,amr?0:1);
        reject([&]{prepareHermiteGalsFirstJetStep(*wrongMesh,source,domain,velocity.previous(),velocity.arrival(),boundary,config);},"mesh");
        const HermiteGalsVelocityView wrongStep(authority.identity(),{2,0,.01,Position::End},[&](Vector2D q){return s.velocity(q,.01);});
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),wrongStep,boundary,config);},"next step");
        const HermiteGalsVelocityView wrongBegin(authority.identity(),{1,.001,.01,Position::End},[&](Vector2D q){return s.velocity(q,.01);});
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),wrongBegin,boundary,config);},"does not begin at source time");
        auto bad=config;bad.unitNormalTolerance=0;
        reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.arrival(),boundary,bad);},"invalid publication tolerance");
        check(source.state().identity()==id && &source.state().geometry().snapshot()==&scalar && source.state().primary().phaseVolumes.empty(),"failed step changed source");
    }
    auto mesh=makeMesh(8);Case s;s.pureFirstJet=true;s.motion=1;const auto source=initial(*mesh,s);
    const HermiteGalsTimeLayer time{1,0,.01,Position::End};HermiteGalsCandidate authority(*mesh,time,source.state().primary());
    const auto velocity=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),time,s.velocity({},0).value);
    const auto domain=traceDomain(*mesh,false);
    const HermiteGalsBoundaryFirstJetQuery boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){
        check(q.x==0 || q.x==1 || q.y==0 || q.y==1,"interior boundary reset");const auto j=s.exact(q,t.time());return HermiteGalsBoundaryFirstJet{j.v,{j.x,j.y}};};
    auto legacyCase=s;legacyCase.pureFirstJet=false;const auto legacy=initial(*mesh,legacyCase);
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,legacy,domain,velocity.previous(),velocity.arrival(),boundary);},"requires FirstJet source without inventory");
    const HermiteGalsVelocityView initialTarget(authority.identity(),{1,0,0,Position::Initial},[&](Vector2D q){return s.velocity(q,0);});
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),initialTarget,boundary);},"target must be Midpoint/End");
    auto wrongDomain=domain;wrongDomain.maxCoord.x=2;
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,wrongDomain,velocity.previous(),velocity.arrival(),boundary);},"domain must match");
    wrongDomain=domain;wrongDomain.cellSize={1,1};
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,wrongDomain,velocity.previous(),velocity.arrival(),boundary);},"cell scale");
    const auto mid=prepareHermiteGalsFirstJetStep(*mesh,source,domain,velocity.previous(),velocity.constantArrivalView(authority.identity(),{1,0,.01,Position::Midpoint}),boundary);
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,mid.publication(),domain,velocity.previous(),velocity.arrival(),boundary);},"source must be accepted Initial/End");
    const auto tooLarge=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),{1,0,1,Position::End},s.velocity({},0).value);
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,tooLarge.previous(),tooLarge.arrival(),boundary);},"sampled speed times full dt exceeds half a cell");
    const auto large=HermiteGalsTemporalVelocitySource::constant(source.state(),authority.identity(),{1,0,.0625,Position::End},{1,0});
    reject([&]{prepareHermiteGalsFirstJetStep(*mesh,source,domain,large.previous(),large.arrival(),boundary);},"interior Jacobian unavailable: boundary event");
    // This shear vanishes on horizontal master edges. Interior vertical edges
    // remain inside, while first-column centers actually hit the inflow wall.
    auto smallMesh=makeMesh(2);const auto smallSource=initial(*smallMesh,s);
    const auto shear=[](Vector2D q){const double y=q.y,f=y*(y-.5)*(y-1),df=3*y*y-3*y+.5;
        return LevelSetBoundaryTrace::VelocityJet{{256*f*f,0},{0,512*f*df,0,0}};};
    const HermiteGalsTemporalVelocitySource centerHit(smallSource.state(),authority.identity(),{1,0,4./9,Position::End},shear,shear);
    reject([&]{prepareHermiteGalsFirstJetStep(*smallMesh,smallSource,traceDomain(*smallMesh,false),centerHit.previous(),centerHit.arrival(),boundary);},
        "center prediction requires unresolved boundary data");
    std::cout<<"I12p3a FirstJet step PASS cases="<<cases<<" guards="<<guards<<'\n';
}
void firstJetMeshTransferContract()
{
    using Kind=HermiteGalsMeshTransferKind;
    int exactCases=0,projectionCases=0,guards=0;
    const auto reject=[&](auto f,const std::string& why){rejects(f,why);++guards;};
    for(bool radial:{false,true})for(bool amr:{false,true})for(bool nonpolynomial:{false,true}){
        auto sourceMesh=makeMesh(8,amr?1:0);Case s;s.pureFirstJet=true;s.radial=radial;s.amr=amr;s.nonpolynomial=nonpolynomial;if(radial)s.radius=.27;
        const auto source=initial(*sourceMesh,s);const auto id=source.state().identity();const auto& scalar=source.state().geometry().snapshot();
        std::optional<HermiteGalsPublication> refined;std::unique_ptr<QuadtreeMesh> refinedMesh;
        for(bool identical:{true,false}){
            auto target=makeMesh(8,identical?(amr?1:0):(amr?2:1));
            const auto result=transferHermiteGalsFirstJetMesh(*sourceMesh,source,*target,Kind::ExactRefinement);const auto& pub=result.publication();
            const auto& actual=pub.state().geometry().snapshot();const auto time=pub.state().timeLayer();
            check(result.sourceIdentity()==id && pub.state().identity()!=id && time.time()==0 && time.step==0
                && pub.state().primary().phaseVolumes.empty(),"exact regrid state/time authority");
            check(&pub.volumes().geometry().snapshot()==&actual && &pub.facePartitions().source()==&actual,"regrid has multiple scalar owners");
            double jet=0,trace=0,curvature=0;std::vector<long double> childVolumes(sourceMesh->cellIndexNum(),0);
            for(int c=0;c<target->cellIndexNum();++c){const auto b=actual.nativeCellBounds(c);const auto h=b[1]-b[0];const auto center=b[0]+h*.5;int parent=-1;
                for(int p:scalar.nativeCellCandidates(center,center)){const auto a=scalar.nativeCellBounds(p);
                    if(center.x>a[0].x && center.x<a[1].x && center.y>a[0].y && center.y<a[1].y){check(parent<0,"ambiguous exact parent");parent=p;}}
                check(parent>=0,"missing exact parent");childVolumes[parent]+=pub.volumes().phaseVolumes()[c].liquid;
                for(double x:{.25,.5,.75})for(double y:{.25,.5,.75}){const auto p=b[0]+Vector2D{x*h.x,y*h.y};
                    const auto a=scalar.scalarJetNative(p),z=actual.scalarJetAtNativeCellChart(c,{{x,false},{y,false}});
                    const double av[]={a.value,a.dx,a.dy,a.dxx,a.dxy,a.dyy},zv[]={z.value,z.dx,z.dy,z.dxx,z.dxy,z.dyy};
                    for(int k=0;k<6;++k)jet=std::max(jet,std::abs(zv[k]-av[k])/(1+std::abs(av[k])));}
                for(double t:{0.,.5,1.})for(Vector2D local:std::array<Vector2D,4>{{{0,t},{1,t},{t,0},{t,1}}}){
                    const auto a=scalar.scalarJetNative(b[0]+Vector2D{local.x*h.x,local.y*h.y}),z=actual.scalarJetAtNativeCellChart(c,
                        {{local.x,false},{local.y,false}});trace=std::max({trace,std::abs(a.value-z.value),std::abs(a.dx-z.dx),std::abs(a.dy-z.dy)});}
            }
            long double volume=0;for(const auto& c:source.volumes().cells())volume=std::max(volume,
                std::abs(childVolumes[c.cell]-source.volumes().phaseVolumes()[c.cell].liquid)/c.capacity);
            for(std::size_t c=0;c<pub.topology().curves.size();++c)for(std::size_t g=0;g<pub.topology().curves[c].graphs().size();++g)
                for(const auto& q:pub.quadrature(c,g)){const auto p=q.sample.frame.point;const auto j=scalar.scalarJetNative(p);
                    double reference=curvatureOf(j);if(radial)reference+=j.dx/(std::hypot(j.dx,j.dy)*p.x);
                    curvature=std::max(curvature,std::abs(q.sample.curvature-reference));}
            check(jet<=1e-10 && trace<=1e-12 && curvature<=1e-9 && volume<=1e-10,"exact regrid scalar/geometry gate");
            const auto& r=result.report();check(r.coarsenedCells==0 && r.refinedCells+r.unchangedCells==std::size_t(target->cellIndexNum()),"exact regrid leaf report");
            check((identical && r.refinedCells==0) || (!identical && r.refinedCells>0),"exact regrid fixture did not refine");
            std::cout<<"FIRST_JET_EXACT_REGRID,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<identical<<','<<jet<<','<<trace<<','<<curvature<<','<<volume<<",PASS\n";
            ++exactCases;
            if(!identical){refined=pub;refinedMesh=std::move(target);}
        }
        reject([&]{transferHermiteGalsFirstJetMesh(*refinedMesh,*refined,*sourceMesh,Kind::ExactRefinement);},"ExactRefinement cannot coarsen");
        const auto projected=transferHermiteGalsFirstJetMesh(*refinedMesh,*refined,*sourceMesh,Kind::FirstJetProjection);
        const auto& pub=projected.publication();const auto e=errors(pub,s);
        check(projected.report().coarsenedCells>0 && pub.state().timeLayer().time()==0 && pub.state().primary().phaseVolumes.empty(),"projection did not coarsen with zero time");
        if(!nonpolynomial)check(std::max({e.position,e.gradient,e.curvature})<=1e-10 && e.cellVolume<=5e-12,"Q2 coarsening reproduction");
        std::cout<<"FIRST_JET_PROJECT_REGRID,"<<radial<<','<<amr<<','<<nonpolynomial<<','<<e.position<<','<<e.gradient<<','<<e.curvature
            <<','<<e.volume<<','<<projected.report().relativeDomainVolumeChange<<",PASS\n";++projectionCases;
        reject([&]{transferHermiteGalsFirstJetMesh(*sourceMesh,source,*sourceMesh,static_cast<Kind>(-1));},"unknown transfer kind");
        reject([&]{transferHermiteGalsFirstJetMesh(*refinedMesh,source,*sourceMesh,Kind::FirstJetProjection);},"mesh");
        HermiteGalsPublicationConfig bad;bad.unitNormalTolerance=0;
        reject([&]{transferHermiteGalsFirstJetMesh(*sourceMesh,source,*sourceMesh,Kind::ExactRefinement,bad);},"invalid publication tolerance");
        auto moved=projected;const auto retained=std::move(moved);
        reject([&]{moved.publication();},"moved-from result");
        check(retained.publication().state().identity()==pub.state().identity() && source.state().identity()==id
            && &source.state().geometry().snapshot()==&scalar,"regrid failure mutated source/owner");
    }
    auto mesh=makeMesh(8);Case s;s.pureFirstJet=true;const auto source=initial(*mesh,s);
    auto nonnested=makeMesh(7);reject([&]{transferHermiteGalsFirstJetMesh(*mesh,source,*nonnested,Kind::FirstJetProjection);},"nested source/target leaves");
    DictionaryMesh d;d.Nx=d.Ny=8;d.maxInitialIterations=0;d.maxCoord={2,1};QuadtreeMesh foreign(d);
    reject([&]{transferHermiteGalsFirstJetMesh(*mesh,source,foreign,Kind::ExactRefinement);},"changes domain");
    auto legacyCase=s;legacyCase.pureFirstJet=false;const auto legacy=initial(*mesh,legacyCase);
    reject([&]{transferHermiteGalsFirstJetMesh(*mesh,legacy,*mesh,Kind::FirstJetProjection);},"requires FirstJet source without inventory");
    // A partial planar publication needs its actual outer boundary contacts.
    auto partial=source.state().primary();partial.geometryCells.pop_back();std::set<int> vertices;std::set<Key> sides;
    for(int c:partial.geometryCells){const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);for(int k=0;k<4;++k){
        const int v=mesh->verts(cell.cellToVert[k]).index();vertices.insert(v);Key key{v,mesh->verts(cell.cellToVert[(k+1)%4]).index()};
        const auto a=point(*mesh,key[0]),b=point(*mesh,key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.insert(key);}}
    partial.vertices.erase(std::remove_if(partial.vertices.begin(),partial.vertices.end(),[&](const auto& v){return !vertices.count(v.vertex);}),partial.vertices.end());
    partial.edges.erase(std::remove_if(partial.edges.begin(),partial.edges.end(),[&](const auto& e){return !sides.count(e.vertices);}),partial.edges.end());
    for(auto& v:partial.vertices){v.phi=point(*mesh,v.vertex).x-.53;v.dx=1;v.dy=0;}
    for(auto& e:partial.edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);e.phi=(a.x+b.x)*.5-.53;e.dn=a.y==b.y?0:1;}
    HermiteGalsCandidate candidate(*mesh,{},partial);const auto mixed=closeHermiteGalsMixedAmr(*mesh,candidate);
    std::vector<double> center;for(int c:partial.geometryCells)center.push_back(mesh->cells(mesh->cellIndexToLabel[c]).centerCoord().x-.53);
    HermiteGalsPublicationConfig partialConfig;partialConfig.volume.boundaryContacts={{.53,0},{.53,1}};
    const auto local=publishHermiteGalsFirstJetGeometry(*mesh,std::move(candidate),mixed,{candidate.identity(),center},partialConfig);
    reject([&]{transferHermiteGalsFirstJetMesh(*mesh,local,*mesh,Kind::ExactRefinement);},"requires full scalar support");
    // Coarsening must not silently lose a resolved compact drop.
    auto fine=makeMesh(32),coarse=makeMesh(4);
    const Profile compact=[](Vector2D p){const D2 x{p.x-.4375,1,0,0,0,0},y{p.y-.4375,0,1,0,0,0};
        const auto j=1-2*exponential((-1/(.04*.04))*(x*x+y*y));return Jet{j.v,j.x,j.y,j.xy};};
    const auto seed=stateFor(*fine,compact,false,{},[](int,Vector2D,Vector2D){return 0.;},{},HermiteGalsTransport::FirstJet);
    HermiteGalsCandidate small(*fine,{},seed.primary());const auto smallMixed=closeHermiteGalsMixedAmr(*fine,small);
    const auto smallDrop=publishHermiteGalsFirstJetGeometry(*fine,std::move(small),smallMixed,{small.identity(),centers(*fine,compact)});
    reject([&]{transferHermiteGalsFirstJetMesh(*fine,smallDrop,*coarse,Kind::FirstJetProjection);},"changed components or endpoint kinds");
    std::cout<<"I12p2b FirstJet mesh transfer PASS exact="<<exactCases<<" projection="<<projectionCases<<" guards="<<guards<<'\n';
}
void firstJetCrossing()
{
    int cases=0;long double maximumVolume=0;double maximumJet=0,maximumCurvature=0;
    for(int configuration=0;configuration<3;++configuration)for(bool amr:{false,true})
        for(int direction:{-1,1})for(int sign:{-1,1})for(double fraction:{.25,.75}){
        Case s;s.pureFirstJet=true;s.planar=true;s.radial=configuration==2;s.axialTranslation=configuration!=0;
        s.amr=amr;s.motion=1;s.speed=direction;s.planeSign=sign;const double dt=.02;s.planeOffset=.5-s.speed*fraction*dt;
        std::string stage="initialize";
        try{
            auto mesh=makeMesh(4,amr?1:0);const auto old=initial(*mesh,s);const auto id=old.state().identity();
            const HermiteGalsTimeLayer endTime{1,0,dt,Position::End};
            HermiteGalsCandidate velocityGeometry(*mesh,endTime,old.state().primary());
            const auto velocity=HermiteGalsTemporalVelocitySource::constant(old.state(),velocityGeometry.identity(),endTime,s.velocity({},0).value);
            const HermiteGalsTemporalCharacteristicSource characteristic(*mesh,traceDomain(*mesh,s.radial));
            const auto mid=advanceFirstJet(*mesh,old,s,{1,0,dt,Position::Midpoint},stage,&velocity,&characteristic),
                end=advanceFirstJet(*mesh,old,s,endTime,stage,&velocity,&characteristic);
            check(old.state().identity()==id && old.state().primary().phaseVolumes.empty(),"pure plane mutated old state");
            long double volumeError=0;double jetError=0,curvature=0;
            for(const auto* pub:{&mid,&end}){const auto t=pub->state().timeLayer().time();const double threshold=s.planeOffset+s.speed*t;
                for(const auto& m:pub->volumes().phaseVolumes()){
                    const auto b=pub->state().geometry().snapshot().nativeCellBounds(m.cell);const double low=s.axial()?b[0].y:b[0].x,high=s.axial()?b[1].y:b[1].x;
                    const long double length=std::max(0.,std::min(threshold,high)-low),cap=pub->domain().cellVolumes()[m.cell];
                    const long double measure=s.radial?pi*(static_cast<long double>(b[1].x)-b[0].x)*(static_cast<long double>(b[1].x)+b[0].x)
                        :s.axial()?b[1].x-b[0].x:b[1].y-b[0].y;
                    const auto reference=sign>0?length*measure:cap-length*measure;
                    volumeError=std::max(volumeError,std::abs(m.liquid-reference)/cap);
                }
                for(const auto& v:pub->state().primary().vertices){const auto j=s.exact(point(*mesh,v.vertex),t);
                    jetError=std::max({jetError,std::abs(v.phi-j.v),std::abs(v.dx-j.x),std::abs(v.dy-j.y)});}
                for(const auto& e:pub->state().primary().edges){const auto a=point(*mesh,e.vertices[0]),b=point(*mesh,e.vertices[1]);const auto j=s.exact(a+(b-a)*.5,t);
                    jetError=std::max({jetError,std::abs(e.phi-j.v),std::abs(e.dn-(a.y==b.y?j.y:j.x))});}
                for(std::size_t c=0;c<pub->topology().curves.size();++c)for(std::size_t g=0;g<pub->topology().curves[c].graphs().size();++g)
                    for(const auto& q:pub->quadrature(c,g)){const auto normal=q.sample.frame.normal;
                        jetError=std::max({jetError,std::abs(normal.x-(s.axial()?0:sign)),std::abs(normal.y-(s.axial()?sign:0))});
                        curvature=std::max(curvature,std::abs(q.sample.curvature));}
            }
            check(mid.state().timeLayer().time()==dt/2 && end.state().timeLayer().time()==dt,"pure plane stage times");
            check(volumeError<=5e-12 && jetError<=1e-11 && curvature<=1e-10,"pure plane analytic gate");
            maximumVolume=std::max(maximumVolume,volumeError);maximumJet=std::max(maximumJet,jetError);maximumCurvature=std::max(maximumCurvature,curvature);++cases;
            std::cout<<"FIRST_JET_CROSSING,"<<configuration<<','<<amr<<','<<direction<<','<<sign<<','<<fraction
                <<','<<volumeError<<','<<jetError<<','<<curvature<<",PASS\n";
        }catch(const std::exception& e){++failures;std::cout<<"FIRST_JET_CROSSING_FAILURE,"<<configuration<<','<<amr<<','<<direction<<','<<sign<<','<<fraction<<','<<stage<<','<<e.what()<<'\n';}
    }
    std::cout<<"I12p2a FirstJet crossing "<<(failures?"FAIL":"PASS")<<" cases="<<cases<<" planned=48 failures="<<failures
        <<" max_volume="<<maximumVolume<<" max_jet="<<maximumJet<<" max_curvature="<<maximumCurvature<<'\n';
}
}
int main(int argc,char** argv) try {
    if(argc!=2) throw std::invalid_argument("publication | step | mesh-transfer | crossing");
    const std::string mode=argv[1];
    if(mode=="publication") firstJetPublicationContract();
    else if(mode=="step") firstJetStepContract();
    else if(mode=="mesh-transfer") firstJetMeshTransferContract();
    else if(mode=="crossing") firstJetCrossing();
    else throw std::invalid_argument("unknown contract");
    return failures?1:0;
} catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}