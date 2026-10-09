#include "Interface/InterfaceTransport.h"
#include "Interface/IO/InterfaceFirstJetOutput.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include "Modules/Interface/HermiteGalsMeshTransfer.h"
#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Core/Geometry/ActiveQuarticHermiteSkeleton.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <set>

using Config=InterfaceTransportConfig;
InterfaceTransportConfigError::InterfaceTransportConfigError(Config::Field f,const char* why)
    :std::invalid_argument(std::string("[Interface FirstJet App] ")+why),field(f){}
void Config::validate() const
{
    const auto check=[](bool ok,Field f,const char* why){if(!ok)throw InterfaceTransportConfigError(f,why);};
    check(resolution>=4 && resolution<=256,Field::Resolution,"resolution must be in [4,256]");
    check(steps>=1 && steps<=1000000,Field::Steps,"steps must be in [1,1000000]");
    check(std::isfinite(dt) && dt>0 && std::isfinite(dt*steps),Field::Dt,"dt and total time must be finite and positive");
    check(std::isfinite(endTime) && endTime>0 && std::abs(dt*steps-endTime)<=1e-12*endTime,
        Field::EndTime,"endTime must equal the configured integer number of dt steps");
    check(!writeFields || (writeEverySteps>=1 && writeEverySteps<=1000000),Field::WriteEverySteps,"writeEverySteps must be in [1,1000000]");
    int previous=0;for(int step:regridSteps){check(step>previous && step<steps,Field::RegridSteps,"Regrid.atSteps must increase strictly before the final step");previous=step;}
    check(std::isfinite(domainMin.x) && std::isfinite(domainMin.y) && std::isfinite(domainMax.x) && std::isfinite(domainMax.y)
        && domainMax.x>domainMin.x && domainMax.y>domainMin.y && std::isfinite(domainMax.x-domainMin.x) && std::isfinite(domainMax.y-domainMin.y)
        && (!radial || domainMin.x==0),Field::Domain,"finite ordered bounds; rz requires radial minimum 0");
    if(initialAmr || !regridSteps.empty())check(refinementLevel>=1 && refinementLevel<=2 && std::isfinite(refineMin.x) && std::isfinite(refineMin.y)
        && std::isfinite(refineMax.x) && std::isfinite(refineMax.y) && refineMin.x>=domainMin.x && refineMin.y>=domainMin.y
        && refineMax.x<=domainMax.x && refineMax.y<=domainMax.y && refineMax.x>refineMin.x && refineMax.y>refineMin.y,
        Field::Refinement,"ordered refinement box inside Domain; level in [1,2]");
}
InterfacePrescribedData::InterfacePrescribedData(const Config& c)
    :initial_(c.initialPhi,c.constants,c.radial?Quantum::ScalarExpression::Coordinates::Meridional:Quantum::ScalarExpression::Coordinates::Planar),
     origin_(c.velocityOrigin),scale_(c.velocityScale)
{
    c.validate();using Expression=Quantum::ScalarExpression;
    if(c.prescribedBoundary)boundary_=std::make_shared<const Expression>(c.boundaryPhi,c.constants,
        c.radial?Expression::Coordinates::MeridionalTime:Expression::Coordinates::PlanarTime);
    if(c.velocity==Config::Velocity::Polynomial){
        if(!(std::isfinite(scale_.x) && std::isfinite(scale_.y) && scale_.x>0 && scale_.y>0))throw std::invalid_argument("invalid polynomial scales");
        for(int k=0;k<16;++k){coefficients_[k][0]=std::make_shared<const Expression>(c.coefficients[k].x,c.constants,Expression::Coordinates::Time);
            coefficients_[k][1]=std::make_shared<const Expression>(c.coefficients[k].y,c.constants,Expression::Coordinates::Time);}
    }
}
Quantum::ScalarExpression::Jet InterfacePrescribedData::initial(Vector2D p) const{return initial_.sample(p.x,p.y);}
Quantum::ScalarExpression::Jet InterfacePrescribedData::boundary(Vector2D p,double time) const
{if(!boundary_)throw std::logic_error("no prescribed boundary data");return boundary_->sample(p.x,p.y,time);}
Quantum::Geometry::PolynomialVelocityData InterfacePrescribedData::polynomial(double previous,double arrival) const
{
    Quantum::Geometry::PolynomialVelocityData data;data.origin=origin_;data.scale=scale_;
    for(int k=0;k<16;++k){if(!coefficients_[k][0])throw std::logic_error("no polynomial velocity data");
        data.previous[k]={coefficients_[k][0]->sample(previous).value,coefficients_[k][1]->sample(previous).value};
        data.arrival[k]={coefficients_[k][0]->sample(arrival).value,coefficients_[k][1]->sample(arrival).value};}
    return data;
}
namespace
{
using namespace Quantum::Interface;
using Position=HermiteGalsTimeLayer::Position;
using Key=std::array<int,2>;
void require(bool ok,const char* why){if(!ok)throw std::invalid_argument(std::string("[Interface FirstJet App] ")+why);}
std::unique_ptr<QuadtreeMesh> makeMesh(const Config& c,bool amr)
{
    DictionaryMesh input;input.minCoord=c.domainMin;input.maxCoord=c.domainMax;input.Nx=input.Ny=c.resolution;input.maxInitialIterations=0;
    input.maxCellNum=std::max(16384,16*c.resolution*c.resolution);
    auto mesh=std::make_unique<QuadtreeMesh>(input);
    if(amr){DictionaryMesh::RefineBox box;box.level=c.refinementLevel;box.min=c.refineMin;box.max=c.refineMax;mesh->applyRefineBoxes({box});}
    return mesh;
}
HermiteGalsPublication initialize(const QuadtreeMesh& mesh,const Config& c,const InterfacePrescribedData& data)
{
    const auto point=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
    Quantum::Geometry::ActiveQuarticHermiteSkeletonInput input;
    input.activeCells.assign(mesh.cellIndexNum(),1);input.nodes.resize(mesh.vertIndexNum());
    std::set<Key> sides;
    for(int c=0;c<mesh.cellIndexNum();++c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        for(int k=0;k<4;++k){Key key{mesh.verts(cell.cellToVert[k]).index(),mesh.verts(cell.cellToVert[(k+1)%4]).index()};
            const auto a=point(key[0]),b=point(key[1]);if(a.y==b.y?a.x>b.x:a.y>b.y)std::swap(key[0],key[1]);sides.insert(key);}}
    for(const auto& key:sides)input.edgeMidpoints.push_back({key,0,0});
    // Use the existing master/hanging owner engine; the zero data are not a
    // geometry or an initial scalar and never leave this topology query.
    const Quantum::Geometry::ActiveQuarticHermiteSkeleton topology(mesh,input,c.radial);
    HermiteGalsPrimaryState primary;primary.transport=HermiteGalsTransport::FirstJet;
    primary.measure=c.radial?HermiteGalsMeasure::AxisymmetricVolume:HermiteGalsMeasure::CartesianArea;
    const Key none{-1,-1};
    for(const auto& v:topology.vertices())if(v.masterSide==none){const auto j=data.initial(point(v.index));primary.vertices.push_back({v.index,j.value,j.derivative,j.partialSecond});}
    for(const auto& e:topology.edges())if(e.parentSide==none){const auto a=point(e.vertices[0]),b=point(e.vertices[1]);const auto j=data.initial(a+(b-a)*.5);
        primary.edges.push_back({e.vertices,j.value,a.y==b.y?j.partialSecond:j.derivative});}
    std::sort(primary.vertices.begin(),primary.vertices.end(),[](const auto& a,const auto& b){return a.vertex<b.vertex;});
    std::sort(primary.edges.begin(),primary.edges.end(),[](const auto& a,const auto& b){return a.vertices<b.vertices;});
    std::vector<double> centers;
    for(int c=0;c<mesh.cellIndexNum();++c){primary.geometryCells.push_back(c);centers.push_back(data.initial(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()).value);}
    HermiteGalsCandidate candidate(mesh,{},primary);const auto mixed=closeHermiteGalsMixedAmr(mesh,candidate);
    const HermiteGalsCenterPrediction center{candidate.identity(),centers};
    return publishHermiteGalsFirstJetGeometry(mesh,std::move(candidate),mixed,center);
}
LevelSetBoundaryTrace::Domain traceDomain(const QuadtreeMesh& mesh,bool radial)
{
    Vector2D h=mesh.maxCoord-mesh.minCoord;for(int c=0;c<mesh.cellIndexNum();++c){const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();h.x=std::min(h.x,size.x);h.y=std::min(h.y,size.y);}
    return {mesh.minCoord,mesh.maxCoord,h,radial};
}
struct Accepted
{
    std::unique_ptr<QuadtreeMesh> mesh;
    HermiteGalsPublication publication;
    bool amr;
};
}

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
        auto mesh=makeMesh(c,c.initialAmr);auto initial=initialize(*mesh,c,data);
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
            const HermiteGalsTemporalCharacteristicSource characteristic(*accepted->mesh,traceDomain(*accepted->mesh,c.radial));
            HermiteGalsBoundaryFirstJetQuery boundary;
            if(c.prescribedBoundary)boundary=[&](Vector2D q,HermiteGalsTimeLayer t,Boundary){const auto j=data.boundary(q,t.time());return HermiteGalsBoundaryFirstJet{j.value,{j.derivative,j.partialSecond}};};
            const auto prepared=prepareHermiteGalsFirstJetStep(*accepted->mesh,accepted->publication,characteristic.domain(),velocity.previous(),velocity.arrival(),boundary);
            auto next=prepared.publication();
            const auto stepSeconds=std::chrono::duration<double>(Clock::now()-stepStart).count();
            double transferSeconds=0;
            std::unique_ptr<Accepted> staged;std::optional<HermiteGalsMeshTransferReport> regrid;
            if(std::binary_search(c.regridSteps.begin(),c.regridSteps.end(),step)){
                const auto transferStart=Clock::now();
                stage="mesh transfer";auto target=makeMesh(c,!accepted->amr);
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
