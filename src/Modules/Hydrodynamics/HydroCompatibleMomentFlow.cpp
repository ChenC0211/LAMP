#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTraction.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMassCurrent.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Hydrodynamics/Detail/HydroCompatibleMomentTransfer.h"
#include "Modules/Hydrodynamics/Detail/HydroCompatibleMomentBlockSolve.h"
#include "Modules/Hydrodynamics/Detail/HydroCompatibleMomentSplineBasis.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/SuperLUSupport>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <iostream>

namespace Quantum::Hydrodynamics {
namespace {
using Matrix=HydroCompatibleMomentSpace::Matrix;
// Optional wall-clock observations only; never enter a physical decision.
struct MomentWorkTimer {
    using Clock=std::chrono::steady_clock;
    bool enabled=std::getenv("QUANTUM_MOMENT_STEP_PROFILE")!=nullptr;
    Clock::time_point last=Clock::now();
    const char* group;
    std::uint64_t step;
    MomentWorkTimer(const char* name,std::uint64_t index):group(name),step(index){}
    void mark(const char* work){
        if(!enabled)return;
        const auto now=Clock::now();std::ostringstream out;
        out<<"[moment_step_profile] step="<<step<<" group="<<group<<" work="<<work
            <<" seconds="<<std::chrono::duration<double>(now-last).count()<<'\n';
        std::cout<<out.str()<<std::flush;last=now;
    }
};
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
bool sameTime(double a,double b){return std::abs(a-b)<=16*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(a),std::abs(b)});}
double divergenceResidual(const Matrix& d,const Eigen::VectorXd& u){
    const Eigen::VectorXd residual=d*u,scale=d.cwiseAbs()*u.cwiseAbs();
    const double denominator=scale.stableNorm(),numerator=residual.stableNorm();
    return denominator>0?numerator/denominator:(numerator==0?0:std::numeric_limits<double>::infinity());
}
std::shared_ptr<const HydroCompatibleMomentSpline::NativeMaterialForms> evaluatedNativeForms(
    HydroCompatibleMomentSpline::NativeMaterialForms forms,const HydroCompatibleMomentSpline::NativeVelocity& velocity){
    const auto& q=velocity.coordinates();const int n=int(q.size());
    forms.kinetic=.5*q.dot(forms.mass.topLeftCorner(n,n)*q);forms.dissipation=q.dot(forms.viscosity.topLeftCorner(n,n)*q);
    const Eigen::Vector2d momentum=forms.momentum.leftCols(n)*q;forms.physicalMomentum={momentum.x(),momentum.y()};
    return std::make_shared<const HydroCompatibleMomentSpline::NativeMaterialForms>(std::move(forms));
}
}
struct HydroCompatibleMomentFlow::PhaseData {
    std::shared_ptr<const InterfaceState> interface;
    HydroCompatibleMomentMaterial material,previousMaterial,midpointMaterial;
    HydroCompatibleMomentViscosity stress;
    std::vector<int> activeTangents;
    std::optional<HydroCompatibleMomentTraction::Result> traction;
    std::shared_ptr<const HydroCompatibleMomentSpline::NativeMaterialForms> nativeMaterial,nativePreviousMaterial,nativeMidpointMaterial;
    double length=0;
};
struct HydroCompatibleMomentFlow::Impl {
    HydroCompatibleMomentSpace space;
    HydroCompatibleMomentViscosity stress;
    HydroCompatibleMomentBoundary boundary;
    std::vector<HydroBoundaryPatch> boundaryPatches;
    std::shared_ptr<const State> state;
    std::shared_ptr<const int> owner=std::make_shared<const int>(0);
    std::uint64_t revision=0;
    double rho=0;
    double surfaceTension=0;
    double initialCorrection=0;
    Matrix mass;
    std::shared_ptr<const PhaseData> phase;
    std::optional<HydroCompatibleMomentComplement> complement;
    Matrix nativeBoundary,allNativeBoundary;
    Eigen::VectorXd nativeBoundaryValues,allNativeBoundaryValues;
    std::optional<bool> fixedPolicy;
    int nb=0,nv=0,np=0;
    std::vector<int> normalBoundary;
    std::vector<bool> fixedBoundary;
    Eigen::VectorXd gauge;
    HydroCompatibleMomentBlockSolve solver;
    std::optional<double> factorDt;
    void validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const{space.validate(mesh,domain);stress.validate(mesh,domain,space);}
    void configureGauge(const QuadtreeMesh& mesh,const NativeLeafDomain& domain){
        gauge=Eigen::VectorXd::Zero(np);
        for(int c=0;c<mesh.cellIndexNum();c++){gauge[3*c]=domain.cellVolumes()[domain.cellRowByIndex()[c]];
            if(domain.axisymmetric()){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);gauge[3*c+1]=gauge[3*c]*cell.size().x/(12*cell.centerCoord().x);}}
    }
    void configureNativeBoundary(){
        const auto& velocity=*state->physicalVelocity;const auto& h=velocity.kinematicConstraints();
        const Matrix normal=velocity.normalBoundaryConstraints(boundary);
        const int count=int(h.rows()+normal.rows()),nq=int(velocity.coordinates().size());
        std::vector<Eigen::Triplet<double>> entries;
        const auto append=[&](const Matrix& source,int offset){
            for(int col=0;col<source.outerSize();col++)for(Matrix::InnerIterator it(source,col);it;++it)
                if(it.value()!=0)entries.emplace_back(offset+it.row(),col,it.value());
        };
        append(h,0);append(normal,int(h.rows()));
        allNativeBoundary.resize(count,nq);allNativeBoundary.setFromTriplets(entries.begin(),entries.end());
        allNativeBoundaryValues=Eigen::VectorXd::Zero(count);
        allNativeBoundaryValues.tail(normal.rows())=boundary.normalConstraintValues();
        Eigen::VectorXd scaling=Eigen::VectorXd::Zero(count);
        for(int col=0;col<nq;col++)for(Matrix::InnerIterator it(allNativeBoundary,col);it;++it)scaling[it.row()]+=it.value()*it.value();
        for(int row=0;row<count;row++)scaling[row]=scaling[row]>0?1/std::sqrt(scaling[row]):1;
        allNativeBoundary=scaling.asDiagonal()*allNativeBoundary;allNativeBoundaryValues=scaling.asDiagonal()*allNativeBoundaryValues;
        // Rank-reveal the original normalized boundary rows, retaining the
        // same pivoting and threshold. Interior field coordinates have zero
        // boundary trace; omit only those zero columns from the dense QR.
        // The evolved full field and the saddle constraints stay sparse.
        std::vector<int> support(nq,-1);int active=0;
        for(int col=0;col<nq;col++)for(Matrix::InnerIterator it(allNativeBoundary,col);it;++it)
            if(it.value()!=0){support[col]=active++;break;}
        std::vector<int> selected(count,-1);int rank=0;
        if(active&&count){
            Eigen::MatrixXd traces=Eigen::MatrixXd::Zero(active,count);
            for(int col=0;col<nq;col++)if(support[col]>=0)
                for(Matrix::InnerIterator it(allNativeBoundary,col);it;++it)traces(support[col],it.row())=it.value();
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(traces);qr.setThreshold(1e-10);
            rank=int(qr.rank());for(int k=0;k<rank;k++)selected[qr.colsPermutation().indices()[k]]=k;
        }
        entries.clear();nativeBoundaryValues.resize(rank);
        for(int row=0;row<count;row++)if(selected[row]>=0)nativeBoundaryValues[selected[row]]=allNativeBoundaryValues[row];
        for(int col=0;col<nq;col++)for(Matrix::InnerIterator it(allNativeBoundary,col);it;++it)
            if(selected[it.row()]>=0)entries.emplace_back(selected[it.row()],col,it.value());
        nativeBoundary.resize(rank,nq);nativeBoundary.setFromTriplets(entries.begin(),entries.end());
        if(std::getenv("QUANTUM_MOMENT_SETUP_PROFILE"))std::cout<<"MOMENT_BOUNDARY_SETUP field_dofs="<<nq
            <<" supported_dofs="<<active<<" rows="<<count<<" rank="<<rank<<" nonzeros="<<allNativeBoundary.nonZeros()<<std::endl;
        require((allNativeBoundary*velocity.coordinates()-allNativeBoundaryValues).norm()<2e-10,"complete initial field changed prescribed boundaries");
    }
    bool movingMaterial() const noexcept{return phase&&bool(phase->interface);}
    void initializeMaterial(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,HydroCompatibleMomentMaterial material,
        std::shared_ptr<const InterfaceState> interface={},double sigma=0){
        auto initial=std::make_shared<State>(*state);
        initial->physicalVelocity=HydroCompatibleMomentSpline::initializeNative(mesh,domain,space,state->velocity,state->time,&boundary);
        initial->previousPhysicalVelocity=initial->physicalVelocity;
        initial->velocity=initial->physicalVelocity->nativeMoments();initial->previousVelocity=initial->velocity;
        state=std::move(initial);complement=HydroCompatibleMomentComplement::build(mesh,domain,space,boundary.fixed());configureNativeBoundary();
        auto current=std::make_shared<PhaseData>();current->interface=std::move(interface);current->material=std::move(material);
        current->previousMaterial=current->midpointMaterial=current->material;current->stress=stress;
        current->nativeMaterial=std::make_shared<const HydroCompatibleMomentSpline::NativeMaterialForms>(
            state->physicalVelocity->mixedMaterialForms(mesh,domain,space,current->material,*complement));
        current->nativePreviousMaterial=current->nativeMidpointMaterial=current->nativeMaterial;
        rho=current->material.uniformDensity().value_or(std::numeric_limits<double>::quiet_NaN());
        if(sigma>0)current->length=HydroCompatibleMomentTraction::measure(*current->material.boundaryHandle(),mesh,domain).length;
        surfaceTension=sigma;phase=std::move(current);
    }
    std::shared_ptr<const PhaseData> preparePhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
        const HydroStepTime& time,const Eigen::VectorXd& wind,Options& options,bool& spaceEnriched) const {
        MomentWorkTimer timer("phase",time.index);
        // Optional failure checkpoint: the exact candidate wind and accepted
        // material geometry permit local replay without another coupled run.
        // It is observational and never supplies physics or acceptance data.
        const auto checkpoint=[&](const std::string& path,const Geometry::Parametric::Curve& curve,double geometryTime,const char* format){
            std::ofstream out(path);require(bool(out),"cannot write requested moment phase checkpoint");
            const auto h=mesh.cells(mesh.cellIndexToLabel[0]).size();
            require(!curve.boundary.fixed(0)&&!curve.boundary.fixed(1)&&curve.boundary.corners.empty(),"V1 phase checkpoint cannot represent supplied endpoints and ports");
            out<<std::setprecision(17)<<format<<'\n'<<domain.axisymmetric()<<' '<<domain.uniformLeaves()<<' '<<mesh.cellIndexNum()<<' '<<wind.size()<<' '
                <<time.index<<' '<<geometryTime<<' '<<time.dt<<'\n'<<mesh.minCoord.x<<' '<<mesh.minCoord.y<<' '<<mesh.maxCoord.x<<' '<<mesh.maxCoord.y<<' '<<h.x<<' '<<h.y<<'\n';
            for(int i=0;i<wind.size();i++)out<<wind[i]<<'\n';out<<curve.boundary.axisPair()<<' '<<curve.arcs.size()<<'\n';
            for(const auto& arc:curve.arcs){out<<arc.a<<' '<<arc.b;for(const auto& p:arc.p)out<<' '<<p.x<<' '<<p.y;out<<'\n';}
            out<<options.interfaceProjection.minimumMaterialKnots.size()<<'\n';for(double s:options.interfaceProjection.minimumMaterialKnots)out<<s<<'\n';
            require(bool(out),"moment phase checkpoint write failed");
        };
        if(const char* path=std::getenv("QUANTUM_MOMENT_PHASE_CHECKPOINT")){
            require(phase->interface->parametricGeometry().componentCount()==1&&phase->interface->parametricGeometry().retiredComponents().empty(),
                "V2 phase checkpoint requires one active component without retired remainders");
            checkpoint(path,phase->interface->parametricGeometry().geometry(),time.begin,"MOMENT_PHASE_COMPLETE_CANDIDATE_V2");
        }
        require(state->physicalVelocity&&complement,"moving phase requires the complete owned velocity");
        const auto velocity=state->physicalVelocity->withCoordinates(wind,time.begin+.5*time.dt);
        timer.mark("velocity");
        const auto map=HydroInterfaceTransport::frozenParametricFlow(velocity);const auto flux=velocity.normalFlux();
        const InterfaceStepTime half{time.index,time.begin,time.begin+.5*time.dt,.5*time.dt},whole{time.index,time.begin,time.begin+time.dt,time.dt};
        auto middle=InterfaceParametricMotion::prepare(mesh,domain,*phase->interface,half,half.dt,map,flux,options.interfaceProjection);
        timer.mark("half_motion");
        auto next=InterfaceParametricMotion::prepare(mesh,domain,*phase->interface,whole,whole.dt,map,flux,options.interfaceProjection);
        timer.mark("full_motion");
        for(;;){
            auto& retained=options.interfaceProjection.minimumMaterialKnots;auto common=retained;
            for(const auto* state:{&middle.state,&next.state})for(const auto& c:state->parametricGeometry().boundaryHandle()->components())for(const auto& arc:c.curve->arcs)common.push_back(arc.a);
            common.push_back(1);std::sort(common.begin(),common.end());common.erase(std::unique(common.begin(),common.end()),common.end());
            if(common==retained)break;
            // Adapt the approximation space outside the midpoint fixed point.
            // Half/end geometry share a nested partition; no knot can toggle
            // between nonlinear evaluations of the same physical time step.
            retained=std::move(common);spaceEnriched=true;
            InterfaceParametricMotion::refineMaterialSpace(middle,retained);
            InterfaceParametricMotion::refineMaterialSpace(next,retained);
        }
        // Interface now supplies domain-restricted projection, outlet phase
        // and surface transport, and explicit retirement. Admission must still
        // validate projected curves against the actual outlet patches and
        // commit retirement with material/history owners; retain rejection
        // until that common boundary acceptance is wired.
        // Check BOTH layers: a midpoint-only traction guard can otherwise
        // admit an unsupported full endpoint that crossed after mid-step.
        const Geometry::Parametric::DomainRectangle physicalBox{mesh.minCoord.x,mesh.maxCoord.x,mesh.minCoord.y,mesh.maxCoord.y};
        for(const auto* layer:{&middle.state,&next.state})for(const auto& c:layer->parametricGeometry().boundaryHandle()->components())
            Geometry::Parametric::validateContainedSurface(*c.curve,physicalBox,domain.axisymmetric());
        auto result=std::make_shared<PhaseData>();const auto properties=phase->material.properties();
        timer.mark("material_space");
        result->previousMaterial=phase->material;
        result->nativePreviousMaterial=phase->nativeMaterial;
        result->midpointMaterial=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,space,middle.state.parametricGeometry().boundaryHandle(),half.end,properties);
        result->material=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,space,next.state.parametricGeometry().boundaryHandle(),whole.end,properties);
        timer.mark("cut_material");
        result->stress=HydroCompatibleMomentViscosity::build(mesh,domain,space,result->midpointMaterial.viscosityMoments());
        timer.mark("viscosity");
        if(surfaceTension>0){
            result->traction=HydroCompatibleMomentTraction::build(mesh,domain,space,velocity,*complement,result->midpointMaterial,
                [&](const Geometry::Parametric::SurfaceSample& q){return HydroCompatibleMomentTraction::Components{surfaceTension*q.curvature,0};},surfaceTension);
            result->length=HydroCompatibleMomentTraction::measure(*result->material.boundaryHandle(),mesh,domain).length;
        }
        timer.mark("traction");
        // Opt-in material-only replay input. These are the actual projected
        // trial layers AFTER the independent inventory constraints and common
        // material-space preparation, not reconstructions of accepted curves.
        if(const char* prefix=std::getenv("QUANTUM_MOMENT_MATERIAL_CHECKPOINT")){
            for(const auto& layer:std::array<std::pair<const HydroCompatibleMomentMaterial*,const char*>,2>{{{&result->midpointMaterial,"half"},{&result->material,"full"}}}){
                require(layer.first->boundaryHandle()->components().size()==1,"material replay snapshot requires one active component");
                const auto& curve=*layer.first->geometryHandle();
                checkpoint(std::string(prefix)+"."+layer.second+".dat",curve,layer.first->geometryTime(),"MOMENT_MATERIAL_COMPLETE_SNAPSHOT_V1");
                std::cout<<"[moment_material_snapshot] step="<<time.index<<" layer="<<layer.second
                    <<" geometry_time="<<std::setprecision(17)<<layer.first->geometryTime()<<" arcs="<<curve.arcs.size()<<std::endl;
            }
        }
        result->nativeMidpointMaterial=std::make_shared<const HydroCompatibleMomentSpline::NativeMaterialForms>(velocity.mixedMaterialForms(mesh,domain,space,result->midpointMaterial,*complement));
        result->nativeMaterial=std::make_shared<const HydroCompatibleMomentSpline::NativeMaterialForms>(velocity.withCoordinates(wind,whole.end).mixedMaterialForms(mesh,domain,space,result->material,*complement));
        if(timer.enabled)std::cout<<"[native_material_cache] step="<<time.index<<" cells="<<mesh.cellIndexNum()
            <<" half_reused="<<result->nativeMidpointMaterial->reusedPureCells<<" whole_reused="<<result->nativeMaterial->reusedPureCells
            <<" half_points="<<result->nativeMidpointMaterial->quadraturePoints<<" whole_points="<<result->nativeMaterial->quadraturePoints
            <<" prepared_norms="<<result->nativeMidpointMaterial->preparedNormCells+result->nativeMaterial->preparedNormCells<<std::endl;
        timer.mark("complete_material");
        result->interface=std::make_shared<const InterfaceState>(std::move(next.state));return result;
    }
};
HydroCompatibleMomentFlow::HydroCompatibleMomentFlow(std::unique_ptr<Impl> p):impl_(std::move(p)){}
HydroCompatibleMomentFlow::~HydroCompatibleMomentFlow()=default;
HydroCompatibleMomentFlow::HydroCompatibleMomentFlow(HydroCompatibleMomentFlow&&) noexcept=default;
HydroCompatibleMomentFlow& HydroCompatibleMomentFlow::operator=(HydroCompatibleMomentFlow&&) noexcept=default;
const HydroCompatibleMomentFlow::State& HydroCompatibleMomentFlow::state() const{require(bool(impl_),"moved moment owner");return *impl_->state;}
double HydroCompatibleMomentFlow::initialBoundaryCorrection() const{require(bool(impl_),"moved moment owner");return impl_->initialCorrection;}
const HydroCompatibleMomentBoundary& HydroCompatibleMomentFlow::boundary() const{require(bool(impl_),"moved moment owner");return impl_->boundary;}
const HydroCompatibleMomentSpace& HydroCompatibleMomentFlow::space() const{require(bool(impl_),"moved moment owner");return impl_->space;}
const HydroCompatibleMomentViscosity& HydroCompatibleMomentFlow::viscosity() const{require(bool(impl_),"moved moment owner");return impl_->phase?impl_->phase->stress:impl_->stress;}
const std::vector<int>& HydroCompatibleMomentFlow::activeTangentialRows() const{require(impl_&&impl_->phase,"moved moment owner");return impl_->phase->activeTangents;}
double HydroCompatibleMomentFlow::density() const{require(bool(impl_)&&std::isfinite(impl_->rho),"moment material has no uniform density");return impl_->rho;}
const Matrix& HydroCompatibleMomentFlow::kineticMass() const{require(bool(impl_),"moved moment owner");return impl_->phase?impl_->phase->material.kineticMass():impl_->mass;}
double HydroCompatibleMomentFlow::kineticEnergy() const{
    require(bool(impl_),"moved moment owner");const auto& p=*impl_;
    require(p.phase&&p.state->physicalVelocity,"missing complete kinetic state");
    const auto& q=p.state->physicalVelocity->coordinates();return .5*q.dot(p.phase->nativeMaterial->mass.topLeftCorner(q.size(),q.size())*q);
}
const InterfaceState& HydroCompatibleMomentFlow::interfaceState() const{require(impl_&&impl_->movingMaterial(),"moment owner has no moving interface");return *impl_->phase->interface;}
HydroCompatibleMomentSpline::NativeCellBalance HydroCompatibleMomentFlow::physicalCellBalance(VelocityLevel level) const{
    require(impl_&&impl_->phase&&impl_->state->physicalVelocity,"physical cell balance requires the complete material state");
    const auto& p=*impl_;const auto& s=*p.state;
    if(level==VelocityLevel::Endpoint)return s.physicalVelocity->cellBalance(*p.phase->nativeMaterial);
    if(level==VelocityLevel::PreviousEndpoint){require(s.acceptedSteps>0&&s.previousPhysicalVelocity&&p.phase->nativePreviousMaterial,"physical cell balance has no previous endpoint");return s.previousPhysicalVelocity->cellBalance(*p.phase->nativePreviousMaterial);}
    require(level==VelocityLevel::Midpoint&&s.stageTime&&s.midpointPhysicalVelocity&&p.phase->nativeMidpointMaterial,"physical cell balance has no solved midpoint");
    return s.midpointPhysicalVelocity->cellBalance(*p.phase->nativeMidpointMaterial);
}
std::array<HydroCompatibleMomentFlow::PhaseIntegral,2> HydroCompatibleMomentFlow::phaseIntegrals(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,VelocityLevel level) const{
    require(impl_&&impl_->phase&&impl_->state->physicalVelocity,"phase integrals require the complete material state");
    const auto& p=*impl_;p.validate(mesh,domain);const auto& state=*p.state;
    const auto* velocity=&*state.physicalVelocity;const auto* forms=p.phase->nativeMaterial.get();
    if(level==VelocityLevel::PreviousEndpoint){
        require(state.acceptedSteps>0&&state.previousPhysicalVelocity&&p.phase->nativePreviousMaterial,"phase integrals have no previous endpoint");
        velocity=&*state.previousPhysicalVelocity;forms=p.phase->nativePreviousMaterial.get();
    }else if(level==VelocityLevel::Midpoint){
        require(state.stageTime&&state.midpointPhysicalVelocity&&p.phase->nativeMidpointMaterial,"phase integrals have no solved midpoint");
        velocity=&*state.midpointPhysicalVelocity;forms=p.phase->nativeMidpointMaterial.get();
    }else require(level==VelocityLevel::Endpoint,"unknown phase integral time level");
    const auto& q=velocity->coordinates();const auto origin=mesh.minCoord;
    std::array<std::array<long double,3>,2> first{};
    for(const auto& f:forms->cells){const auto& cell=mesh.cells(mesh.cellIndexToLabel[f.index]);const auto h=cell.size(),c=cell.centerCoord()-origin;
        for(int phase=0;phase<2;phase++){const auto& m=phase==0?f.liquid:f.gas;
            first[phase][0]+=m[0];first[phase][1]+=c.x*m[0]+h.x*m[1];first[phase][2]+=c.y*m[0]+h.y*m[2];}
    }
    std::array<PhaseIntegral,2> result;std::array<Vector2D,2> localMean;
    const double undefined=std::numeric_limits<double>::quiet_NaN();
    for(int phase=0;phase<2;phase++){
        auto& r=result[phase];r.volume=double(first[phase][0]);r.time=velocity->time();
        require(std::isfinite(r.volume)&&r.volume>=0,"invalid phase volume integral");
        const auto& map=forms->phaseVelocityIntegral[phase];require(map.rows()==2&&map.cols()>=q.size(),"phase velocity integral has a stale layout");
        const Eigen::Vector2d u=map.leftCols(q.size())*q;require(u.allFinite(),"nonfinite phase velocity integral");
        r.velocityIntegral={u.x(),u.y()};r.present=r.volume>0;
        if(r.present){localMean[phase]={double(first[phase][1]/first[phase][0]),double(first[phase][2]/first[phase][0])};
            r.coordinateMean=origin+localMean[phase];r.meanVelocity={u.x()/r.volume,u.y()/r.volume};
        }else{r.coordinateMean=r.meanVelocity={undefined,undefined};r.centralSecondMoment.fill(undefined);}
    }
    // A second pass about the measured centroid avoids subtracting two
    // large global raw moments, especially in a translated physical domain.
    std::array<std::array<long double,3>,2> second{};
    for(const auto& f:forms->cells){const auto& cell=mesh.cells(mesh.cellIndexToLabel[f.index]);const auto h=cell.size(),c=cell.centerCoord()-origin;
        for(int phase=0;phase<2;phase++)if(result[phase].present){const auto& m=phase==0?f.liquid:f.gas;const auto d=c-localMean[phase];
            second[phase][0]+=d.x*d.x*m[0]+2*d.x*h.x*m[1]+h.x*h.x*m[3];
            second[phase][1]+=d.x*d.y*m[0]+d.y*h.x*m[1]+d.x*h.y*m[2]+h.x*h.y*m[4];
            second[phase][2]+=d.y*d.y*m[0]+2*d.y*h.y*m[2]+h.y*h.y*m[5];}
    }
    for(int phase=0;phase<2;phase++)if(result[phase].present)for(int k=0;k<3;k++)result[phase].centralSecondMoment[k]=double(second[phase][k]);
    return result;
}
const HydroCompatibleMomentMaterial& HydroCompatibleMomentFlow::material() const{require(impl_&&impl_->phase,"moment owner has no moving material");return impl_->phase->material;}
HydroCompatibleMomentFlow::NormalPolynomial HydroCompatibleMomentFlow::normalPolynomial(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    int face,VelocityLevel level) const{
    require(bool(impl_),"moved moment owner");const auto& p=*impl_;p.validate(mesh,domain);const auto& s=*p.state;
    require(face>=0&&face<int(p.space.faces().size()),"invalid moment normal face");
    const Eigen::VectorXd* velocity=&s.velocity;double time=s.time;
    if(level==VelocityLevel::PreviousEndpoint){require(s.acceptedSteps>0,"initial moment state has no previous endpoint");velocity=&s.previousVelocity;time=s.previousTime;}
    else if(level==VelocityLevel::Midpoint){require(bool(s.stageTime),"initial moment state has no solved midpoint");velocity=&s.midpointVelocity;time=*s.stageTime;}
    else require(level==VelocityLevel::Endpoint,"invalid moment velocity layer");
    const auto& f=p.space.faces()[face];
    if(f.first<0)return {0,0,domain.face(mesh,face).center.y,1,0,time};
    return {(*velocity)[f.first],(*velocity)[f.first+1],f.tangentCenter,f.tangentScale,f.area,time};
}
const HydroCompatibleMomentFlow::State& HydroCompatibleMomentFlow::Candidate::state() const{require(bool(state_),"consumed moment candidate");return *state_;}
const InterfaceState& HydroCompatibleMomentFlow::Candidate::interfaceState() const{require(state_&&phase_&&phase_->interface,"moment candidate has no moving interface");return *phase_->interface;}
HydroCompatibleMomentFlow HydroCompatibleMomentFlow::create(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    HydroCompatibleMomentSpace space,HydroCompatibleMomentViscosity stress,double rho,Eigen::VectorXd velocity,double time,const std::vector<HydroBoundaryPatch>& patches){
    stress.validate(mesh,domain,space);
    require(bool(stress.cellViscosity()),"complete flow needs the physical viscosity: use buildCellwise or createFrozenMaterial, not six stress moments");
    auto material=HydroCompatibleMomentMaterial::buildSinglePhase(mesh,domain,space,rho,*stress.cellViscosity(),time);
    auto result=createWithMass(mesh,domain,std::move(space),std::move(stress),rho,material.kineticMass(),std::move(velocity),time,patches);
    result.impl_->initializeMaterial(mesh,domain,std::move(material));return result;
}
HydroCompatibleMomentFlow HydroCompatibleMomentFlow::createWithMass(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    HydroCompatibleMomentSpace space,HydroCompatibleMomentViscosity stress,double rho,Matrix mass,Eigen::VectorXd velocity,double time,const std::vector<HydroBoundaryPatch>& patches,bool projectOpenInitialVelocity){
    space.validate(mesh,domain);stress.validate(mesh,domain,space);
    require(std::isfinite(rho)&&rho>0&&std::isfinite(time),"invalid moment material or initial clock");
    require(velocity.size()==space.mass().rows()&&velocity.allFinite(),"invalid initial moment velocity");
    auto p=std::make_unique<Impl>();p->space=std::move(space);p->stress=std::move(stress);p->rho=rho;
    p->mass=std::move(mass);
    p->nb=int(p->space.mass().rows());p->nv=int(p->stress.matrix().rows());p->np=int(p->space.divergence().rows());
    p->boundary=HydroCompatibleMomentBoundary::build(mesh,domain,p->space,p->stress,patches);
    p->boundaryPatches=patches;
    const auto& fixed=p->boundary.fixed();
    for(int row=0;row<p->nb;row++)if(fixed[row])p->normalBoundary.push_back(row);
    p->fixedBoundary=fixed;
    if(p->boundary.open()&&projectOpenInitialVelocity){
        // Explicit constrained initialization, independent of viscosity/dt.
        // Pressure here is only a multiplier of the initial kinetic projection
        // and is not published as a solved physical pressure stage.
        const auto& d=p->space.divergence();const Eigen::VectorXd prescribed=p->boundary.values().head(p->nb);
        std::vector<int> rows(p->nb,-1);int count=0;for(int i=0;i<p->nb;i++)if(!fixed[i])rows[i]=count++;
        std::vector<Eigen::Triplet<double>> entries;
        for(int j=0;j<p->mass.outerSize();j++)for(Matrix::InnerIterator it(p->mass,j);it;++it)
            if(rows[it.row()]>=0&&rows[it.col()]>=0)entries.emplace_back(rows[it.row()],rows[it.col()],it.value());
        for(int j=0;j<d.outerSize();j++)for(Matrix::InnerIterator it(d,j);it;++it)if(rows[it.col()]>=0){
            entries.emplace_back(rows[it.col()],count+it.row(),-it.value());entries.emplace_back(count+it.row(),rows[it.col()],-it.value());}
        Matrix matrix(count+p->np,count+p->np);matrix.setFromTriplets(entries.begin(),entries.end());
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(matrix.rows()),weighted=p->mass*(velocity-prescribed);
        for(int i=0;i<p->nb;i++)if(rows[i]>=0)rhs[rows[i]]=weighted[i];rhs.tail(p->np)=d*prescribed;
        Eigen::SuperLU<Matrix> solve;solve.compute(matrix);require(solve.info()==Eigen::Success,"moment boundary initialization factorization failed");
        const Eigen::VectorXd solution=solve.solve(rhs);require(solve.info()==Eigen::Success&&solution.allFinite(),"moment boundary initialization failed");
        require((matrix*solution-rhs).norm()<=1e-11*std::max(rhs.norm(),std::numeric_limits<double>::min()),"moment boundary initialization residual rejected");
        Eigen::VectorXd initialized=prescribed;for(int i=0;i<p->nb;i++)if(rows[i]>=0)initialized[i]=solution[rows[i]];
        const Eigen::VectorXd change=initialized-velocity;p->initialCorrection=std::sqrt(std::max(0.,change.dot(p->mass*change)));velocity=std::move(initialized);
    }
    for(int row:p->normalBoundary)require(velocity[row]==p->boundary.values()[row],"initial moment velocity violates prescribed normal flux");
    require(divergenceResidual(p->space.divergence(),velocity)<=4096*std::numeric_limits<double>::epsilon(),"initial moment velocity is not incompressible; no hidden projection");
    p->configureGauge(mesh,domain);
    auto initial=std::make_shared<State>();initial->velocity=std::move(velocity);initial->previousVelocity=initial->velocity;initial->time=initial->previousTime=time;
    // No pressure/trace exists at an initial state with no solved stage.
    p->state=std::move(initial);return HydroCompatibleMomentFlow(std::move(p));
}
HydroCompatibleMomentFlow HydroCompatibleMomentFlow::createFrozenMaterial(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    HydroCompatibleMomentSpace space,HydroCompatibleMomentMaterial material,Eigen::VectorXd velocity,double time,const std::vector<HydroBoundaryPatch>& patches){
    material.validate(mesh,domain,space);require(sameTime(material.geometryTime(),time),"frozen material and initial geometry times differ");
    auto stress=HydroCompatibleMomentViscosity::build(mesh,domain,space,material.viscosityMoments());
    auto result=createWithMass(mesh,domain,std::move(space),std::move(stress),1.,material.kineticMass(),std::move(velocity),time,patches);
    result.impl_->initializeMaterial(mesh,domain,std::move(material));return result;
}
HydroCompatibleMomentFlow HydroCompatibleMomentFlow::createTwoPhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    HydroCompatibleMomentSpace space,InterfaceState interface,HydroCompatibleMomentMaterial::Properties properties,Eigen::VectorXd velocity,double sigma,const std::vector<HydroBoundaryPatch>& patches){
    InterfaceParametricMotion::validateDomain(mesh,domain,interface);
    require(std::isfinite(sigma)&&sigma>=0,"invalid moment surface tension");
    require(interface.progress().completedSteps==0&&interface.progress().endTime==0,"moving moment owner requires an initial interface");
    auto material=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,space,interface.parametricGeometry().boundaryHandle(),0,properties);
    auto stress=HydroCompatibleMomentViscosity::build(mesh,domain,space,material.viscosityMoments());
    for(const auto& patch:patches)if(patch.type!=HydroBoundaryType::FreeSlip)
        require(patch.incomingPhase==-1||patch.incomingPhase==1,"two-phase port needs an explicit exterior phase");
    auto result=createWithMass(mesh,domain,std::move(space),std::move(stress),1.,material.kineticMass(),std::move(velocity),0,patches);
    result.impl_->initializeMaterial(mesh,domain,std::move(material),std::make_shared<const InterfaceState>(std::move(interface)),sigma);return result;
}
HydroCompatibleMomentFlow::Load HydroCompatibleMomentFlow::makeLoad(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,Eigen::VectorXd force,double time) const{
    require(bool(impl_),"moved moment owner");impl_->validate(mesh,domain);
    require(force.size()==impl_->nb&&force.allFinite()&&std::isfinite(time),"invalid moment force load");
    Load result;result.owner_=impl_->owner;result.coefficients_=std::move(force);result.time_=time;return result;
}
HydroCompatibleMomentFlow::Candidate HydroCompatibleMomentFlow::prepareStokesStep(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroStepTime& time,const Load& load){return prepareStokesStep(mesh,domain,time,load,Options{});}
HydroCompatibleMomentFlow::Load HydroCompatibleMomentFlow::makeBodyForceLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::function<Vector2D(Vector2D)>& force,double time) const{
    require(bool(impl_)&&bool(force)&&std::isfinite(time),"invalid physical body force");const auto& p=*impl_;p.validate(mesh,domain);
    require(bool(p.state->physicalVelocity),"body force needs the complete owned velocity");
    require(bool(p.complement),"complete body force lacks its physical pressure tests");
    const MomentSplineBasis::GaussTwelve rule;std::vector<HydroCompatibleMomentSpline::PointForce> points;
    for(const auto& cell:p.complement->cells()){
        const auto& xs=cell.partition[0];const auto& ys=cell.partition[1];
        for(std::size_t x=1;x<xs.size();x++)for(std::size_t y=1;y<ys.size();y++)
        for(int half=0;half<2;half++)for(int i=0;i<12;i++)for(int j=0;j<12;j++){
            const double hx=xs[x]-xs[x-1],hy=ys[y]-ys[y-1];
            const Vector2D at{xs[x-1]+.5*hx*(half+.5+.5*rule.x[i]),ys[y-1]+.5*hy*(1+rule.x[j])};
            const double weight=.125*hx*hy*rule.w[i]*rule.w[j]*(domain.axisymmetric()?CartesianMeasure::twoPi*at.x:1.);
            const auto value=force(at);require(std::isfinite(value.x)&&std::isfinite(value.y),"nonfinite physical body force");
            points.push_back({at,value*weight});
        }
    }
    return makeQuadratureLoad(mesh,domain,points,time);
}
HydroCompatibleMomentFlow::Load HydroCompatibleMomentFlow::makeQuadratureLoad(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::vector<HydroCompatibleMomentSpline::PointForce>& points,double time) const{
    require(bool(impl_)&&std::isfinite(time),"invalid physical quadrature load");const auto& p=*impl_;p.validate(mesh,domain);
    require(p.state->physicalVelocity&&p.complement,"quadrature load needs complete velocity and pressure tests");
    for(const auto& point:points)require(std::isfinite(point.point.x)&&std::isfinite(point.point.y),"nonfinite quadrature location");
    auto result=makeLoad(mesh,domain,Eigen::VectorXd::Zero(p.nb),time);
    const auto q=p.state->physicalVelocity->transposePointForces(points),eta=p.complement->transposePointForces(points);
    result.physicalCoefficients_.resize(q.size()+eta.size());result.physicalCoefficients_<<q,eta;return result;
}
HydroCompatibleMomentFlow::Load HydroCompatibleMomentFlow::makeUniformAccelerationLoad(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,Vector2D acceleration,double time) const{
    require(bool(impl_)&&std::isfinite(acceleration.x)&&std::isfinite(acceleration.y),"invalid uniform moment acceleration");const auto& p=*impl_;
    auto load=makeLoad(mesh,domain,Eigen::VectorXd::Zero(p.nb),time);load.acceleration_=acceleration;return load;
}
HydroCompatibleMomentFlow::Candidate HydroCompatibleMomentFlow::prepareStokesStep(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroStepTime& time,const Load& load,const Options& options){return prepareStep(mesh,domain,time,load,options,false);}
HydroCompatibleMomentFlow::Candidate HydroCompatibleMomentFlow::prepareNavierStokesStep(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroStepTime& time,const Load& load){return prepareNavierStokesStep(mesh,domain,time,load,Options{});}
HydroCompatibleMomentFlow::Candidate HydroCompatibleMomentFlow::prepareNavierStokesStep(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroStepTime& time,const Load& load,const Options& options){return prepareStep(mesh,domain,time,load,options,true);}
HydroCompatibleMomentFlow::Candidate HydroCompatibleMomentFlow::prepareStep(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const HydroStepTime& time,const Load& load,const Options& options,bool convection){
    require(bool(impl_),"moved moment owner");auto& p=*impl_;p.validate(mesh,domain);const auto& old=*p.state;
    require(p.phase&&old.physicalVelocity&&p.complement,"flow has no complete physical velocity/material state");
    require(!convection||p.movingMaterial()||std::isfinite(p.rho),"variable-density convection requires a consistent phase mass flux; frozen material admits Stokes only");
    require(!p.movingMaterial()||convection,"moving material requires its coupled momentum transport");
    require(!convection||options.convectionFlux==HydroCompatibleMomentTransport::Flux::Central,"complete C1 momentum uses central transport; legacy traceResidual is not a complete-field model");
    require(std::isfinite(time.begin)&&std::isfinite(time.dt)&&time.dt>0&&std::isfinite(time.begin+time.dt)&&time.begin+time.dt>time.begin,"invalid moment step interval");
    require(time.index==old.acceptedSteps+1&&sameTime(time.begin,old.time),"moment step does not continue accepted clock");
    const bool fixed=time.policy==HydroTimePolicy::Fixed;
    require(fixed||time.policy==HydroTimePolicy::ConservativeOneStep,"invalid moment step policy");
    require(!p.fixedPolicy||*p.fixedPolicy==fixed,"moment time policy changed after acceptance");
    require(!fixed||old.acceptedSteps==0||time.dt==old.lastDt,"fixed moment time step changed");
    const double stage=time.begin+.5*time.dt;
    require(load.owner_==p.owner&&sameTime(load.time_,stage),"moment force source/time differs from midpoint");
    require(std::isfinite(options.relativeResidualTolerance)&&options.relativeResidualTolerance>0,"invalid moment residual tolerance");
    require(options.maximumIterations>0,"invalid moment iteration limit");
    require(options.accelerationDepth>=0&&options.accelerationDepth<=8,"invalid moment acceleration history");
    return prepareNativeStep(mesh,domain,time,load,options,convection);
}
HydroCompatibleMomentFlow::Candidate HydroCompatibleMomentFlow::prepareNativeStep(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroStepTime& time,const Load& load,const Options& options,bool convection){
    auto& p=*impl_;const auto& old=*p.state;const auto& field=*old.physicalVelocity;
    const auto& complement=*p.complement;const int nq=int(field.coordinates().size()),ne=int(complement.divergence().cols()),nv=nq+ne;
    const int np=p.np,nb=int(p.nativeBoundary.rows()),ng=p.boundary.hasPressureOutlet()?0:1,total=nv+np+nb+ng;
    const double dt=time.dt,stage=time.begin+.5*dt;
    const Matrix& before=p.phase->nativeMaterial->mass;
    Eigen::VectorXd initial=Eigen::VectorXd::Zero(nv);initial.head(nq)=field.coordinates();
    const auto lift=[&](const Eigen::VectorXd& force){Eigen::VectorXd result(nv);result.head(nq)=field.momentMap().transpose()*force;result.tail(ne)=complement.nativeMoments().transpose()*force;return result;};
    const Eigen::VectorXd external=load.physicalCoefficients_.size()?load.physicalCoefficients_:lift(load.coefficients_);
    require(external.size()==nv&&external.allFinite(),"complete body force changed the physical test space");
    // The two rows of momentum already integrate rho times every physical
    // test velocity, including the pressure complement. Reuse the same cut
    // geometry and positive quadrature at the actual nonlinear midpoint.
    const auto externalAt=[&](const PhaseData& phase){
        Eigen::VectorXd force=external;
        if(load.acceleration_){const auto g=*load.acceleration_;
            force+=phase.nativeMidpointMaterial->momentum.transpose()*Eigen::Vector2d(g.x,g.y);}
        require(force.allFinite(),"nonfinite phase-integrated acceleration load");return force;
    };
    const Eigen::VectorXd boundaryForce=lift(p.boundary.traction().head(p.nb));
    std::vector<Eigen::Triplet<double>> constraints;
    const auto& c=complement.pressureLoad();
    for(int col=0;col<c.outerSize();col++)for(Matrix::InnerIterator it(c,col);it;++it){constraints.emplace_back(nq+it.row(),nv+col,it.value());constraints.emplace_back(nv+col,nq+it.row(),it.value());}
    for(int col=0;col<p.nativeBoundary.outerSize();col++)for(Matrix::InnerIterator it(p.nativeBoundary,col);it;++it){
        constraints.emplace_back(col,nv+np+it.row(),it.value());constraints.emplace_back(nv+np+it.row(),col,it.value());}
    if(ng)for(int row=0;row<np;row++)if(p.gauge[row]){constraints.emplace_back(nv+row,total-1,p.gauge[row]);constraints.emplace_back(total-1,nv+row,p.gauge[row]);}
    Matrix constraint(total,total);constraint.setFromTriplets(constraints.begin(),constraints.end());
    const auto system=[&](const PhaseData& phase,const HydroCompatibleMomentSpline::NativeMomentumTransport& transport){
        Matrix a=(2/dt)*phase.nativeMaterial->mass+phase.nativeMidpointMaterial->viscosity+(1/dt)*transport.central;
        std::vector<Eigen::Triplet<double>> entries=constraints;
        for(int col=0;col<a.outerSize();col++)for(Matrix::InnerIterator it(a,col);it;++it)entries.emplace_back(it.row(),it.col(),it.value());
        Matrix result(total,total);result.setFromTriplets(entries.begin(),entries.end());return result;
    };
    const auto rhs=[&](const PhaseData& phase,const HydroCompatibleMomentSpline::NativeMomentumTransport& transport){
        Eigen::VectorXd value=Eigen::VectorXd::Zero(total);
        value.head(nv)=((before+phase.nativeMaterial->mass)*initial+transport.boundaryLoad)/dt+externalAt(phase)+boundaryForce;
        if(phase.traction)value.head(nv)+=phase.traction->variation;
        value.segment(nv+np,nb)=p.nativeBoundaryValues;return value;
    };
    const bool moving=p.movingMaterial();Options coupled=options;
    Matrix capillary(nq,nq),capillarySystem(total,total);
    if(moving&&p.surfaceTension>0){
        MomentWorkTimer timer("capillary_iteration",time.index);
        capillary=(.5*dt)*HydroCompatibleMomentTraction::iterationMetric(mesh,domain,p.space,field,p.phase->material,p.surfaceTension);
        std::vector<Eigen::Triplet<double>> entries;entries.reserve(capillary.nonZeros());
        for(int col=0;col<capillary.outerSize();col++)for(Matrix::InnerIterator it(capillary,col);it;++it)entries.emplace_back(it.row(),it.col(),it.value());
        capillarySystem.setFromTriplets(entries.begin(),entries.end());timer.mark("surface_metric");
    }
    if(moving){auto& knots=coupled.interfaceProjection.minimumMaterialKnots;
        for(int i=0;i<=coupled.interfaceProjection.materialSegments;i++)knots.push_back(double(i)/coupled.interfaceProjection.materialSegments);
        std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());}
    bool enriched=false;std::shared_ptr<const PhaseData> phase;
    HydroCompatibleMomentSpline::NativeMomentumTransport transport;double inventory=0,geometric=0;
    std::vector<HydroCompatibleMomentMassCurrent::FacePhase> faceVolumes(domain.nativeFaceCount());
    const InterfaceStepTime interval{time.index,time.begin,time.begin+dt,dt};
    if(!moving){
        auto fixed=std::make_shared<PhaseData>(*p.phase);fixed->previousMaterial=p.phase->material;
        fixed->midpointMaterial=p.phase->material.atTime(stage);fixed->material=p.phase->material.atTime(interval.end);
        fixed->nativePreviousMaterial=p.phase->nativeMaterial;fixed->nativeMidpointMaterial=fixed->nativeMaterial=p.phase->nativeMaterial;
        phase=std::move(fixed);
    }
    const auto evaluate=[&](const Eigen::VectorXd& q,Options& trialOptions){
        // Each trial is rebuilt from the accepted p.phase, never from the
        // preceding trial. Release its large forms before building new ones;
        // accepted history and the current nonlinear secants stay untouched.
        transport=HydroCompatibleMomentSpline::NativeMomentumTransport{};
        if(moving){phase.reset();phase=p.preparePhase(mesh,domain,time,q,trialOptions,enriched);}
        const auto velocity=field.withCoordinates(q,stage);
        if(convection){
            MomentWorkTimer timer("native_transport",time.index);
            const auto current=moving
                ?HydroCompatibleMomentMassCurrent::build(mesh,domain,p.space,velocity,*p.phase->interface,*phase->interface,interval,phase->midpointMaterial,2e-12,&complement)
                :HydroCompatibleMomentMassCurrent::buildConstantDensity(mesh,domain,p.space,velocity,interval,phase->midpointMaterial,2e-12,&complement);
            timer.mark("mass_current");
            inventory=current.inventoryClosure();geometric=0;for(const auto& cell:current.cells())geometric=std::max(geometric,cell.constantDefect);
            faceVolumes=current.facePhase();transport=velocity.nativeMomentumTransport(mesh,domain,p.space,current,p.boundary,&complement);
            timer.mark("momentum_transport");
        }else{
            transport.central=Matrix(nv,nv);transport.boundaryMetric=Matrix(nv,nv);transport.geometricMetric=Matrix(nv,nv);
            transport.boundaryLoad=Eigen::VectorXd::Zero(nv);
            if(p.boundary.open()){
                const auto trace=velocity.normalFlux();const auto volumes=InterfaceParametricMotion::backgroundFaceVolumes(mesh,domain,trace,dt);
                for(std::size_t i=0;i<volumes.size();i++){const auto& v=volumes[i];faceVolumes[i]={v.total,v.totalPositive,v.totalNegative,0,0};}
            }
        }
    };
    Candidate result;result.owner_=p.owner;result.revision_=p.revision;result.fixedPolicy_=time.policy==HydroTimePolicy::Fixed;auto& report=result.report_;
    Eigen::VectorXd iterate=Eigen::VectorXd::Zero(total);iterate.head(nv)=initial;
    std::vector<Eigen::VectorXd> changes,images;std::vector<double> residuals;
    // Extrapolate the COMPLETE accepted velocity to the new midpoint, only
    // as the initial nonlinear guess. The endpoint on the physical RHS and
    // the final midpoint equations remain unchanged. In particular, an
    // instantaneous neck projection is not a time derivative.
    Eigen::VectorXd seed=field.coordinates();bool predicted=old.continuousHistory;
    if(predicted){
        require(old.previousPhysicalVelocity&&old.time>old.previousTime
            &&old.previousPhysicalVelocity->coordinates().size()==nq,"incomplete continuous velocity history");
        seed+=(.5*dt/(old.time-old.previousTime))*(field.coordinates()-old.previousPhysicalVelocity->coordinates());
        predicted=seed.allFinite();
    }
    if(predicted){
        Options trialOptions=coupled;
        const auto rejectPrediction=[&](const std::exception& error){
            // A predictor has no accepted physical status. If its geometry
            // cannot be prepared, start the SAME solve at the old endpoint;
            // any failure there still propagates. Failed knots are discarded.
            predicted=false;enriched=false;
            if(std::getenv("QUANTUM_MOMENT_STEP_PROFILE"))std::cout<<"[native_midpoint_predictor_rejected] step="<<time.index<<" reason="<<error.what()<<std::endl;
        };
        try{evaluate(seed,trialOptions);coupled=std::move(trialOptions);}
        catch(const std::invalid_argument& error){rejectPrediction(error);}
        catch(const std::runtime_error& error){rejectPrediction(error);}
    }
    if(!predicted){seed=field.coordinates();evaluate(seed,coupled);}
    iterate.head(nq)=seed;Eigen::VectorXd solution;
    const auto fullResidual=[&](const Eigen::VectorXd& value){
        const Matrix actual=system(*phase,transport);const Eigen::VectorXd actualLoad=rhs(*phase,transport);
        double scale=actualLoad.norm();if(phase->traction)scale+=(c*phase->traction->pressurePotential).norm();
        return (actual*value-actualLoad).norm()/std::max(scale,std::numeric_limits<double>::min());
    };
    double currentResidual=fullResidual(iterate);
    require(std::isfinite(currentResidual),"nonfinite initial complete midpoint residual");
    if(std::getenv("QUANTUM_MOMENT_STEP_PROFILE"))std::cout<<"[native_midpoint_initial] step="<<time.index
        <<" history="<<predicted<<" residual="<<currentResidual<<std::endl;
    if(convection)p.factorDt.reset();
    report.reusedFactorization=true;
    for(int iteration=0;iteration<options.maximumIterations;iteration++){
        // (A(q_k)+S) q_* = b(q_k)+S q_k. The surface-gradient metric S
        // anticipates midpoint capillary response but cancels at the fixed
        // point. fullResidual still assembles the ORIGINAL physical system.
        {
        Matrix matrix=system(*phase,transport);Eigen::VectorXd loadVector=rhs(*phase,transport);
        if(capillary.nonZeros()){matrix+=capillarySystem;loadVector.head(nq)+=capillary*iterate.head(nq);}
        double scale=loadVector.norm();
        if(phase->traction)scale+=(c*phase->traction->pressurePotential).norm();
        scale=std::max(scale,std::numeric_limits<double>::min());
        MomentWorkTimer timer("native_linear",time.index);auto& solver=p.solver;
        const int factorsBefore=solver.statistics().velocityFactorizations;
        const bool reused=!convection&&!capillary.nonZeros()&&p.factorDt&&*p.factorDt==dt;
        if(!reused){
            timer.mark("block_factor_begin");if(!convection)p.factorDt.reset();
            solver.compute(matrix,nq,ne,np,nb,ng,field.cellInteriorCoordinates());if(!convection&&!capillary.nonZeros())p.factorDt=dt;
        }
        timer.mark(reused?"reuse_factor":"factor");
        solution=solver.solve(loadVector,convection?.01*options.relativeResidualTolerance*scale:0);timer.mark("solve");
        const auto& linear=solver.statistics();report.reusedFactorization&=linear.velocityFactorizations==factorsBefore;
        if(timer.enabled)std::cout<<"[native_block] step="<<time.index<<" iteration="<<iteration+1
            <<" velocity_factors="<<linear.velocityFactorizations<<" pressure_factors="<<linear.pressureFactorizations
            <<" krylov="<<linear.krylovIterations<<std::endl;
        require((matrix*solution-loadVector).norm()/scale<=options.relativeResidualTolerance,"complete midpoint linear residual rejected");
        // eta=0 is imposed by exact block elimination before factorization.
        // Pressure is recovered from the remaining momentum equations. This
        // linear check includes S; physical acceptance below excludes it.
        report.complementResidual=solution.segment(nq,ne).norm();
        } // The assembled linear system is not needed during trial geometry.
        const Eigen::VectorXd picard=solution;bool accelerated=false;
        if(options.accelerationDepth>0){
            changes.push_back(solution-iterate);images.push_back(solution);
            if(int(images.size())>options.accelerationDepth+1){images.erase(images.begin());changes.erase(changes.begin());}
            const int count=int(images.size())-1;
            if(count>0){Eigen::MatrixXd difference(total,count),mapped(total,count);
                for(int k=0;k<count;k++){difference.col(k)=changes[k+1]-changes[k];mapped.col(k)=images[k+1]-images[k];}
                const auto metric=[&](const Eigen::MatrixXd& values){
                    Eigen::MatrixXd weighted=Eigen::MatrixXd::Zero(total,values.cols());
                    // Use dt/2 times the positive velocity part of the SAME
                    // current Picard matrix: M^{n+1} + dt*K_mu/2 + dt*S/2.
                    // A mass-only norm would underweight short surface waves
                    // whose corrections the capillary preconditioner reduces.
                    // Advection is excluded from this positive metric; eta is
                    // zero and pressure/multipliers remain algebraic responses.
                    weighted.topRows(nv)=phase->nativeMaterial->mass*values.topRows(nv)
                        +(.5*dt)*(phase->nativeMidpointMaterial->viscosity*values.topRows(nv));
                    if(capillary.nonZeros())weighted.topRows(nq).noalias()+=(.5*dt)*(capillary*values.topRows(nq));
                    return weighted;
                };
                Eigen::MatrixXd gram=difference.transpose()*metric(difference);Eigen::VectorXd scales(count);
                for(int k=0;k<count;k++)scales[k]=gram(k,k)>0?1/std::sqrt(gram(k,k)):0;
                gram=scales.asDiagonal()*gram*scales.asDiagonal();gram=.5*(gram+gram.transpose()).eval();
                const Eigen::VectorXd right=scales.asDiagonal()*difference.transpose()*metric(changes.back());
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(gram);require(eigen.info()==Eigen::Success,"complete midpoint acceleration failed");
                Eigen::VectorXd coordinates=eigen.eigenvectors().transpose()*right;const double rank=4096*std::numeric_limits<double>::epsilon()*eigen.eigenvalues().maxCoeff();
                for(int k=0;k<count;k++)coordinates[k]=eigen.eigenvalues()[k]>rank?coordinates[k]/eigen.eigenvalues()[k]:0;
                solution-=mapped*(scales.asDiagonal()*eigen.eigenvectors()*coordinates);
                require(solution.allFinite(),"nonfinite complete midpoint acceleration");
                accelerated=true;
            }
        }
        // Histories live in the fixed complete [q,eta,p,lambda] space. Material
        // enrichment changes the approximate nonlinear map, not these physical
        // coordinates; use old secants only as a proposal, guarded by the
        // newly assembled full residual. Rejected trials never commit knots.
        const Eigen::VectorXd proposal=solution;bool accepted=false;double acceptedScale=0;
        int trials=0;std::string resolutionFailure;
        for(int direction=0;direction<(accelerated?2:1)&&!accepted;direction++){
            if(direction){changes.clear();images.clear();changes.push_back(picard-iterate);images.push_back(picard);}
            const Eigen::VectorXd update=(direction?picard:proposal)-iterate;
            double stepScale=1;
            for(int backtrack=0;backtrack<16;backtrack++,stepScale*=.5){
                solution=iterate+stepScale*update;Options trialOptions=coupled;
                enriched=false;++trials;double trialResidual=std::numeric_limits<double>::infinity();
                try{
                    // The entire material, traction and flux problem belongs
                    // to this trial q. A previous frozen-wind residual is not
                    // an acceptance test, even for a damped Picard update.
                    evaluate(solution.head(nq),trialOptions);trialResidual=fullResidual(solution);
                }catch(const InterfaceParametricProjection::MaterialResolutionLimit& error){
                    resolutionFailure=error.what();
                }
                if(std::getenv("QUANTUM_MOMENT_STEP_PROFILE"))std::cout<<"[native_midpoint_trial] step="<<time.index
                    <<" iteration="<<iteration+1<<" direction="<<direction<<" scale="<<stepScale
                    <<" residual="<<trialResidual<<" previous="<<currentResidual<<std::endl;
                if(std::isfinite(trialResidual)&&(trialResidual<=options.relativeResidualTolerance||
                    trialResidual<=(1-1e-4*stepScale)*currentResidual)){
                    coupled=std::move(trialOptions);report.equationResidual=trialResidual;
                    accepted=true;acceptedScale=stepScale;break;
                }
            }
        }
        if(!accepted){std::ostringstream reason;reason<<"complete midpoint nonlinear update rejected; iteration="<<iteration+1
            <<" residual="<<currentResidual<<" trials="<<trials;
            if(!resolutionFailure.empty())reason<<"; last geometry rejection: "<<resolutionFailure;
            throw std::invalid_argument(reason.str());}
        currentResidual=report.equationResidual;
        report.nonlinearIterations=iteration+1;residuals.push_back(report.equationResidual);iterate=solution;
        if(std::getenv("QUANTUM_MOMENT_STEP_PROFILE"))std::cout<<"[native_midpoint] step="<<time.index<<" iteration="<<iteration+1
            <<" residual="<<report.equationResidual<<" update_scale="<<acceptedScale<<" trials="<<trials<<std::endl;
        if(!enriched&&report.equationResidual<=options.relativeResidualTolerance)break;
    }
    if(!(report.equationResidual<=options.relativeResidualTolerance)){
        std::ostringstream reason;reason<<"complete midpoint full equation rejected; residuals:";for(double r:residuals)reason<<' '<<r;throw std::invalid_argument(reason.str());
    }
    const Eigen::VectorXd middle=solution.head(nv),endpoint=2*middle-initial;
    auto next=std::make_shared<State>();next->physicalVelocity=field.withCoordinates(endpoint.head(nq),time.begin+dt);
    next->previousPhysicalVelocity=old.physicalVelocity;next->midpointPhysicalVelocity=field.withCoordinates(middle.head(nq),stage);
    next->velocity=next->physicalVelocity->nativeMoments();next->previousVelocity=old.velocity;next->midpointVelocity=next->midpointPhysicalVelocity->nativeMoments();
    next->pressure=solution.segment(nv,np);
    if(phase->traction)next->pressure+=phase->traction->pressurePotential;
    next->acceptedSteps=time.index;next->time=time.begin+dt;next->previousTime=old.time;next->lastDt=dt;next->stageTime=stage;
    next->continuousHistory=true;
    report.boundaryConstraintResidual=(p.allNativeBoundary*endpoint.head(nq)-p.allNativeBoundaryValues).norm();
    report.divergenceResidual=divergenceResidual(p.space.divergence(),next->velocity);
    require(report.boundaryConstraintResidual<=2e-9&&report.divergenceResidual<=options.relativeResidualTolerance,"complete midpoint boundary or divergence rejected");
    const Matrix& after=phase->nativeMaterial->mass;const Matrix delta=after-before;
    report.kineticBefore=.5*initial.dot(before*initial);report.kineticAfter=.5*endpoint.dot(after*endpoint);
    report.viscousDissipation=dt*middle.dot(phase->nativeMidpointMaterial->viscosity*middle);
    report.externalWork=dt*middle.dot(externalAt(*phase));report.boundaryTractionWork=dt*middle.dot(boundaryForce);
    report.boundaryReactionWork=-dt*solution.segment(nv+np,nb).dot(p.nativeBoundaryValues);
    report.pressureWork=-dt*middle.tail(ne).dot(c*next->pressure);
    report.centralConvectionWork=middle.dot(transport.central*middle-transport.boundaryLoad);
    report.temporalMassWork=.125*(endpoint-initial).dot(delta*(endpoint-initial));
    report.materialMetricWork=report.centralConvectionWork+.5*middle.dot(delta*middle);
    report.geometricMassWork=.5*middle.dot(transport.geometricMetric*middle);
    report.boundaryKineticFlux=transport.boundaryKineticFlux;report.boundaryInflowDissipation=transport.boundaryInflowDissipation;
    report.phaseMassClosure=inventory;report.geometricMassClosure=geometric;
    report.massMetricDefect=(transport.central+Matrix(transport.central.transpose())+delta-transport.boundaryMetric-transport.geometricMetric).norm()/std::max(before.norm(),std::numeric_limits<double>::min());
    for(const auto face:domain.faces(mesh))if(face.outwardSign&&face.area>0){const auto& volume=faceVolumes.at(face.nativeRow);
        report.boundaryVolumeRate+=face.outwardSign*volume.total/dt;
        report.boundaryInflowRate+=(face.outwardSign>0?volume.totalNegative:volume.totalPositive)/dt;
        report.boundaryOutflowRate+=(face.outwardSign>0?volume.totalPositive:volume.totalNegative)/dt;}
    if(phase->traction){report.capillaryWork=dt*middle.dot(phase->traction->variation)+dt*middle.tail(ne).dot(c*phase->traction->pressurePotential);
        report.surfaceEnergyBefore=p.surfaceTension*p.phase->length;report.surfaceEnergyAfter=p.surfaceTension*phase->length;
        const auto& surface=phase->interface->parametricFluxes().surface;
        report.surfaceAdvection=p.surfaceTension*surface.advection;report.surfaceEndpointWork=p.surfaceTension*surface.endpointWork;
        report.surfaceFluxError=p.surfaceTension*surface.quadratureError;report.surfaceFluxBudget=p.surfaceTension*surface.budget;report.surfaceCriticalError=p.surfaceTension*surface.criticalError;
        report.surfaceWorkDefect=report.surfaceEnergyAfter-report.surfaceEnergyBefore+report.capillaryWork+report.surfaceAdvection-report.surfaceEndpointWork;}
    report.energyDefect=report.kineticAfter-report.kineticBefore+report.viscousDissipation+report.materialMetricWork-report.temporalMassWork
        -report.externalWork-report.pressureWork-report.boundaryTractionWork-report.boundaryReactionWork-report.capillaryWork;
    require(std::isfinite(report.energyDefect),"nonfinite complete midpoint energy");result.state_=std::move(next);result.phase_=std::move(phase);return result;
}
void HydroCompatibleMomentFlow::accept(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,Candidate&& candidate){
    require(bool(impl_),"moved moment owner");auto& p=*impl_;p.validate(mesh,domain);
    require(candidate.owner_==p.owner&&candidate.state_&&candidate.revision_==p.revision,"moment candidate is consumed, foreign or stale");
    p.state=std::move(candidate.state_);p.phase=std::move(candidate.phase_);p.fixedPolicy=candidate.fixedPolicy_;++p.revision;candidate.owner_.reset();
}
const HydroCompatibleMomentFlow::State& HydroCompatibleMomentFlow::NeckCandidate::state() const{
    require(bool(state_),"consumed moment neck candidate");return *state_;
}
const InterfaceState& HydroCompatibleMomentFlow::NeckCandidate::interfaceState() const{
    require(state_&&phase_&&phase_->interface,"consumed moment neck candidate");return *phase_->interface;
}
HydroCompatibleMomentFlow::NeckCandidate HydroCompatibleMomentFlow::prepareNeckSplit(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceParametricProjection::NeckSplitOptions& options) const{
    auto assessment=assessNeckSplit(mesh,domain,options);require(bool(assessment.candidate),"neck event has not entered the coupled feasible set");return std::move(*assessment.candidate);
}
HydroCompatibleMomentFlow::NeckAssessment HydroCompatibleMomentFlow::assessNeckSplit(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceParametricProjection::NeckSplitOptions& options) const{
    require(impl_&&impl_->movingMaterial(),"moment neck event requires moving two-phase ownership");const auto& p=*impl_;p.validate(mesh,domain);
    auto geometry=InterfaceParametricMotion::assessNeckSplit(mesh,domain,*p.phase->interface,options);
    NeckAssessment assessment;assessment.report.geometry=geometry.report;
    if(!geometry.state){assessment.readiness=geometry.readiness==InterfaceParametricProjection::NeckReadiness::ScalePending?NeckReadiness::ScalePending:NeckReadiness::CapPending;return assessment;}
    auto interface=std::move(*geometry.state);
    auto phase=std::make_shared<PhaseData>(*p.phase);auto state=std::make_shared<State>(*p.state);
    phase->material=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,p.space,interface.parametricGeometry().boundaryHandle(),state->time,p.phase->material.properties());
    phase->length=HydroCompatibleMomentTraction::measure(*phase->material.boundaryHandle(),mesh,domain).length;
    NeckCandidate result;result.owner_=p.owner;result.revision_=p.revision;auto& report=assessment.report;
    report.geometry=interface.parametricNeckEvents().back().report;
    report.surfaceBefore=p.surfaceTension*HydroCompatibleMomentTraction::measure(*p.phase->material.boundaryHandle(),mesh,domain).length;
    report.surfaceAfter=p.surfaceTension*phase->length;
    MomentTransfer::Feasibility feasible;
    require(bool(state->physicalVelocity),"neck event needs the complete accepted velocity");
    {
        auto projection=state->physicalVelocity->projectionForms(mesh,domain,p.space,phase->material,*state->physicalVelocity,&*p.complement);
        const auto& cells=InterfaceParametricMotion::cellIndices(interface);InterfaceParametricProjection::MeshTransferReport overlay;
        for(int k=0;k<int(cells.size());k++)overlay.overlaps.push_back({k,k,{},0});
        auto transferred=MomentTransfer::nativeVelocity(*state->physicalVelocity,*p.phase->nativeMaterial,*state->physicalVelocity,projection,
            p.space,p.boundary,overlay,cells,cells,report.surfaceAfter-report.surfaceBefore,report.momentum,feasible);
        if(transferred){state->physicalVelocity=std::move(*transferred);state->velocity=state->physicalVelocity->nativeMoments();
            phase->nativeMaterial=evaluatedNativeForms(std::move(projection.material),*state->physicalVelocity);}
    }
    if(feasible!=MomentTransfer::Feasibility::Admissible){assessment.readiness=feasible==MomentTransfer::Feasibility::Momentum?NeckReadiness::MomentumPending:NeckReadiness::EnergyPending;return assessment;}
    report.kineticChange=report.momentum.kineticAfter-report.momentum.kineticBefore;
    report.totalEnergyChange=report.kineticChange+report.surfaceAfter-report.surfaceBefore;report.modelDissipation=-report.totalEnergyChange;
    report.relativeMassMatrixChange=(phase->nativeMaterial->mass-p.phase->nativeMaterial->mass).norm()/p.phase->nativeMaterial->mass.norm();
    if(!state->stageTime){
        state->previousVelocity=state->velocity;phase->previousMaterial=phase->midpointMaterial=phase->material;
        if(state->physicalVelocity){state->previousPhysicalVelocity=state->physicalVelocity;phase->nativePreviousMaterial=phase->nativeMidpointMaterial=phase->nativeMaterial;}
        phase->stress=HydroCompatibleMomentViscosity::build(mesh,domain,p.space,phase->material.viscosityMoments());
        phase->activeTangents.clear();
    }
    // Previous endpoint, pressure, viscous trace and midpoint material belong
    // to the already solved physical step; never rewrite that history as if
    // the instantaneous event had happened at its midpoint.
    state->continuousHistory=false;
    phase->interface=std::make_shared<const InterfaceState>(std::move(interface));
    result.report_=report;result.state_=std::move(state);result.phase_=std::move(phase);
    assessment.readiness=NeckReadiness::Ready;assessment.candidate=std::move(result);return assessment;
}
void HydroCompatibleMomentFlow::acceptNeckSplit(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,NeckCandidate&& candidate){
    require(impl_&&impl_->movingMaterial(),"moment neck event has no accepted two-phase owner");auto& p=*impl_;p.validate(mesh,domain);
    require(candidate.owner_==p.owner&&candidate.state_&&candidate.phase_&&candidate.revision_==p.revision,"moment neck candidate is consumed, foreign or stale");
    InterfaceParametricMotion::validateDomain(mesh,domain,*candidate.phase_->interface);
    require(candidate.state_->time==p.state->time&&candidate.state_->acceptedSteps==p.state->acceptedSteps,"neck event changed physical time");
    p.state=std::move(candidate.state_);p.phase=std::move(candidate.phase_);p.factorDt.reset();++p.revision;candidate.owner_.reset();
}
HydroCompatibleMomentFlow::MeshTransferCandidate::MeshTransferCandidate()=default;
HydroCompatibleMomentFlow::MeshTransferCandidate::~MeshTransferCandidate()=default;
HydroCompatibleMomentFlow::MeshTransferCandidate::MeshTransferCandidate(MeshTransferCandidate&&) noexcept=default;
HydroCompatibleMomentFlow::MeshTransferCandidate& HydroCompatibleMomentFlow::MeshTransferCandidate::operator=(MeshTransferCandidate&&) noexcept=default;
HydroCompatibleMomentFlow::MeshTransferCandidate HydroCompatibleMomentFlow::prepareMeshTransfer(
    const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Options& options) const{
    require(impl_&&impl_->movingMaterial(),"joint moment mesh transfer requires moving two-phase ownership");
    require(std::isfinite(options.relativeResidualTolerance)&&options.relativeResidualTolerance>0,"invalid mesh transfer residual tolerance");
    impl_->validate(oldMesh,oldDomain);
    auto assessment=prepareNativeMeshTransfer(oldMesh,oldDomain,mesh,domain,options,std::numeric_limits<double>::infinity());
    if(!assessment.candidate){const auto& r=assessment.report.levels.at(assessment.rejectedLevel);std::ostringstream message;message.precision(17);
        message<<"complete mesh transfer infeasible: level="<<assessment.rejectedLevel<<" constraint="<<r.constraintResidual<<" kinetic_before="<<r.kineticBefore<<" minimum="<<r.minimumKinetic;
        throw std::invalid_argument(message.str());}
    return std::move(*assessment.candidate);
}
HydroCompatibleMomentFlow::MeshTransferAssessment HydroCompatibleMomentFlow::assessMeshTransfer(
    const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Options& options,double relativeVelocityError) const{
    require(impl_&&impl_->movingMaterial()&&impl_->state->physicalVelocity,"optional mesh assessment requires the complete two-phase physical velocity");
    impl_->validate(oldMesh,oldDomain);
    require(std::isfinite(relativeVelocityError)&&relativeVelocityError>0&&relativeVelocityError<1,"invalid coarsening velocity error budget");
    require(std::isfinite(options.relativeResidualTolerance)&&options.relativeResidualTolerance>0,"invalid mesh transfer residual tolerance");
    return prepareNativeMeshTransfer(oldMesh,oldDomain,mesh,domain,options,relativeVelocityError);
}
HydroCompatibleMomentFlow::MeshTransferAssessment HydroCompatibleMomentFlow::prepareNativeMeshTransfer(
    const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Options& options,double relativeVelocityError) const{
    MeshTransferAssessment assessment;
    const auto& old=*impl_;const auto& state=*old.state;
    auto space=HydroCompatibleMomentSpace::build(mesh,domain);const auto properties=old.phase->material.properties();
    auto phase=std::make_shared<PhaseData>();
    phase->material=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,space,old.phase->material.boundaryHandle(),state.time,properties);
    phase->previousMaterial=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,space,old.phase->previousMaterial.boundaryHandle(),state.previousTime,properties);
    phase->midpointMaterial=HydroCompatibleMomentMaterial::buildFrozen(mesh,domain,space,old.phase->midpointMaterial.boundaryHandle(),state.stageTime.value_or(state.time),properties);
    phase->stress=HydroCompatibleMomentViscosity::build(mesh,domain,space,phase->midpointMaterial.viscosityMoments());
    auto boundary=HydroCompatibleMomentBoundary::build(mesh,domain,space,phase->stress,old.boundaryPatches).inheritNormalConstraints(mesh,domain,space,old.boundary);
    auto transferred=std::make_shared<State>(state);
    auto interface=InterfaceParametricMotion::prepareMeshTransfer(oldMesh,oldDomain,*old.phase->interface,mesh,domain,options.interfaceProjection);
    const auto& sourceCells=InterfaceParametricMotion::cellIndices(*old.phase->interface);const auto& targetCells=InterfaceParametricMotion::cellIndices(interface.state);
    MeshTransferCandidate result;result.sourceOwner_=old.owner;result.sourceRevision_=old.revision;result.report_.interface=interface.report;
    result.report_.levels.resize(state.stageTime?3:1);
    std::vector<int> parent(mesh.cellIndexNum(),-1);bool exactRestriction=true;
    for(const auto& overlap:interface.report.overlaps){const int from=sourceCells.at(overlap.source),to=targetCells.at(overlap.target);
        if(parent[to]>=0&&parent[to]!=from)exactRestriction=false;parent[to]=from;}
    require(std::all_of(parent.begin(),parent.end(),[](int c){return c>=0;}),"complete remap lost an overlap patch");
    // An embedded field is the exact zero-distance minimizer; retain its
    // algebraic restriction. Otherwise solve the same physical projection
    // on the full target space, not on its lower-order native observations.
    const auto target=exactRestriction?state.physicalVelocity->refine(mesh,domain,space,&boundary)
        :HydroCompatibleMomentSpline::initializeNative(mesh,domain,space,Eigen::VectorXd::Zero(space.mass().rows()),state.time,&boundary);
    auto complement=HydroCompatibleMomentComplement::build(mesh,domain,space,boundary.fixed());
    using Forms=HydroCompatibleMomentSpline::NativeMaterialForms;
    const auto remap=[&](const auto& source,const Forms& before,const auto& material,int level,std::shared_ptr<const Forms>& forms)
        -> std::optional<HydroCompatibleMomentSpline::NativeVelocity>{
        if(exactRestriction){auto value=target.restrictFrom(source);forms=std::make_shared<const Forms>(value.mixedMaterialForms(mesh,domain,space,material,complement));return value;}
        auto layout=target.withCoordinates(target.coordinates(),source.time());
        auto projection=layout.projectionForms(mesh,domain,space,material,source,&complement);MomentTransfer::Feasibility feasible;
        auto value=MomentTransfer::nativeVelocity(source,before,layout,projection,space,boundary,interface.report,sourceCells,targetCells,0,result.report_.levels[level],feasible);
        if(!value){assessment.readiness=feasible==MomentTransfer::Feasibility::Momentum?MeshTransferReadiness::MomentumPending:MeshTransferReadiness::EnergyPending;
            assessment.rejectedLevel=level;assessment.report=result.report_;return std::nullopt;}
        if(result.report_.levels[level].relativeProjectionError>relativeVelocityError){
            assessment.readiness=MeshTransferReadiness::AccuracyPending;assessment.rejectedLevel=level;assessment.report=result.report_;return std::nullopt;}
        forms=evaluatedNativeForms(std::move(projection.material),*value);return std::move(*value);
    };
    transferred->physicalVelocity=remap(*state.physicalVelocity,*old.phase->nativeMaterial,phase->material,0,phase->nativeMaterial);
    if(!transferred->physicalVelocity)return assessment;
    transferred->velocity=transferred->physicalVelocity->nativeMoments();
    if(state.stageTime){
        transferred->previousPhysicalVelocity=remap(*state.previousPhysicalVelocity,*old.phase->nativePreviousMaterial,phase->previousMaterial,1,phase->nativePreviousMaterial);
        if(!transferred->previousPhysicalVelocity)return assessment;
        transferred->midpointPhysicalVelocity=remap(*state.midpointPhysicalVelocity,*old.phase->nativeMidpointMaterial,phase->midpointMaterial,2,phase->nativeMidpointMaterial);
        if(!transferred->midpointPhysicalVelocity)return assessment;
        transferred->midpointVelocity=transferred->midpointPhysicalVelocity->nativeMoments();transferred->remappedStage=true;
    }else{transferred->previousPhysicalVelocity=transferred->physicalVelocity;phase->nativePreviousMaterial=phase->nativeMidpointMaterial=phase->nativeMaterial;}
    transferred->previousVelocity=transferred->previousPhysicalVelocity->nativeMoments();
    if(state.stageTime)transferred->pressure=MomentTransfer::pressure(oldMesh,state.pressure,mesh,domain,interface.report,sourceCells,targetCells);
    auto prepared=std::make_unique<Impl>();auto& p=*prepared;p.space=std::move(space);p.stress=phase->stress;p.boundary=std::move(boundary);p.boundaryPatches=old.boundaryPatches;
    p.nb=int(p.space.mass().rows());p.nv=int(p.stress.matrix().rows());p.np=int(p.space.divergence().rows());p.fixedBoundary=p.boundary.fixed();
    p.state=transferred;p.configureGauge(mesh,domain);p.configureNativeBoundary();
    p.complement=std::move(complement);
    const auto balance=[&](const auto& beforeField,const auto& afterField,const auto& beforeForms,const auto& afterForms,const auto& material,int level){
        if(!exactRestriction)return; // projected levels already checked the same owned forms
        auto& r=result.report_.levels[level];r.time=beforeField.time();
        const auto a=beforeField.cellBalance(beforeForms),b=afterField.cellBalance(afterForms);
        Eigen::MatrixXd sums=Eigen::MatrixXd::Zero(2,oldMesh.cellIndexNum());Eigen::VectorXd energy=Eigen::VectorXd::Zero(oldMesh.cellIndexNum());
        for(int cell=0;cell<mesh.cellIndexNum();cell++){sums.col(parent[cell])+=b.momentum.col(cell);energy[parent[cell]]+=b.kinetic[cell];}
        for(int cell=0;cell<oldMesh.cellIndexNum();cell++){
            const double scale=std::sqrt(std::max(0.,2*material.cells()[cell].density(0,0)*std::max(a.kinetic[cell],energy[cell])));
            const double difference=(sums.col(cell)-a.momentum.col(cell)).norm();
            r.momentumResidual=std::max(r.momentumResidual,scale>0?difference/scale:(difference==0?0:std::numeric_limits<double>::infinity()));
        }
        r.kineticBefore=a.kinetic.sum();r.kineticAfter=b.kinetic.sum();r.minimumKinetic=std::numeric_limits<double>::quiet_NaN();
        r.unconstrainedKinetic=r.kineticAfter;r.kineticCeiling=r.kineticBefore;r.exactFieldRestriction=true;
        r.patches=oldMesh.cellIndexNum();r.independentConstraints=int(p.nativeBoundary.rows());r.constraintResidual=(p.allNativeBoundary*afterField.coordinates()-p.allNativeBoundaryValues).norm();
        r.divergenceResidual=divergenceResidual(p.space.divergence(),afterField.nativeMoments());
        const double scale=std::max({r.kineticBefore,r.kineticAfter,std::numeric_limits<double>::min()});
        require(r.momentumResidual<=3.2e-11&&std::abs(r.kineticAfter-r.kineticBefore)<=3.2e-11*scale
            &&r.constraintResidual<=2e-9&&r.divergenceResidual<=options.relativeResidualTolerance,"complete refinement lost physical momentum, energy or original boundary constraints");
    };
    balance(*state.physicalVelocity,*transferred->physicalVelocity,*old.phase->nativeMaterial,*phase->nativeMaterial,old.phase->material,0);
    if(state.stageTime){
        balance(*state.previousPhysicalVelocity,*transferred->previousPhysicalVelocity,*old.phase->nativePreviousMaterial,*phase->nativePreviousMaterial,old.phase->previousMaterial,1);
        balance(*state.midpointPhysicalVelocity,*transferred->midpointPhysicalVelocity,*old.phase->nativeMidpointMaterial,*phase->nativeMidpointMaterial,old.phase->midpointMaterial,2);
    }
    phase->interface=std::make_shared<const InterfaceState>(std::move(interface.state));phase->length=old.phase->length;
    p.initialCorrection=old.initialCorrection;p.mass=phase->material.kineticMass();p.rho=phase->material.uniformDensity().value_or(std::numeric_limits<double>::quiet_NaN());
    p.surfaceTension=old.surfaceTension;p.fixedPolicy=old.fixedPolicy;p.phase=std::move(phase);result.prepared_=std::move(prepared);
    assessment.report=result.report_;assessment.readiness=MeshTransferReadiness::Ready;assessment.candidate=std::move(result);return assessment;
}
void HydroCompatibleMomentFlow::acceptMeshTransfer(const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,MeshTransferCandidate&& candidate){
    require(bool(impl_),"moved moment owner");impl_->validate(oldMesh,oldDomain);
    require(candidate.prepared_&&candidate.sourceOwner_==impl_->owner&&candidate.sourceRevision_==impl_->revision,"joint moment mesh candidate is consumed, foreign or stale");
    candidate.prepared_->validate(mesh,domain);InterfaceParametricMotion::validateDomain(mesh,domain,*candidate.prepared_->phase->interface);
    impl_.swap(candidate.prepared_);candidate.prepared_.reset();candidate.sourceOwner_.reset();
}
}
