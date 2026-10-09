#include "Electrohydrodynamics/Coupling/DielectricCoupling.h"
#include "Hydrodynamics/IO/HydroCaseInput.h"
#include "Modules/Electricity/FittedChargeEvolution.h"
#include "Modules/Hydrodynamics/HydroUniformLoad.h"
#include "Modules/Interface/Compatibility/InterfaceSurfaceTransfer.h"
#include "Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/OutputLayout.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <algorithm>
#include <cmath>

namespace {
using namespace Quantum::Electricity;
using Solver=PolarizationPoissonFVSolverCellCentered;
void require(bool valid,const char* message) {if(!valid)throw std::invalid_argument(std::string("[ConservativeEHD] ")+message);}
class ConservativeColumnCoupling final:public HydroCaseCoupling {
    struct Completed {
        FittedChargeCollectionCandidate electric;
        FittedSharpDielectricTractions::Diagnostics traction;
        Quantum::Hydrodynamics::ConservativeTractionReport tailForce;
        bool afterSectionEvent=false;
        Quantum::Hydrodynamics::HydroUniformLoad force;
        explicit Completed(const QuadtreeMesh& mesh):force(mesh){}
    };
    class Candidate final:public HydroConservativeCouplingCandidate {
        ConservativeColumnCoupling& owner_;
        std::unique_ptr<Completed> prepared_;
    public:
        Candidate(ConservativeColumnCoupling& owner,std::unique_ptr<Completed> prepared)
            :owner_(owner),prepared_(std::move(prepared)){}
        const Quantum::Hydrodynamics::HydroUniformLoad& load() const override {return prepared_->force;}
        void commit() noexcept override {owner_.current_.swap(prepared_);}
    };
    DielectricConfig config_;
    bool axisCapModel_;
    std::shared_ptr<const HydroCaseSchedule> schedule_;
    FaceSegmentGraph faces_;
    std::vector<Solver::FittedBoundary> boundaries_;
    std::unique_ptr<Completed> current_;
    const QuadtreeMesh* mesh_=nullptr;
    std::ofstream history_;
    std::ofstream componentHistory_;
    int published_=-1;
    FittedChargeProblem problem() const {
        FittedChargeParameters parameters{config_,{config_.conductivityNegative,config_.conductivityPositive},
            config_.chargeClosureTolerance.value_or(16e-9)};
        parameters.field.chargeClosureTolerance.reset();
        return {faces_,boundaries_,std::move(parameters)};
    }
    void write(int step,double time,double dt,bool appendHistory=true,bool forceFields=false) {
        const auto& e=current_->electric;long double total=0;
        std::size_t refinedPanels=0,maxIntervals=0,rationalArc=0,rationalArea=0;double measureError=0;
        bool detached=false;double sectionResidual=0;
        for(const auto& report:e.chargeSteps){detached|=report.sectionDetached;sectionResidual+=report.sectionRemapConservationResidual;}
        for(std::size_t s=0;s<e.state.componentCount();++s) {
            const auto charge=e.state.charge(s);long double child=0;
            for(double value:charge->transport().panelCharge())child+=value;total+=child;
            rationalArc+=charge->measureInterpolation().rationalArcIntervals;
            rationalArea+=charge->measureInterpolation().rationalAreaIntervals;
            for(const auto& panel:charge->geometry().panels()) {
                refinedPanels+=panel.measureQuadratureIntervals>1;maxIntervals=std::max(maxIntervals,panel.measureQuadratureIntervals);
                for(double error:panel.measureQuadratureRelativeError)measureError=std::max(measureError,error);
            }
            if(appendHistory)componentHistory_<<std::setprecision(17)<<step<<','<<time<<','<<e.state.surfaceId(s)<<','
                <<child<<','<<charge->transport().panelCharge().size()<<'\n';
        }
        if(appendHistory) {history_<<std::setprecision(17)<<step<<','<<time<<','<<dt<<','<<total<<','
            <<e.report.relativeChargeClosure<<','<<detached<<','<<sectionResidual<<','<<e.state.electric().report.relativeResidual<<','
            <<current_->force.maxForce<<','<<current_->traction.axisCapExteriorSamples<<','
            <<current_->traction.maxAxisCapRelativeResidual<<','<<rationalArc<<','<<rationalArea
            <<','<<refinedPanels<<','<<maxIntervals<<','<<measureError
            <<','<<current_->tailForce.tailAlternatives<<','<<current_->tailForce.maxScaledTailForceUncertainty
            <<','<<current_->tailForce.maxTailForceUncertainty<<','<<e.state.electric().report.axisNeckHarmonicFaces
            <<','<<current_->traction.axisNeckHarmonicSamples<<','<<current_->traction.maxAxisNeckRelativeResidual<<'\n';history_.flush();componentHistory_.flush();}
        const bool fields=step==0?schedule_->writeStepFields(0):
            (schedule_->controls.stepRetry?schedule_->writeVariableStepFields(step,time-dt,time):schedule_->writeStepFields(step));
        if(schedule_->controls.writeFields && (forceFields || fields || detached || e.interiorSplit || current_->afterSectionEvent)) {
            std::ofstream panels(OutputLayout::fields()/("electric_charge_"+std::to_string(step)+".csv"));
            panels<<"step,time,surface_id,panel,r,z,area,charge,density\n"<<std::setprecision(17);
            for(std::size_t s=0;s<e.state.componentCount();++s) {
                const auto charge=e.state.charge(s);const auto& q=charge->transport().panelCharge();
                const auto& geometry=charge->geometry().panels();const auto& measure=charge->transport().panelMeasure();
                for(std::size_t i=0;i<q.size();++i)panels<<step<<','<<time<<','<<e.state.surfaceId(s)<<','<<i<<','
                    <<geometry[i].frame.point.x<<','<<geometry[i].frame.point.y<<','<<measure[i]<<','<<q[i]<<','<<q[i]/measure[i]<<'\n';
            }
            std::ofstream potential(OutputLayout::fields()/("electric_potential_"+std::to_string(step)+".csv"));
            potential<<"step,time,r,z,potential\n"<<std::setprecision(17);const auto& phi=e.state.electric().potential;
            for(int i=0;i<phi.size();++i) {
                const auto p=mesh_->cells(mesh_->cellIndexToLabel[i]).centerCoord();
                potential<<step<<','<<time<<','<<p.x<<','<<p.y<<','<<phi[i]<<'\n';
            }
            require(bool(panels) && bool(potential),"electric field output failed after physical acceptance");
        }
        require(bool(history_) && bool(componentHistory_),"electric history output failed after physical acceptance");
    }
public:
    ConservativeColumnCoupling(DielectricConfig config,bool model):config_(std::move(config)),axisCapModel_(model){}
    bool supportsConservative() const override {return true;}
    void writeConservativeEventSource(int step,double time) override {
        require(current_ && current_->electric.state.completedStep()==step && current_->electric.state.time()==time,
            "pre-event electric output requires the actual accepted source");
        write(step,time,current_->electric.report.dt,false,true);
    }
    void preflight(const HydrodynamicsConfig& flow,const InterfaceNodalConfig& interface,
        const DictionaryGeometry&,const DictionaryMesh&,int interfaceIndex,std::shared_ptr<const HydroCaseSchedule> schedule) override {
        require(schedule && flow.coordinateSystem=="axisymmetricRZ" && flow.materialModel=="twoPhase"
            && interface.transport=="conservative" && interface.advect && interface.reinitializeEvery==0
            && flow.interfaceForceGeometry=="conservative" && flow.interfaceForceLayer=="inventory"
            && interfaceIndex<0 && interface.neckBreakup.has_value(),
            "requires the conservative RZ reflection-column breakup case");
        require(config_.model=="leakyDielectric" && config_.interfaceDiscretization=="sharpHermiteP3"
            && config_.chargeConvection && config_.lowerPotential==0 && config_.upperPotential==0,
            "grounded radial-wall column requires sharp leaky dielectric, moving charge and zero supplied electrode potentials");
        config_.axisymmetric=true;config_.nativeAxisNeckHarmonics=axisCapModel_;schedule_=std::move(schedule);
        preflightFittedChargeStep(problem().parameters,schedule_->controls.dt);
        std::cout<<"[ConservativeEHD] Grounded radial wall; axis and axial reflection sections: zero normal electric current.\n"
            <<"[ConservativeEHD] Drift interface/material, old conduction current once, new P3/Maxwell, NS, joint acceptance.\n";
        if(axisCapModel_)std::cout<<"[ConservativeEHD] EXPERIMENTAL regular harmonic axis-cap traction enabled; residual is reported, model accuracy is not accepted.\n";
        if(config_.nativeAxisNeckHarmonics)std::cout<<"[ConservativeEHD] EXPERIMENTAL thin-neck harmonic transmission: certified liquid-axis window, regular interior/logarithmic exterior, conservative FV/current and Maxwell.\n";
    }
    void initializeConservative(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state) override {
        require(!current_ && domain.axisymmetric() && mesh.minCoord.x==0 && state.topology().surfaces().size()==1,
            "initialization requires one retained RZ column on the shared mesh");
        const auto& panels=state.topology().surfaces()[0].geometry;
        for(bool end:{false,true})require(panels.nativeChain().endpointSymmetryPlane(end)!=nullptr,
            "initial column must meet its two declared reflection sections");
        faces_.build(mesh);boundaries_.resize(faces_.segments().size());
        for(std::size_t i=0;i<boundaries_.size();++i) {
            const auto& f=faces_.segments()[i];if(f.L>=0)continue;
            boundaries_[i].type=f.boundary==Boundary::east?Solver::FittedBoundaryType::Dirichlet:Solver::FittedBoundaryType::HomogeneousNeumann;
            boundaries_[i].value=0;
        }
        CellField<double> seed(mesh,0.);auto ready=std::make_unique<Completed>(mesh);
        ready->electric.state=initializeFittedChargeState(mesh,problem(),state.geometry(),
            std::vector<FittedChargeComponentInitial>{{state.topology().surfaces()[0].id,&panels,config_.initialSurfaceChargeDensity}},seed);
        current_=std::move(ready);mesh_=&mesh;
        history_.open(OutputLayout::series()/"conservative_electric.csv");
        history_<<"step,time,dt,total_charge,charge_closure,section_detached,section_remap_residual,field_residual,max_force,axis_model_samples,axis_model_residual,rational_arc_intervals,rational_area_intervals,refined_measure_panels,max_measure_intervals,measure_relative_error_estimate,tail_alternatives,max_scaled_tail_force_uncertainty,max_tail_force_uncertainty,axis_neck_harmonic_faces,axis_neck_harmonic_samples,axis_neck_model_residual\n";
        componentHistory_.open(OutputLayout::series()/"conservative_electric_components.csv");
        componentHistory_<<"step,time,surface_id,total_charge,panel_count\n";
        write(0,0,0);published_=0;
    }
    std::unique_ptr<HydroConservativeCouplingCandidate> prepareConservative(const QuadtreeMesh& mesh,
        const NativeLeafDomain& domain,const InterfaceState& old,const InterfaceState& next,const InterfaceSurfaceTransfer& physical,
        const InterfaceNeckCap::MaterialTransfer* event,const InterfaceStepTime& time,
        const Quantum::Hydrodynamics::HydroConservativePressure* pressure,
        const Quantum::Hydrodynamics::ConservativeTractionGeometry& geometry) override {
        require(current_ && mesh_==&mesh && current_->electric.state.completedStep()==time.index-1
            && current_->electric.state.time()==time.begin && published_==time.index-1,"coupling source clock/publication mismatch");
        auto prepared=std::make_unique<Completed>(mesh);
        prepared->afterSectionEvent=current_->electric.interiorSplit || std::any_of(current_->electric.chargeSteps.begin(),
            current_->electric.chargeSteps.end(),[](const auto& q){return q.sectionDetached;});
        const auto& interior=next.conservativeGeometry().interiorMaterialTransfer;
        const auto bind=[&](const InterfaceSurfaceTransfer::Surface& map,const LevelSetInterfacePanels* target) {
            std::size_t source=0;while(source<current_->electric.state.componentCount()
                && current_->electric.state.surfaceId(source)!=map.sourceId)++source;
            require(source<current_->electric.state.componentCount(),"physical material surface has no charged identity");
            SurfaceChargePreparedTransfer transfer;transfer.source=&current_->electric.state.charge(source)->geometry();
            transfer.target=target;transfer.step=time.index;transfer.sourceTime=time.begin;transfer.targetTime=time.end;transfer.dt=time.dt;
            transfer.maxMaterialDisplacement=map.maxPhysicalDisplacement;transfer.maxProjectionDistance=map.maxNumericalDisplacement;
            for(const auto& row:map.rows) {
                transfer.weights.emplace_back();for(const auto& w:row.overlaps)transfer.weights.back().push_back({w.targetPanel,w.share});
            }
            return transfer;
        };
        if(interior) {
            require(!event && old.topology().surfaces().size()==1 && physical.surfaces().size()==1
                && next.topology().surfaces().size()==2 && interior->time.index==time.index
                && interior->time.begin==time.begin && interior->time.end==time.end && interior->time.dt==time.dt,
                "interior event does not match the coupled material step");
            const auto& connected=interior->sourceState->topology().surfaces().front().geometry;
            physical.validateFor(old.geometry(),interior->sourceState->geometry(),time);
            auto transfer=bind(physical.surfaces().front(),&connected);
            FittedChargeMotion motion{interior->sourceState->geometry(),&connected,{},{},time.dt,time.begin,time.end,&transfer};
            FittedChargeInteriorSplit split;split.geometry=next.geometry();split.cells=interior->cells;split.section=interior->section;
            for(const auto& target:interior->targets) {
                const auto found=std::find_if(next.topology().surfaces().begin(),next.topology().surfaces().end(),
                    [&](const auto& surface){return surface.id==target.surfaceId && surface.liquidId==target.liquidId;});
                require(found!=next.topology().surfaces().end() && target.panels->usesSource(*next.geometry()),
                    "interior event child is not in the coupled candidate");
                split.targets.push_back({target.surfaceId,&found->geometry});
            }
            for(const auto& row:interior->weights) {
                split.weights.emplace_back();for(const auto& w:row)split.weights.back().push_back({w.targetSurface,w.targetPanel,w.share});
            }
            prepared->electric=prepareFittedChargeInteriorSplitStep(mesh,problem(),current_->electric.state,motion,split);
        } else if(next.topology().surfaces().size()>1) {
            require(!event && physical.surfaces().size()==next.topology().surfaces().size(),"incomplete child physical transfer");
            physical.validateFor(old.geometry(),next.geometry(),time);
            std::vector<SurfaceChargePreparedTransfer> transfers;transfers.reserve(physical.surfaces().size());
            std::vector<FittedChargeComponentMotion> motions;
            for(const auto& map:physical.surfaces()) {
                const auto target=std::find_if(next.topology().surfaces().begin(),next.topology().surfaces().end(),
                    [&](const auto& surface){return surface.id==map.targetId && surface.liquidId==map.liquidId;});
                require(target!=next.topology().surfaces().end() && map.sourceId==map.targetId,"child physical identity changed");
                transfers.push_back(bind(map,&target->geometry));
                motions.push_back({map.targetId,{next.geometry(),&target->geometry,{},{},time.dt,time.begin,time.end,&transfers.back()}});
            }
            prepared->electric=prepareFittedChargeCollectionStep(mesh,problem(),current_->electric.state,motions);
        } else {
            require(old.topology().surfaces().size()==1 && next.topology().surfaces().size()==1 && physical.surfaces().size()==1,
                "column coupling requires a complete single-surface material leg");
            const auto& target=next.topology().surfaces()[0].geometry;const auto* physicalTarget=event?event->source.get():&target;
            physical.validateFor(old.geometry(),physicalTarget->nativeChain().source(),time);
            auto transfer=bind(physical.surfaces()[0],physicalTarget);SurfaceChargeSectionDetachment section;
            if(event) {
                section.source=event->source.get();section.target=&target;section.cells=event->cells;
                section.step=time.index;section.time=time.end;section.beginZ=event->beginZ;section.upperSection=event->upperSection;
                for(const auto& row:event->weights) {
                    section.weights.emplace_back();for(const auto& w:row)section.weights.back().push_back({w.targetPanel,w.share});
                }
            }
            FittedChargeMotion motion{next.geometry(),&target,{},{},time.dt,time.begin,time.end,&transfer,event?&section:nullptr};
            motion.targetSurfaceId=next.topology().surfaces()[0].id;
            auto single=prepareFittedChargeStep(mesh,problem(),current_->electric.state,motion);
            prepared->electric.state=std::move(single.state);prepared->electric.report=std::move(single.report);
            prepared->electric.chargeSteps.push_back(std::move(single.chargeStep));
        }
        const auto source=prepared->electric.state.source();
        const auto addTraction=[&](const Quantum::Hydrodynamics::ConservativeTractionGeometry& receiving,
            const LevelSetInterfacePanels& surface) {
            FittedAxisCapTractionModel model{next.geometry(),{}};
            if(axisCapModel_)for(int cell:next.eventCurvatureCells()) {
                const auto& c=mesh.cells(mesh.cellIndexToLabel[cell]);
                if(mesh.verts(c.cellToVert[0]).coord().x==0 && !surface.nativeChain().cellToGraphs()[cell].empty())model.cells.push_back(cell);
            }
            const auto loads=reconstructNativeFittedSharpDielectricTractions(mesh,receiving.trace,surface,problem().parameters.field,
                prepared->electric.state.electric().potential,source,model.cells.empty()?nullptr:&model);
            auto& quality=prepared->traction;const auto& d=loads.diagnostics;
            quality.samplePoints+=d.samplePoints;quality.localTraceFits+=d.localTraceFits;quality.axisCapExteriorSamples+=d.axisCapExteriorSamples;
            quality.axisNeckHarmonicSamples+=d.axisNeckHarmonicSamples;
            quality.maxAxisNeckRelativeResidual=std::max(quality.maxAxisNeckRelativeResidual,d.maxAxisNeckRelativeResidual);
            quality.maxTraceCondition=std::max(quality.maxTraceCondition,d.maxTraceCondition);
            quality.maxTraceReproductionError=std::max(quality.maxTraceReproductionError,d.maxTraceReproductionError);
            quality.maxNormalTraction=std::max(quality.maxNormalTraction,d.maxNormalTraction);
            quality.maxTangentialTraction=std::max(quality.maxTangentialTraction,d.maxTangentialTraction);
            quality.maxAxisCapRelativeResidual=std::max(quality.maxAxisCapRelativeResidual,d.maxAxisCapRelativeResidual);
            Quantum::Hydrodynamics::ConservativeTractionReport tail;
            const auto force=Quantum::Hydrodynamics::buildConservativeInterfaceTractionLoad(mesh,domain,next,pressure,
                receiving,loads.geometry,loads.tractionComponents(),&tail);
            prepared->tailForce.tailAlternatives+=tail.tailAlternatives;
            prepared->tailForce.maxScaledTailForceUncertainty=std::max(prepared->tailForce.maxScaledTailForceUncertainty,tail.maxScaledTailForceUncertainty);
            prepared->tailForce.maxTailForceUncertainty=std::max(prepared->tailForce.maxTailForceUncertainty,tail.maxTailForceUncertainty);
            for(int i=0;i<mesh.vEdgeIndexNum();++i)prepared->force.faceForceX.value[i]+=force.faceForceX.value[i];
            for(int i=0;i<mesh.hEdgeIndexNum();++i)prepared->force.faceForceY.value[i]+=force.faceForceY.value[i];
            for(int i=0;i<mesh.cellIndexNum();++i){prepared->force.forceX.value[i]+=force.forceX.value[i];prepared->force.forceY.value[i]+=force.forceY.value[i];}
            prepared->force.totalForceX+=force.totalForceX;prepared->force.totalForceY+=force.totalForceY;
        };
        if(geometry.components.empty()) {
            require(next.topology().surfaces().size()==1,"multi-child Maxwell receiving support is missing");
            addTraction(geometry,next.topology().surfaces().front().geometry);
        } else {
            require(geometry.components.size()==next.topology().surfaces().size(),"incomplete child Maxwell supports");
            for(const auto& child:geometry.components) {
                const auto surface=std::find_if(next.topology().surfaces().begin(),next.topology().surfaces().end(),
                    [&](const auto& s){return s.liquidId==child->liquidId;});
                require(surface!=next.topology().surfaces().end(),"Maxwell support has no owning child");addTraction(*child,surface->geometry);
            }
        }
        for(int i=0;i<mesh.cellIndexNum();++i)prepared->force.maxForce=std::max(prepared->force.maxForce,
            std::hypot(prepared->force.forceX.value[i],prepared->force.forceY.value[i]));
        return std::make_unique<Candidate>(*this,std::move(prepared));
    }
    void publishCompletedStep(int step,double begin,double dt) override {
        require(current_ && step==published_+1 && current_->electric.state.completedStep()==step
            && current_->electric.state.time()==begin+dt,"electric publication does not match accepted flow clock");
        published_=step;write(step,current_->electric.state.time(),dt);
        if(step%schedule_->controls.logEverySteps==0 || current_->electric.interiorSplit
            || std::any_of(current_->electric.chargeSteps.begin(),current_->electric.chargeSteps.end(),[](const auto& q){return q.sectionDetached;}))
            std::cout<<"[ConservativeEHD] accepted "<<step<<" time="<<current_->electric.state.time()
                <<" charge_closure="<<current_->electric.report.relativeChargeClosure
                <<" axis_model_samples="<<current_->traction.axisCapExteriorSamples
                <<" axis_model_residual="<<current_->traction.maxAxisCapRelativeResidual<<'\n';
    }
    void initialize(const QuadtreeMesh&,const ImmersedGeometry<QuadtreeMesh>*,int,const HydrodynamicsConfig&) override
    {throw std::logic_error("conservative column requires shared conservative initialization");}
    void computeForce(const QuadtreeMesh&,const ImmersedGeometry<QuadtreeMesh>*,int,int,double,double,CellField<double>&,CellField<double>&) override
    {throw std::logic_error("conservative column requires a joint candidate");}
};
}
std::unique_ptr<HydroCaseCoupling> makeConservativeColumnCoupling(DielectricConfig config,bool model)
{return std::make_unique<ConservativeColumnCoupling>(std::move(config),model);}
