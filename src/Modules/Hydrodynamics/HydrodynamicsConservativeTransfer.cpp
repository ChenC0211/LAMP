#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"

std::unique_ptr<Hydrodynamics> Hydrodynamics::prepareConservativeMeshTransfer(
    const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,
    const InterfaceState& current,const InterfaceState& previous,
    QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& targetCurrent,const InterfaceState& targetPrevious,
    double liquidDensity,double gasDensity,double liquidViscosity,double gasViscosity,
    Quantum::Hydrodynamics::HydroConservativeTransport::StateRemapReport* output)
{
    using Report=Quantum::Hydrodynamics::HydroConservativeTransport::StateRemapReport;
    const auto require=[](bool ok,const char* why) {
        if(!ok)throw std::invalid_argument(std::string("[HydroConservativeTransfer] ")+why);
    };
    require(stepProgress_.completedSteps>0 && pressureMode_==PressureMode::Incremental
        && activeReport_.success && activeOptions_.boundaryPatches.empty()
        && frozenBoundaryPatches_.empty(),"transfer requires accepted closed conservative NS history");
    require(std::isfinite(liquidViscosity) && liquidViscosity>=0
        && std::isfinite(gasViscosity) && gasViscosity>=0,"invalid phase viscosity");
    const auto checkState=[&](const QuadtreeMesh& m,const NativeLeafDomain& d,const InterfaceState& s) {
        s.validateMesh(m);
        require(s.stage()==InterfaceState::Stage::Ready && s.hasConservativePhase(),"missing ready Interface time layer");
        s.conservativePhase().validate(m,d);
    };
    checkState(oldMesh,oldDomain,current);checkState(oldMesh,oldDomain,previous);
    checkState(mesh,domain,targetCurrent);checkState(mesh,domain,targetPrevious);
    validateStateGeometry(oldMesh,HydroStepGeometry::complete(oldDomain.axisymmetric()));
    const auto& old=current.conservativePhase();const auto& pre=previous.conservativePhase();
    const auto& next=targetCurrent.conservativePhase();const auto& nextPre=targetPrevious.conservativePhase();
    require(old.follows(pre,{stepProgress_.completedSteps,pre.progress().endTime,stepProgress_.endTime,fixedDt_}),
        "previous C is not the actual source of the most recently accepted physical step");
    require(next.transferredFrom(old) && nextPre.transferredFrom(pre),"missing current/previous C transfer lineage");
    require(next.progress().completedSteps==old.progress().completedSteps && next.progress().endTime==old.progress().endTime
        && nextPre.progress().completedSteps==pre.progress().completedSteps && nextPre.progress().endTime==pre.progress().endTime,
        "zero-time transfer changed Interface history clocks");
    for(int c:oldDomain.cellIndexByRow()) {
        const long double fraction=old.liquidVolumes()[c]/old.capacities()[c];
        const auto expected=fraction*liquidViscosity+(1-fraction)*gasViscosity;
        require(std::isfinite(mu.value[c]) && std::abs(mu.value[c]-expected)
            <=128*std::numeric_limits<double>::epsilon()*expected,"source viscosity differs from supplied C mixture");
    }
    Report report;
    // This also checks the accepted source, its current material clock and
    // actual density against old C. Its temporary solve/cache is restored.
    report.current=projectConservativeMomentumRemap(oldMesh,oldDomain,old,mesh,domain,next,liquidDensity,gasDensity);
    std::vector<double> historical(oldDomain.nativeFaceCount());
    for(int row=0;row<oldDomain.nativeFaceCount();++row)historical[row]=row<oldMesh.vEdgeIndexNum()
        ?uxPre.value[row]:uyPre.value[row-oldMesh.vEdgeIndexNum()];
    report.previous=projectConservativeMomentumLayer(oldMesh,oldDomain,pre,historical,mesh,domain,nextPre,liquidDensity,gasDensity);
    CellField<double> pressure(mesh,0);VEdgeField<double> vx(mesh,0);HEdgeField<double> vy(mesh,0);
    long double pressureScale=0;
    for(int c:oldDomain.cellIndexByRow()) {
        const auto value=static_cast<long double>(p.value[c])*oldDomain.cellVolumes()[oldDomain.cellRowByIndex()[c]];
        require(std::isfinite(value),"nonfinite source physical pressure");
        report.sourcePressureIntegral+=value;pressureScale+=std::abs(value);
    }
    // Refinement-only nested cells: preserve the source cell average and its
    // volume integral. No analytic Laplace jump, regauging, or remap multiplier.
    for(int c:domain.cellIndexByRow()) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const int label=oldMesh.findCellLabel(cell.centerCoord());
        require(label>=0 && label<oldMesh.cellLabelCapacity(),"pressure transfer lost its source parent");
        const int parent=oldMesh.cells(label).index();
        require(parent>=0 && parent<oldMesh.cellIndexNum(),"pressure source is not a native leaf");
        pressure.value[c]=p.value[parent];
        report.targetPressureIntegral+=static_cast<long double>(pressure.value[c])*domain.cellVolumes()[domain.cellRowByIndex()[c]];
    }
    require(std::abs(report.targetPressureIntegral-report.sourcePressureIntegral)
        <=4096*std::numeric_limits<double>::epsilon()*pressureScale,"physical pressure integral changed beyond measure roundoff");
    for(int row=0;row<domain.nativeFaceCount();++row) {
        if(row<mesh.vEdgeIndexNum())vx.value[row]=report.current.velocity[row];
        else vy.value[row-mesh.vEdgeIndexNum()]=report.current.velocity[row];
    }
    auto result=std::make_unique<Hydrodynamics>(mesh);
    const auto geometry=HydroStepGeometry::complete(domain.axisymmetric());
    result->initializeState(mesh,geometry,{&pressure,&vx,&vy});
    for(int row=0;row<domain.nativeFaceCount();++row) {
        if(row<mesh.vEdgeIndexNum())result->uxPre.value[row]=report.previous.velocity[row];
        else result->uyPre.value[row-mesh.vEdgeIndexNum()]=report.previous.velocity[row];
    }
    result->stepProgress_=stepProgress_;result->fixedDt_=fixedDt_;
    result->pressureMode_=pressureMode_;result->conservativeTransportActive_=true;
    result->activeOptions_=activeOptions_;result->boundaryConfigurationFixed_=boundaryConfigurationFixed_;
    result->materialState_=materialState_;
    result->updateMaterials(mesh,{geometry,stepProgress_.endTime,HydroMaterialEvolution::TimeDependent,
        HydroConservativeMaterial{&domain,&next,liquidDensity,gasDensity,liquidViscosity,gasViscosity}});
    result->maxDivURes_=report.current.divergenceAfter;result->uMagMax_=0;
    for(double u:report.current.velocity)result->uMagMax_=std::max(result->uMagMax_,std::abs(u));
    for(int c:domain.cellIndexByRow())result->uMagMax_=std::max(result->uMagMax_,result->uMag.value[c]);
    result->activeReport_=HydroSolveReport{};result->activeReport_.success=true;
    result->activeReport_.stage="conservative-mesh-transfer";
    if(output)*output=std::move(report);
    return result;
}
