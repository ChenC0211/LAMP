#include "Modules/Hydrodynamics/HydroParametricStep.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/Detail/HydroMidpointIteration.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <optional>

namespace Quantum::Hydrodynamics {
HydroParametricStepResult advanceParametricCapillaryStep(QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    ::Hydrodynamics& hydro,const InterfaceState& old,const InterfaceStepTime& time,const HydroParametricStepOptions& options,
    const HydroParametricStepObserver& observer)
{
    const auto require=[](bool ok,const char* why){if(!ok)throw std::invalid_argument(why);};
    InterfaceParametricMotion::validateDomain(mesh,domain,old);
    require(!domain.axisymmetric() && options.maximumIterations>0 && options.maximumIterations<=100,
        "parameter midpoint requires XY and a finite iteration budget");
    require(std::isfinite(options.liquidDensity) && options.liquidDensity>0 && std::isfinite(options.gasDensity) && options.gasDensity>0
        && std::isfinite(options.liquidViscosity) && options.liquidViscosity>=0 && std::isfinite(options.gasViscosity) && options.gasViscosity>=0
        && std::isfinite(options.surfaceTension) && options.surfaceTension>=0,"invalid midpoint phase properties");
    require(std::isfinite(options.absoluteVelocityTolerance) && options.absoluteVelocityTolerance>0
        && std::isfinite(options.relativeVelocityTolerance) && options.relativeVelocityTolerance>0,"invalid midpoint convergence tolerance");
    const auto beginning=hydro.state();
    require(beginning.progress.completedSteps==old.progress().completedSteps &&
        std::abs(beginning.progress.endTime-old.progress().endTime)<=64*std::numeric_limits<double>::epsilon()
            *std::max({time.dt,std::abs(beginning.progress.endTime),std::abs(old.progress().endTime)}),
        "midpoint fluid and interface do not share the accepted physical layer");
    const auto kind=HydroStepGeometry::complete(false);const auto pressure=HydroConservativePressure::buildDivergenceAdjoint(mesh,domain);
    auto solver=options.solver;solver.compatibleNativePressure=true;solver.conservativeMidpoint=true;
    const InterfaceStepTime half{time.index,time.begin,time.begin+.5*time.dt,.5*time.dt};
    std::optional<InterfaceParametricMotion::Result> next;
    std::optional<HydroConservativeTransport> mass;ParametricCapillaryReport surface;
    const HydroMidpointIterationOptions iterationOptions{options.maximumIterations,
        options.relativeVelocityTolerance,options.absoluteVelocityTolerance};
    const auto result=Detail::solveMidpointIterations(mesh,hydro,iterationOptions,[&](int iteration,const auto& guessX,const auto& guessY){
        if(observer)observer(iteration,guessX,guessY);
        const HydroSplineVelocity velocity(mesh,domain,guessX,guessY);const auto flux=HydroInterfaceTransport::normalFlux(velocity);
        const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);
        auto middle=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,old,half,flow,flux,flux,options.geometry);
        next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,old,time,flow,flux,flux,options.geometry);
        mass=HydroConservativeTransport::build(mesh,domain,old,next->state,time,options.liquidDensity,options.gasDensity,&middle.state);
        const auto load=buildParametricCapillaryLoad(mesh,domain,pressure,middle.state,options.surfaceTension,&surface);
        hydro.updateMaterials(mesh,{kind,time.begin,HydroMaterialEvolution::TimeDependent,
            HydroParametricMaterial{&domain,&old,options.liquidDensity,options.gasDensity,options.liquidViscosity,options.gasViscosity}});
        HydroForceDensity force;force.faceX=&load.faceForceX;force.faceY=&load.faceForceY;force.evaluationTime=half.end;
        force.timeLevel=HydroForceTimeLevel::StepMidpoint;force.source=hydro.captureParametricForceSource(mesh,domain,middle.state);
        HydroMidpointData midpoint{&domain,&middle.state,&guessX,&guessY,options.liquidViscosity,options.gasViscosity};
        return hydro.prepareStep(mesh,{kind,{time.index,time.begin,time.dt,HydroTimePolicy::ConservativeOneStep},solver,
            force,nullptr,&*mass,&middle.state,&midpoint});
    },[&]{
        hydro.updateMaterials(mesh,{kind,time.end,HydroMaterialEvolution::TimeDependent,
            HydroParametricMaterial{&domain,&next->state,options.liquidDensity,options.gasDensity,options.liquidViscosity,options.gasViscosity}});
    });
    long double kinetic=0;const auto end=hydro.state();
    for(int row=0;row<domain.nativeFaceCount();++row){const double u=row<mesh.vEdgeIndexNum()?end.velocityX[row]:end.velocityY[row-mesh.vEdgeIndexNum()];
        kinetic+=.5L*mass->massAfter()[row]*u*u;}
    return {std::move(next->state),result.hydro,surface,result.iterations,result.midpointVelocityDefect,next->inventoryBalance,static_cast<double>(kinetic)};
}
}
