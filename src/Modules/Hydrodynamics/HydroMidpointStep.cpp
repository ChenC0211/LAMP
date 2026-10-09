#include "Modules/Hydrodynamics/HydroMidpointStep.h"
#include "Modules/Hydrodynamics/Detail/HydroMidpointIteration.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <iomanip>
#include <sstream>

namespace Quantum::Hydrodynamics {
namespace Detail {
HydroMidpointStepResult solveMidpointIterations(QuadtreeMesh& mesh,::Hydrodynamics& hydro,
    const HydroMidpointIterationOptions& options,const MidpointTrial& trial,const std::function<void()>& publish)
{
    if(options.maximumIterations<=0 || options.maximumIterations>100
        || !std::isfinite(options.absoluteVelocityTolerance) || options.absoluteVelocityTolerance<=0
        || !std::isfinite(options.relativeVelocityTolerance) || options.relativeVelocityTolerance<=0 || !trial)
        throw std::invalid_argument("invalid midpoint iteration controls");
    const auto beginning=hydro.state();
    VEdgeField<double> oldX(mesh,0),guessX(mesh,0);HEdgeField<double> oldY(mesh,0),guessY(mesh,0);
    for(int row=0;row<mesh.vEdgeIndexNum();++row)oldX.value[row]=guessX.value[row]=beginning.velocityX[row];
    for(int row=0;row<mesh.hEdgeIndexNum();++row)oldY.value[row]=guessY.value[row]=beginning.velocityY[row];
    double defect=0,threshold=0;
    for(int iteration=1;iteration<=options.maximumIterations;++iteration){
        auto transaction=hydro.beginStepTransaction(mesh);
        auto candidate=trial(iteration,guessX,guessY);const auto report=candidate.report();
        if(!report.success)throw std::runtime_error("midpoint solve: "+report.stage+": "+report.reason);
        const auto end=hydro.state();defect=0;double scale=0;
        for(int row=0;row<mesh.vEdgeIndexNum()+mesh.hEdgeIndexNum();++row){
            const bool horizontal=row>=mesh.vEdgeIndexNum();const int index=horizontal?row-mesh.vEdgeIndexNum():row;
            const double first=horizontal?oldY.value[index]:oldX.value[index],last=horizontal?end.velocityY[index]:end.velocityX[index];
            auto& guess=horizontal?guessY.value[index]:guessX.value[index];const double center=.5*(first+last);
            defect=std::max(defect,std::abs(center-guess));scale=std::max({scale,std::abs(center),std::abs(first),std::abs(last)});
            guess=center;
        }
        threshold=options.absoluteVelocityTolerance+options.relativeVelocityTolerance*scale;
        if(defect<=threshold){candidate.accept(mesh);if(publish)publish();transaction.commit();
            return {report,iteration,defect};}
        transaction.rollback();
    }
    std::ostringstream message;message<<std::setprecision(17)<<"midpoint did not converge in "<<options.maximumIterations
        <<" iterations: velocity defect "<<defect<<", tolerance "<<threshold;throw std::runtime_error(message.str());
}
}
HydroMidpointStepResult advanceSinglePhaseMidpointStep(QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    ::Hydrodynamics& hydro,const HydroStepTime& time,const HydroSolveOptions& inputOptions,const HydroForceDensity& force,
    const HydroMidpointIterationOptions& iteration)
{
    domain.validateMeshIdentity(mesh);const auto beginning=hydro.state();
    const auto kind=HydroStepGeometry::complete(domain.axisymmetric());
    if(!domain.uniformLeaves() || !beginning.initialStateSet || !beginning.materials
        || beginning.materials->evolution!=HydroMaterialEvolution::Frozen || beginning.progress.geometry!=kind.kind
        || time.policy!=HydroTimePolicy::ConservativeOneStep)
        throw std::invalid_argument("single-phase midpoint requires initialized uniform XY/RZ and frozen owned materials");
    const double density=beginning.density[0];
    for(int c=0;c<beginning.density.size();++c)if(beginning.density[c]!=density)
        throw std::invalid_argument("single-phase midpoint requires constant density");
    auto solver=inputOptions;solver.compatibleNativePressure=true;solver.conservativeMidpoint=true;
    return Detail::solveMidpointIterations(mesh,hydro,iteration,[&](int,const auto& x,const auto& y){
        std::vector<double> velocity(domain.nativeFaceCount());
        for(int row=0;row<domain.nativeFaceCount();++row)velocity[row]=row<mesh.vEdgeIndexNum()?x.value[row]:y.value[row-mesh.vEdgeIndexNum()];
        const auto mass=HydroConservativeTransport::singlePhase(mesh,domain,time,density,velocity);
        const HydroMidpointData midpoint{&domain,nullptr,&x,&y,0,0};
        return hydro.prepareStep(mesh,{kind,time,solver,force,nullptr,&mass,nullptr,&midpoint});
    },{});
}
}
