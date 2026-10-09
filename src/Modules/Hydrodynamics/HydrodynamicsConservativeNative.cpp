#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativePressure.h"
#include "Modules/Hydrodynamics/HydroConservativeStress.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/QuadtreeMesh.h"

void Hydrodynamics::solveConservativeNativeStep(QuadtreeMesh& mesh,double dt,bool rz,
    const Quantum::Hydrodynamics::HydroConservativeTransport& transport,
    const VEdgeField<double>* faceX,const HEdgeField<double>* faceY,const HydroMidpointData* midpointData)
{
    using namespace Quantum::Hydrodynamics;
    const auto require=[](bool ok,const char* why){if(!ok)throw std::invalid_argument(why);};
    require(activeOptions_.boundaryPatches.empty(),"mixed conservative boundary patches are not yet bound");
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
    const bool compatible=activeOptions_.compatibleNativePressure;
    const bool midpointScheme=activeOptions_.conservativeMidpoint;
    require(midpointScheme==bool(midpointData),"midpoint scheme lost its iterate");
    const auto pressure=compatible?HydroConservativePressure::buildDivergenceAdjoint(mesh,domain):HydroConservativePressure::build(mesh,domain);
    const int count=domain.nativeFaceCount(),cells=mesh.cellIndexNum(),nx=mesh.vEdgeIndexNum();
    std::vector<HydroFittedDomain::Boundary> boundary(count);std::vector<unsigned char> free(count,1);
    std::vector<double> cellMu(cells),vertexMu(mesh.vertIndexNum()),oldPressure(cells),exterior(count,0),volume(count);
    Eigen::VectorXd old(count),body=Eigen::VectorXd::Zero(count);
    for(int c=0;c<cells;++c){cellMu[c]=mu.value[c];oldPressure[c]=p.value[c];}
    for(int v=0;v<mesh.vertIndexNum();++v)vertexMu[v]=vertMu.value[v];
    for(const auto f:domain.faces(mesh)) {
        const int row=f.nativeRow;old[row]=row<nx?ux.value[row]:uy.value[row-nx];
        const auto& dual=pressure.duals().volumes()[row];volume[row]=dual.measure;
        require(dual.measure==transport.duals().volumes()[row].measure,"mixed momentum and pressure dual metrics differ");
        if(f.outwardSign!=0) {
            // Existing complete-domain closed boundary: zero normal velocity
            // and zero tangential stress, including the true RZ axis.
            boundary[row].type=HydroFittedDomain::BoundaryType::Symmetry;free[row]=0;
            require(old[row]==0,"mixed closed normal velocity is not constrained");
        }
        for(int owner:{f.negativeCell,f.positiveCell})if(owner>=0) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[owner]);
            const auto a=mesh.verts(cell.cellToVert[0]).coord(),b=mesh.verts(cell.cellToVert[2]).coord();
            const Vector2D lo{std::max(a.x,dual.lower.x),std::max(a.y,dual.lower.y)},
                hi{std::min(b.x,dual.upper.x),std::min(b.y,dual.upper.y)};
            const auto measure=CartesianMeasure::cell((hi.x-lo.x)*(hi.y-lo.y),(hi.x+lo.x)/2,rz);
            body[row]+=measure*(f.horizontal?externalFy.value[owner]:externalFx.value[owner]);
        }
        if(faceX)body[row]+=dual.measure*(f.horizontal?faceY->value[row-nx]:faceX->value[row]);
        if(midpointScheme)body[row]+=static_cast<double>(.5L*(transport.massAfter()[row]-transport.massBefore()[row]))
            *(f.horizontal?activeOptions_.accelerationY:activeOptions_.accelerationX);
    }
    activeReport_.stage="physical-native-assembly";
    Eigen::SparseMatrix<double> viscous;
    if(midpointScheme && midpointData->interface){HydroConservativeStress::EnergyReport report;
        viscous=-HydroConservativeStress::assembleParametricEnergy(mesh,*midpointData->domain,*midpointData->interface,
            midpointData->viscosityNegative,midpointData->viscosityPositive,&report);
        activeReport_.viscousMeasureFractionAdjustment=report.maximumMeasureFractionAdjustment;
    } else if(compatible)viscous=-HydroConservativeStress::assembleCellEnergy(mesh,domain,cellMu);
    else {const auto stress=HydroConservativeStress::build(mesh,domain,pressure,boundary);
        viscous=stress.assemble(mesh,domain,cellMu,vertexMu);}
    const Eigen::SparseMatrix<double> divergence=pressure.integratedDivergence();
    const auto historical=pressure.apply(mesh,domain,oldPressure,exterior);
    Eigen::VectorXd physicalRhs=.5*(viscous*old)+body;
    std::vector<long double> centeredOutward;
    if(midpointScheme){std::vector<double> center(count);
        for(int row=0;row<count;++row)center[row]=row<nx?midpointData->velocityX->value[row]:midpointData->velocityY->value[row-nx];
        centeredOutward=transport.centered(center).outward;
        for(int row=0;row<count;++row)physicalRhs[row]-=static_cast<double>(centeredOutward[row]/dt);}
    std::vector<Eigen::Triplet<double>> entries;
    for(int column=0;column<viscous.outerSize();++column)
        for(Eigen::SparseMatrix<double>::InnerIterator it(viscous,column);it;++it)
            if(free[it.row()])entries.emplace_back(static_cast<int>(it.row()),static_cast<int>(it.col()),-.5*it.value());
    for(int row=0;row<count;++row) {
        if(free[row]) {
            entries.emplace_back(row,row,static_cast<double>(transport.massAfter()[row]/dt));
            physicalRhs[row]+=static_cast<double>(transport.massBefore()[row]*old[row]/dt);
        } else {entries.emplace_back(row,row,1);physicalRhs[row]=0;}
    }
    const auto& segments=transport.duals().segments();
    if(!midpointScheme)for(std::size_t i=0;i<segments.size();++i) {
        const auto& s=segments[i];const auto& flux=transport.fluxes()[i];
        if(s.negativeNative<0 || s.positiveNative<0)
            require(flux.positive==0 && flux.negative==0,"closed mixed momentum received external mass flux");
        for(int side=0;side<2;++side) {
            const int row=side?s.positiveNative:s.negativeNative,donor=side?s.negativeNative:s.positiveNative;
            if(row<0 || !free[row])continue;
            const auto out=side?flux.negative:flux.positive,in=side?flux.positive:flux.negative;
            if(out>0)entries.emplace_back(row,row,static_cast<double>(out/dt));
            if(in>0){require(donor>=0,"missing closed momentum donor");entries.emplace_back(row,donor,-static_cast<double>(in/dt));}
        }
    }
    Eigen::SparseMatrix<double> momentum(count,count);momentum.setFromTriplets(entries.begin(),entries.end());
    if(compatible)momentum.prune([&](int row,int column,double){return row==column || (free[row] && free[column]);});
    Eigen::VectorXd rhs=physicalRhs;
    for(int row=0;row<count;++row)if(free[row])rhs[row]-=historical[row];
    // CCH solves the complete block directly from the old velocity. The old
    // physical-trace strategy keeps its predictor; both consume the SAME
    // assembled momentum equation and coupled residual/linear solver.
    const Eigen::VectorXd predicted=compatible?old:checkedLinearSolve(momentum,rhs,"momentum",activeReport_.momentumRelativeResidual);
    activeReport_.divergenceBeforeProjection=HydroContinuityBalance::LocalFlux(divergence,predicted).maxDivergence(domain.cellVolumes());
    HydroCoupledBalance::Assembly block(momentum,cells,pressure.gradientLoad().nonZeros()+divergence.nonZeros(),0,compatible);
    for(int row=0;row<count;++row)if(free[row])
        for(HydroConservativePressure::Matrix::InnerIterator it(pressure.gradientLoad(),row);it;++it)
            block.addPressureGradient(row,static_cast<int>(it.col()),it.value());
    for(int column=0;column<divergence.outerSize();++column)
        for(Eigen::SparseMatrix<double>::InnerIterator it(divergence,column);it;++it)
            if(!compatible || free[it.col()])block.addDivergence(static_cast<int>(it.row()),static_cast<int>(it.col()),it.value());
    const auto coupled=block.finish();
    Eigen::VectorXd velocity=predicted,increment=Eigen::VectorXd::Zero(cells),load=Eigen::VectorXd::Zero(count);
    activeReport_.stage="physical-native-coupled";activeReport_.coupledCorrectionUsed=true;
    for(int refinement=0;refinement<4;++refinement) {
        const Eigen::VectorXd defect=divergence*velocity;
        const auto coupledRhs=coupled.rightHandSide(rhs,velocity,load,defect);
        const auto correction=checkedLinearSolve(coupled.matrix(),coupledRhs,"coupled",activeReport_.coupledRelativeResidual);
        velocity+=correction.head(count);increment+=correction.tail(cells);
        for(int row=0;row<count;++row)if(!free[row])velocity[row]=0;
        load=pressure.gradientLoad()*increment;
        for(int row=0;row<count;++row)if(!free[row])load[row]=0;
        const auto balance=HydroCoupledBalance::momentumResidual(momentum,velocity,load,rhs,[&](int row){return free[row]?volume[row]:0.;});
        const HydroContinuityBalance::LocalFlux continuity(divergence,velocity);
        activeReport_.momentumBalanceRelativeResidual=balance.relative;
        activeReport_.momentumBalanceMaxForceDensity=balance.maxForceDensity;
        if(balance.relative<=activeOptions_.completeConservation->momentumBalance
            && continuity.relativeResidual()<=activeOptions_.completeConservation->initialContinuity)break;
        require(refinement<3,"mixed physical momentum/continuity failed after coupled residual refinement");
    }
    require(velocity.allFinite() && increment.allFinite(),"nonfinite mixed conservative state");
    std::vector<double> finalPressure(cells);
    for(int c=0;c<cells;++c)finalPressure[c]=oldPressure[c]+increment[c];
    if(compatible) {
        long double integral=0,measure=0;
        for(int c=0;c<cells;++c){const double v=domain.cellVolumes()[domain.cellRowByIndex()[c]];integral+=v*static_cast<long double>(finalPressure[c]);measure+=v;}
        const double mean=static_cast<double>(integral/measure);
        for(int c=0;c<cells;++c){finalPressure[c]-=mean;increment[c]=finalPressure[c]-oldPressure[c];}
    }
    const auto finalLoad=pressure.apply(mesh,domain,finalPressure,exterior);
    Eigen::VectorXd actualLoad(count);
    for(int row=0;row<count;++row)actualLoad[row]=free[row]?finalLoad[row]:0;
    const auto physical=HydroCoupledBalance::momentumResidual(momentum,velocity,actualLoad,physicalRhs,
        [&](int row){return free[row]?volume[row]:0.;});
    require(physical.relative<=activeOptions_.completeConservation->momentumBalance,"mixed final physical pressure/momentum balance failed");
    activeReport_.momentumBalanceRelativeResidual=physical.relative;
    activeReport_.momentumBalanceMaxForceDensity=physical.maxForceDensity;
    const Eigen::VectorXd midpoint=.5*(old+velocity),viscousForce=viscous*midpoint;
    activeReport_.physicalNativeMomentumUsed=true;
    activeReport_.compatibleNativePressureUsed=compatible;
    activeReport_.midpointNativeMomentumUsed=midpointScheme;
    activeReport_.physicalPressurePower=midpoint.dot(actualLoad);
    activeReport_.physicalViscousPower=midpoint.dot(viscousForce);
    long double energyRate=0;
    for(int row=0;row<count;++row)energyRate+=(transport.massAfter()[row]*velocity[row]*velocity[row]
        -transport.massBefore()[row]*old[row]*old[row])/(2*dt);
    activeReport_.physicalEnergyRate=static_cast<double>(energyRate);
    if(midpointScheme){long double convection=0,massDefect=0;
        for(int row=0;row<count;++row){convection+=midpoint[row]*centeredOutward[row]/dt;
            massDefect+=(transport.massAfter()[row]-transport.massBefore()[row])*(velocity[row]-old[row])*(velocity[row]-old[row])/(8*dt);}
        activeReport_.physicalConvectionPower=static_cast<double>(convection);
        activeReport_.physicalExternalPower=midpoint.dot(body);
        activeReport_.midpointMassEnergyDefect=static_cast<double>(massDefect);}
    require(std::isfinite(activeReport_.physicalPressurePower) && std::isfinite(activeReport_.physicalViscousPower)
        && std::isfinite(activeReport_.physicalEnergyRate),"mixed power diagnostic overflow");
    activeReport_.conservativeMomentumUsed=true;activeReport_.conservativeMassRelativeClosure=transport.maxRelativeClosure();
    for(int row=0;row<count;++row) {
        if(row<nx){uxPre.value[row]=old[row];uxMid.value[row]=predicted[row];ux.value[row]=velocity[row];}
        else {uyPre.value[row-nx]=old[row];uyMid.value[row-nx]=predicted[row];uy.value[row-nx]=velocity[row];}
    }
    for(int c=0;c<cells;++c){p.value[c]=finalPressure[c];pressurePotential_.value[c]=increment[c];}
}
