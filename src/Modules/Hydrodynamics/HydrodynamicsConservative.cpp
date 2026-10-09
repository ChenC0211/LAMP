#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/Detail/HydroMomentumAssembly.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Hydrodynamics/HydroContinuityBalance.h"
#include <iomanip>
#include <sstream>

void Hydrodynamics::addConservativeMomentum(const QuadtreeMesh& mesh,double dt,bool axisymmetric,
    const Quantum::Hydrodynamics::HydroConservativeTransport& transport,HydroMomentumAssembly& assembly)
{
    const int nx=mesh.vEdgeIndexNum(),count=nx+mesh.hEdgeIndexNum();
    std::vector<double> velocity(count);
    for(int row=0;row<count;++row) velocity[row]=row<nx?ux.value[row]:uy.value[row-nx];
    const auto require=[](bool valid,const char* reason) {if(!valid) throw std::invalid_argument(reason);};
    const auto& segments=transport.duals().segments();
    for(std::size_t i=0;i<segments.size();++i) {
        const auto& s=segments[i];if(s.negativeNative>=0 && s.positiveNative>=0) continue;
        const auto& f=transport.fluxes()[i];
        if(boundaryFaces_.empty()) {
            require(f.positive==0 && f.negative==0,"closed Hydro boundary has a conservative mass flux");continue;
        }
        const auto& b=boundaryFaces_[s.primalNative];
        if(!b.outlet) {
            const auto normal=b.sign*b.value;
            require((normal>0 || f.positive==0) && (normal<0 || f.negative==0),
                "conservative boundary mass direction disagrees with prescribed normal velocity");
        }
    }
    const auto incoming=[&](const auto& segment,int component) {
        require(!boundaryFaces_.empty(),"missing physical conservative momentum inlet");
        const auto& b=boundaryFaces_[segment.primalNative];
        require(b.owner>=0 && !b.outlet && b.value<0,"conservative momentum cannot enter through a wall or pressure outlet");
        if(component==segment.normalAxis) return velocity[segment.primalNative];
        // Reuse the existing tau_nt=0 inlet trace for tangential momentum.
        const int row=segment.negativeNative>=0?segment.negativeNative:segment.positiveNative;
        const bool horizontal=row>=nx;
        const auto& edge=horizontal?mesh.hEdges(mesh.hEdgeIndexToLabel[row-nx]):mesh.vEdges(mesh.vEdgeIndexToLabel[row]);
        if(edge.masterCell<0 || edge.slaveCell<0) {
            // At a domain corner use the same one-sided boundary derivative
            // as the existing outlet half-volume, including constrained rows.
            const auto& owner=boundaryFaces_[row];
            const auto& cross=horizontal?mesh.vEdges(mesh.vEdgeIndexToLabel[segment.primalNative]):
                mesh.hEdges(mesh.hEdgeIndexToLabel[segment.primalNative-nx]);
            const int end=owner.sign>0?0:1;const auto& vertex=mesh.verts(cross.edgeToVert[end]);
            const int label=horizontal?vertex.vertToVEdge[end]:vertex.vertToHEdge[end];
            require(label>=0,"conservative corner inflow lacks an inward boundary trace");
            const int inward=horizontal?mesh.vEdges(label).index():nx+mesh.hEdges(label).index();
            require(boundaryFaces_[inward].owner>=0 && boundaryFaces_[inward].sign==b.sign,
                "conservative corner trace leaves its physical boundary");
            const double derivative=owner.sign*(velocity[segment.primalNative]-velocity[inward])/owner.normalSpacing;
            return HydroBoundaryDetail::zeroShearTrace(velocity[row],owner.tangentSpacing,b.sign,derivative);
        }
        const auto& high=mesh.cells(edge.masterCell);const auto& low=mesh.cells(edge.slaveCell);
        const int side=b.sign>0?1:0;
        const int hi=horizontal?mesh.vEdges(high.cellToVEdge[side][0]).index():nx+mesh.hEdges(high.cellToHEdge[side][0]).index();
        const int lo=horizontal?mesh.vEdges(low.cellToVEdge[side][0]).index():nx+mesh.hEdges(low.cellToHEdge[side][0]).index();
        const double hn=horizontal?high.size().x:high.size().y,ht=horizontal?high.size().y:high.size().x;
        return HydroBoundaryDetail::zeroShearTrace(velocity[row],hn,b.sign,velocity[hi]-velocity[lo],ht);
    };
    std::vector<unsigned char> physical(count,0);
    for(int row=0;row<count;++row) {
        const bool horizontal=row>=nx;
        const auto& edge=horizontal?mesh.hEdges(mesh.hEdgeIndexToLabel[row-nx]):mesh.vEdges(mesh.vEdgeIndexToLabel[row]);
        double legacyMass=0;
        if(edge.masterCell>=0 && edge.slaveCell>=0) {
            const auto& high=mesh.cells(edge.masterCell);const auto& low=mesh.cells(edge.slaveCell);
            const double volume=high.area()*(axisymmetric?2.*Pi*edge.midpointCoord().x:1.);
            legacyMass=volume*.5*(rho.value[high.index()]+rho.value[low.index()]);
        } else {
            if(boundaryFaces_.empty() || !boundaryFaces_[row].outlet) continue; // Constraint, not a momentum row.
            const auto& b=boundaryFaces_[row];legacyMass=rho.value[b.owner]*b.volume;
        }
        // The legacy kernels retain viscosity/force and pressure geometry.
        // Replace ONLY their temporal mass and already-disabled advection.
        physical[row]=1;
        assembly.entries.emplace_back(row,row,static_cast<double>((transport.massAfter()[row]-legacyMass)/dt));
        assembly.rhs[row]+=static_cast<double>((transport.massBefore()[row]-legacyMass)*velocity[row]/dt);
    }
    // Implicit donor velocities for the already integrated phase-mass flux.
    // Interior off-diagonals are nonpositive; the temporal+advective row sum
    // is M_old/dt by dual continuity, not a new explicit mass-outflow CFL.
    // Physical incoming boundary traces remain frozen at the source time.
    for(std::size_t i=0;i<segments.size();++i) {
        const auto& s=segments[i];const auto& f=transport.fluxes()[i];
        for(int side=0;side<2;++side) {
            const int row=side?s.positiveNative:s.negativeNative;
            if(row<0 || !physical[row]) continue;
            const int donor=side?s.negativeNative:s.positiveNative;
            const auto out=side?f.negative:f.positive,in=side?f.positive:f.negative;
            if(out>0) assembly.entries.emplace_back(row,row,static_cast<double>(out/dt));
            if(in==0) continue;
            if(donor>=0) assembly.entries.emplace_back(row,donor,-static_cast<double>(in/dt));
            else assembly.rhs[row]+=static_cast<double>(in*incoming(s,row<nx?0:1)/dt);
        }
    }
    activeReport_.conservativeMomentumUsed=true;
    activeReport_.conservativeMassRelativeClosure=transport.maxRelativeClosure();
}

Quantum::Hydrodynamics::HydroConservativeTransport::ProjectedMomentumRemap
Hydrodynamics::projectConservativeMomentumRemap(const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,
    const InterfaceConservativePhase& source,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& target,double liquidDensity,double gasDensity)
{
    using Transport=Quantum::Hydrodynamics::HydroConservativeTransport;
    const auto require=[](bool ok,const char* why) {if(!ok)throw std::invalid_argument(std::string("[HydroMomentumRemap] ")+why);};
    require(initialStateSet_ && conservativeTransportActive_ && stepMeshIdentity_ && materialState_
        && !stepProgress_.stopped && !pendingStep_ && !stepTransaction_ && !stepProgress_.pendingAcceptance,
        "projection requires an accepted conservative Hydro source, without a pending transaction");
    stepMeshIdentity_->validate(oldMesh);source.validate(oldMesh,oldDomain);
    require(stepProgress_.geometry==HydroStepGeometry::complete(oldDomain.axisymmetric()).kind
        && stepProgress_.completedSteps==source.progress().completedSteps
        && stepProgress_.acceptedSteps==stepProgress_.completedSteps
        && stepProgress_.endTime==source.progress().endTime
        && materialState_->sourceTime==source.progress().endTime,"source geometry/material/accepted clocks differ");
    for(int c:oldDomain.cellIndexByRow()) {
        const long double fraction=source.liquidVolumes()[c]/source.capacities()[c];
        const auto density=fraction*liquidDensity+(1-fraction)*gasDensity;
        require(std::isfinite(density) && density>0 && std::abs(rho.value[c]-density)
            <=128*std::numeric_limits<double>::epsilon()*density,"source density is not the supplied C mixture");
    }
    std::vector<double> oldVelocity(oldDomain.nativeFaceCount());
    for(int row=0;row<oldDomain.nativeFaceCount();++row)oldVelocity[row]=row<oldMesh.vEdgeIndexNum()?
        ux.value[row]:uy.value[row-oldMesh.vEdgeIndexNum()];
    return projectConservativeMomentumLayer(oldMesh,oldDomain,source,oldVelocity,mesh,domain,target,liquidDensity,gasDensity);
}

Quantum::Hydrodynamics::HydroConservativeTransport::ProjectedMomentumRemap
Hydrodynamics::projectConservativeMomentumLayer(const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,
    const InterfaceConservativePhase& source,const std::vector<double>& oldVelocity,
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceConservativePhase& target,
    double liquidDensity,double gasDensity)
{
    using Transport=Quantum::Hydrodynamics::HydroConservativeTransport;
    const auto require=[](bool ok,const char* why) {if(!ok)throw std::invalid_argument(std::string("[HydroMomentumRemap] ")+why);};
    const auto conservation=activeOptions_.completeConservation.value_or(HydroCompleteConservation{});
    conservation.validate();
    Transport::ProjectedMomentumRemap result;
    result.unprojected=Transport::remapMomentum(oldMesh,oldDomain,source,oldVelocity,mesh,domain,target,liquidDensity,gasDensity);
    const auto& raw=result.unprojected;const int count=domain.nativeFaceCount(),cells=mesh.cellIndexNum();
    const bool rz=domain.axisymmetric();const auto& patches=activeOptions_.boundaryPatches;
    validateHydroBoundaryPatches(patches,mesh.minCoord,mesh.maxCoord,
        {(mesh.maxCoord.x-mesh.minCoord.x)/mesh.Nx,(mesh.maxCoord.y-mesh.minCoord.y)/mesh.Ny},rz);
    const auto patchAt=[&](const NativeLeafFace& face)->const HydroBoundaryPatch* {
        const auto side=face.horizontal?(face.outwardSign<0?HydroBoundarySide::South:HydroBoundarySide::North)
            :(face.outwardSign<0?HydroBoundarySide::West:HydroBoundarySide::East);
        const auto coordinate=face.horizontal?face.center.x:face.center.y;
        for(const auto& patch:patches)if(patch.side==side && coordinate>patch.start && coordinate<patch.end)return &patch;
        return nullptr;
    };
    std::vector<HydroBoundaryProjection::Face> faces(count);bool open=false;
    std::vector<Eigen::Triplet<double>> incidence;
    Eigen::VectorXd predicted(count);long double componentMass[2]{};
    for(const auto face:domain.faces(mesh)) {
        auto& f=faces[face.nativeRow];f.area=face.area;f.inverseMass=static_cast<double>(1/raw.mass[face.nativeRow]);
        require(std::isfinite(f.inverseMass) && f.inverseMass>0,"unrepresentable inverse target mass");
        const int negative=face.negativeCell<0?-1:domain.cellRowByIndex()[face.negativeCell];
        const int positive=face.positiveCell<0?-1:domain.cellRowByIndex()[face.positiveCell];
        if(negative>=0)incidence.emplace_back(negative,face.nativeRow,f.area);
        if(positive>=0)incidence.emplace_back(positive,face.nativeRow,-f.area);
        if(face.outwardSign==0) {f.high=positive;f.low=negative;}
        else {
            const auto* patch=patchAt(face);f.high=positive>=0?positive:negative;f.low=-1;f.sign=face.outwardSign;
            f.fixed=!patch || patch->type==HydroBoundaryType::NormalVelocity;
            f.fixedVelocity=f.fixed && patch?f.sign*patch->value:0;open=open || !f.fixed;
        }
        predicted[face.nativeRow]=f.fixed?f.fixedVelocity:raw.velocity[face.nativeRow];
        componentMass[face.horizontal?1:0]+=raw.mass[face.nativeRow];
    }
    // Source boundary data are unchanged through zero-time refinement. An
    // arbitrary change of inlet/wall values cannot masquerade as remapping.
    double sourceSpeed=0;for(double u:oldVelocity)sourceSpeed=std::max(sourceSpeed,std::abs(u));
    for(const auto face:oldDomain.faces(oldMesh))if(face.outwardSign!=0) {
        const auto* patch=patchAt(face);
        if(patch && patch->type==HydroBoundaryType::PressureOutlet)continue;
        const double fixed=patch?face.outwardSign*patch->value:0;
        require(std::abs(oldVelocity[face.nativeRow]-fixed)<=512*std::numeric_limits<double>::epsilon()*sourceSpeed,
            "accepted source disagrees with its prescribed normal boundary");
    }
    Eigen::SparseMatrix<double> divergence(cells,count);divergence.setFromTriplets(incidence.begin(),incidence.end());
    const HydroBoundaryProjection projection(cells,std::move(faces),[&](std::size_t row){return predicted[row];});
    auto matrix=projection.equations().matrix;
    if(!open) {
        // The prescribed multiplier is zero. Eliminate its column as well as
        // its row to retain the symmetry of the remaining mass projection.
        matrix.prune([](int row,int column,double){return row!=0 && column!=0;});matrix.coeffRef(0,0)=1;
    }
    matrix.makeCompressed();
    // One temporary factor may be reused for the predictor, two constraint
    // directions and residual refinement. Restore the exact physical report
    // and old pressure cache on every exit, including exceptions.
    struct Restore {
        HydroSolveReport& report;HydroSolveReport saved;
        std::shared_ptr<LinearSolveCache>& cache;std::shared_ptr<LinearSolveCache> savedCache;
        Restore(HydroSolveReport& r,std::shared_ptr<LinearSolveCache>& c):report(r),cache(c)
        {std::swap(report,saved);std::swap(cache,savedCache);}
        ~Restore(){std::swap(report,saved);std::swap(cache,savedCache);}
    } restore(activeReport_,pressureLinearCache_);
    // This velocity is the next conservative phase-flux source. A pressure
    // solve accepted at the general NS tolerance can leave volume divergence
    // above that consumer's roundoff budget, especially at high density ratio.
    // Correct the actual returned velocity with the same mass projection;
    // never repair cell inventories or relax the phase transport gate.
    const double remapContinuity=std::min(conservation.initialContinuity,
        32*std::numeric_limits<double>::epsilon());
    const auto project=[&](const Eigen::VectorXd& value) {
        Eigen::VectorXd corrected=value;
        for(int pass=0;pass<4;++pass) {
            Eigen::VectorXd rhs=-(divergence*corrected);if(!open)rhs[0]=0;
            double residual=0;
            const auto multiplier=checkedLinearSolve(matrix,rhs,"pressure",residual,{},3);
            result.linearRelativeResidual=std::max(result.linearRelativeResidual,residual);
            ++result.projectionSolves;result.usedSuperLU=true;
            projection.apply(multiplier,[&](std::size_t row){return corrected[row];},
                [&](std::size_t row,double u){corrected[row]=u;},true,"nonfinite conservative remap projection");
            if(HydroContinuityBalance::LocalFlux(divergence,corrected).satisfiesScaledFlux(remapContinuity))break;
        }
        // A projected component mode can be mathematically zero: its tiny
        // floating remainder need not meet a relative gate. SVD below removes
        // dependent modes; the final physical velocity must pass the gate.
        return corrected;
    };
    const Eigen::VectorXd rawVelocity=Eigen::Map<const Eigen::VectorXd>(raw.velocity.data(),count);
    result.divergenceBefore=HydroContinuityBalance::LocalFlux(divergence,rawVelocity).maxDivergence(domain.cellVolumes());
    Eigen::VectorXd velocity=project(predicted);
    // z_j=P(e_j) spans precisely the mass-orthogonal divergence-free response
    // to a component momentum multiplier. Fixed boundary entries of e_j are 0.
    // Normalize BEFORE forming the small Gram system. SVD of sqrt(M)*Z avoids
    // subtracting nearly equal masses to discover redundant closed-box modes.
    Eigen::MatrixXd modes(count,2),weighted(count,2);double massScale[2]{};
    for(int component=0;component<2;++component) {
        massScale[component]=std::sqrt(static_cast<double>(componentMass[component]));
        Eigen::VectorXd basis=Eigen::VectorXd::Zero(count);
        for(int row=0;row<count;++row)if(!projection.face(row).fixed && (row<mesh.vEdgeIndexNum()?0:1)==component)
            basis[row]=1/massScale[component];
        modes.col(component)=project(basis);
        for(int row=0;row<count;++row)weighted(row,component)=std::sqrt(static_cast<double>(raw.mass[row]))*modes(row,component);
    }
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(weighted,Eigen::ComputeThinV);
    const double rankTolerance=4096*std::numeric_limits<double>::epsilon();
    for(int j=0;j<2;++j)if(svd.singularValues()[j]>rankTolerance)++result.independentMomentumConstraints;
    // Scale conservation by the original extensive absolute momenta. Retain
    // zero components without introducing a dimensional absolute tolerance.
    const auto constrainMomenta=[&](Eigen::VectorXd& velocity) {
    bool compatible=false;
    for(int pass=0;pass<4;++pass) {
        long double momentum[2]{};
        for(int row=0;row<count;++row)momentum[row<mesh.vEdgeIndexNum()?0:1]+=raw.mass[row]*velocity[row];
        Eigen::Vector2d wanted,coefficient=Eigen::Vector2d::Zero();
        for(int j=0;j<2;++j)wanted[j]=static_cast<double>((raw.sourceMomentum[j]-momentum[j])/massScale[j]);
        for(int j=0;j<2;++j) {
            const auto direction=svd.matrixV().col(j);const double singular=svd.singularValues()[j];
            if(singular>rankTolerance)coefficient+=direction*(direction.dot(wanted)/(singular*singular));
            // Dependent modes are never silently declared satisfied: the
            // final physical momentum check below also covers these directions.
        }
        velocity+=modes*coefficient;
        const HydroContinuityBalance::LocalFlux continuity(divergence,velocity);
        result.divergenceAfter=continuity.maxDivergence(domain.cellVolumes());
        result.localContinuity=continuity.relativeResidual();
        result.momentum[0]=result.momentum[1]=0;long double absolute[2]{};
        result.kineticEnergy=0;result.momentumRelativeClosure=0;bool momentumOK=true;
        for(int row=0;row<count;++row) {
            const int component=row<mesh.vEdgeIndexNum()?0:1;const auto p=raw.mass[row]*velocity[row];
            result.momentum[component]+=p;absolute[component]+=std::abs(p);result.kineticEnergy+=.5L*p*velocity[row];
        }
        for(int component=0;component<2;++component) {
            const auto error=std::abs(result.momentum[component]-raw.sourceMomentum[component]);
            const auto scale=std::max(raw.sourceAbsoluteMomentum[component],absolute[component]);
            const double relative=scale>0?static_cast<double>(error/scale):error==0?0:std::numeric_limits<double>::infinity();
            result.momentumRelativeClosure=std::max(result.momentumRelativeClosure,relative);
            momentumOK=momentumOK && relative<=conservation.momentumBalance+1024*std::numeric_limits<double>::epsilon();
        }
        if(momentumOK && continuity.satisfiesScaledFlux(remapContinuity)) {compatible=true;break;}
        if(pass<3)velocity=project(velocity);
    }
    require(compatible,"joint remap projection has incompatible or unresolved continuity/momentum constraints");
    };
    constrainMomenta(velocity);
    const auto sourceEnergy=raw.sourceKineticEnergy[0]+raw.sourceKineticEnergy[1];
    result.kineticEnergyBeforeConstraint=result.kineticEnergy;
    const auto energyBudget=sourceEnergy*(1+4096*std::numeric_limits<double>::epsilon());
    if(std::isfinite(result.kineticEnergy) && result.kineticEnergy>energyBudget) {
        // Intersect the affine continuity/momentum/boundary constraints with
        // the original kinetic-energy ball. The minimum-energy feasible field
        // is its centre in the mass metric. Contracting only the orthogonal
        // fluctuation is the nearest feasible remap, unlike scaling all
        // velocities (which would destroy momentum and prescribed inflow).
        Eigen::VectorXd minimum=Eigen::VectorXd::Zero(count);
        for(int row=0;row<count;++row)if(projection.face(row).fixed)
            minimum[row]=projection.face(row).fixedVelocity;
        minimum=project(minimum);constrainMomenta(minimum);
        result.minimumKineticEnergy=result.kineticEnergy;
        require(std::isfinite(result.minimumKineticEnergy) && result.minimumKineticEnergy<=sourceEnergy,
            "joint remap constraints have no field within the source kinetic-energy budget");
        const Eigen::VectorXd fluctuation=velocity-minimum;
        long double quadratic=0,linear=0;
        for(int row=0;row<count;++row) {
            quadratic+=.5L*raw.mass[row]*fluctuation[row]*fluctuation[row];
            linear+=raw.mass[row]*minimum[row]*fluctuation[row];
        }
        const auto gap=sourceEnergy-result.minimumKineticEnergy;
        const auto discriminant=linear*linear+4*quadratic*gap;
        require(quadratic>0 && std::isfinite(discriminant) && discriminant>=0,
            "invalid constrained remap energy quadratic");
        const auto root=std::sqrt(discriminant);
        const auto scale=gap==0?0.L:linear>=0?2*gap/(linear+root):(root-linear)/(2*quadratic);
        require(std::isfinite(scale) && scale>=0 && scale<1,"invalid constrained remap energy contraction");
        result.energyContraction=static_cast<double>(scale);result.energyContractionApplied=true;
        velocity=minimum+result.energyContraction*fluctuation;
        // Check/correct the returned field itself; roundoff is not permission
        // to omit either continuity or extensive momentum from acceptance.
        constrainMomenta(velocity);
    }
    if(!(std::isfinite(result.kineticEnergy) && result.kineticEnergy<=sourceEnergy*(1+4096*std::numeric_limits<double>::epsilon()))) {
        std::ostringstream message;message<<std::setprecision(17)
            <<"joint remap projection exceeds the source kinetic-energy budget; source="<<sourceEnergy
            <<", projected="<<result.kineticEnergy<<", relative_change="<<(result.kineticEnergy/sourceEnergy-1)
            <<", momentum_closure="<<result.momentumRelativeClosure<<", continuity="<<result.localContinuity;
        throw std::invalid_argument("[HydroMomentumRemap] "+message.str());
    }
    double maxSpeed=velocity.cwiseAbs().maxCoeff();
    for(int row=0;row<count;++row) {
        const auto& face=projection.face(row);
        require(!face.fixed || velocity[row]==face.fixedVelocity,"projection changed a fixed boundary velocity");
        require(face.low>=0 || face.fixed || face.sign*velocity[row]>=-(conservation.outletAdmissibility
            +128*std::numeric_limits<double>::epsilon())*std::max(sourceSpeed,maxSpeed),"remap requires unsupported outlet backflow");
        result.maxVelocityCorrection=std::max(result.maxVelocityCorrection,std::abs(velocity[row]-raw.velocity[row]));
    }
    result.velocity.assign(velocity.data(),velocity.data()+velocity.size());
    return result;
}
