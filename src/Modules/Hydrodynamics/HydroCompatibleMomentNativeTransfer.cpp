#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentTransfer.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include <Eigen/SparseQR>
#include <Eigen/SuperLUSupport>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

namespace Quantum::Hydrodynamics::MomentTransfer {
namespace {
using Matrix=HydroCompatibleMomentSpace::Matrix;
void require(bool valid,const char* message){if(!valid)throw std::invalid_argument(message);}
Matrix gatherMomentum(const Matrix& cells,const std::vector<int>& patch,int count){
    std::vector<Eigen::Triplet<double>> entries;
    for(int col=0;col<cells.outerSize();col++)for(Matrix::InnerIterator it(cells,col);it;++it)
        entries.emplace_back(2*patch.at(it.row()/2)+it.row()%2,col,it.value());
    Matrix result(2*count,cells.cols());result.setFromTriplets(entries.begin(),entries.end());return result;
}
}
std::optional<HydroCompatibleMomentSpline::NativeVelocity> nativeVelocity(
    const HydroCompatibleMomentSpline::NativeVelocity& source,const HydroCompatibleMomentSpline::NativeMaterialForms& before,
    const HydroCompatibleMomentSpline::NativeVelocity& target,const HydroCompatibleMomentSpline::NativeProjectionForms& projection,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentBoundary& boundary,
    const InterfaceParametricProjection::MeshTransferReport& overlay,const std::vector<int>& sourceCells,const std::vector<int>& targetCells,
    double surfaceChange,HydroCompatibleMomentFlow::TransferLevel& report,Feasibility& feasibility){
    require(source.time()==target.time()&&std::isfinite(surfaceChange),"complete transfer changed time or energy definition");
    const int no=int(before.cells.size()),nn=int(projection.material.cells.size()),nq=int(target.coordinates().size());
    require(int(sourceCells.size())==no&&int(targetCells.size())==nn&&projection.load.size()==nq,"complete transfer layout mismatch");
    const Matrix mass=projection.material.mass.topLeftCorner(nq,nq);
    const auto oldBalance=source.cellBalance(before);
    const Matrix oldCellMomentum=source.cellMomentumMap(before),newCellMomentum=target.cellMomentumMap(projection.material);
    std::vector<int> parents(no+nn);std::iota(parents.begin(),parents.end(),0);
    const auto root=[&](int k){while(parents[k]!=k){parents[k]=parents[parents[k]];k=parents[k];}return k;};
    std::vector<bool> oldCovered(no,false),newCovered(nn,false);
    for(const auto& o:overlay.overlaps){const int a=sourceCells.at(o.source),b=targetCells.at(o.target);
        require(a>=0&&a<no&&b>=0&&b<nn,"invalid physical transfer overlap");
        parents[root(a)]=root(no+b);oldCovered[a]=newCovered[b]=true;}
    require(std::all_of(oldCovered.begin(),oldCovered.end(),[](bool b){return b;})&&std::all_of(newCovered.begin(),newCovered.end(),[](bool b){return b;}),"incomplete physical transfer overlay");
    std::map<int,int> ids;std::vector<int> oldPatch(no),patch(nn);
    for(int i=0;i<no+nn;i++){const int group=ids.emplace(root(i),int(ids.size())).first->second;if(i<no)oldPatch[i]=group;else patch[i-no]=group;}
    const int patches=int(ids.size());const Matrix oldMomentum=gatherMomentum(oldCellMomentum,oldPatch,patches),momentum=gatherMomentum(newCellMomentum,patch,patches);
    const Eigen::VectorXd desired=oldMomentum*source.coordinates();
    report={};report.time=source.time();report.patches=patches;report.kineticBefore=oldBalance.kinetic.sum();report.kineticCeiling=report.kineticBefore-surfaceChange;
    require(std::isfinite(report.kineticCeiling)&&projection.load.allFinite()&&std::isfinite(projection.sourceNormSquared),"nonfinite complete transfer data");
    // Exact zero is the optimizer, not a velocity cutoff or alternate model.
    if(source.coordinates().isZero(0)&&boundary.normalConstraintValues().isZero(0)){
        report.independentConstraints=-1;
        if(surfaceChange>0){feasibility=Feasibility::Energy;return std::nullopt;}
        feasibility=Feasibility::Admissible;return target.withCoordinates(Eigen::VectorXd::Zero(nq),source.time());
    }
    const Matrix h=target.kinematicConstraints(),normal=target.normalBoundaryConstraints(boundary);
    const int nh=int(h.rows()),nb=int(normal.rows()),nr=nh+nb+2*patches;
    Eigen::VectorXd scale(nq);for(int k=0;k<nq;k++){require(mass.coeff(k,k)>0,"complete transfer mass has a nonpositive diagonal");scale[k]=1/std::sqrt(mass.coeff(k,k));}
    std::vector<Eigen::Triplet<double>> entries;
    const auto append=[&](const Matrix& a,int first){for(int col=0;col<a.outerSize();col++)for(Matrix::InnerIterator it(a,col);it;++it)entries.emplace_back(first+it.row(),col,it.value()*scale[col]);};
    append(h,0);append(normal,nh);append(momentum,nh+nb);
    Matrix all(nr,nq);all.setFromTriplets(entries.begin(),entries.end());Eigen::VectorXd rhs=Eigen::VectorXd::Zero(nr);
    rhs.segment(nh,nb)=boundary.normalConstraintValues();rhs.tail(2*patches)=desired;
    Eigen::VectorXd rowScale=Eigen::VectorXd::Zero(nr);
    for(int col=0;col<all.outerSize();col++)for(Matrix::InnerIterator it(all,col);it;++it)rowScale[it.row()]+=it.value()*it.value();
    for(int k=0;k<nr;k++)rowScale[k]=rowScale[k]>0?1/std::sqrt(rowScale[k]):1;
    all=rowScale.asDiagonal()*all;rhs=rowScale.asDiagonal()*rhs;
    Eigen::SparseQR<Matrix,Eigen::COLAMDOrdering<int>> qr;qr.setPivotThreshold(1e-10);qr.compute(Matrix(all.transpose()));
    require(qr.info()==Eigen::Success,"complete transfer constraint rank failed");const int rank=int(qr.rank());report.independentConstraints=rank;
    std::vector<int> selected(nr,-1);for(int k=0;k<rank;k++)selected[qr.colsPermutation().indices()[k]]=k;
    entries.clear();
    for(int col=0;col<nq;col++)for(Matrix::InnerIterator it(mass,col);it;++it)entries.emplace_back(it.row(),col,it.value()*scale[it.row()]*scale[col]);
    for(int col=0;col<nq;col++)for(Matrix::InnerIterator it(all,col);it;++it)if(selected[it.row()]>=0){entries.emplace_back(nq+selected[it.row()],col,it.value());entries.emplace_back(col,nq+selected[it.row()],it.value());}
    Matrix system(nq+rank,nq+rank);system.setFromTriplets(entries.begin(),entries.end());Eigen::SuperLU<Matrix> solver;solver.compute(system);
    require(solver.info()==Eigen::Success,"complete transfer constrained mass factorization failed");
    Eigen::MatrixXd loads=Eigen::MatrixXd::Zero(nq+rank,2);loads.col(0).head(nq)=scale.asDiagonal()*projection.load;
    for(int row=0;row<nr;row++)if(selected[row]>=0)loads.row(nq+selected[row]).setConstant(rhs[row]);
    // One factorization, two right-hand sides: the closest physical field
    // and the affine conservation space's minimum-energy center.
    const Eigen::MatrixXd solution=solver.solve(loads);
    require(solver.info()==Eigen::Success&&solution.allFinite(),"complete transfer constrained solve failed");
    for(int col=0;col<2;col++){
        const double denominator=std::max({rhs.norm(),solution.col(col).head(nq).norm(),std::numeric_limits<double>::min()});
        report.constraintResidual=std::max(report.constraintResidual,(all*solution.col(col).head(nq)-rhs).norm()/denominator);
        const double solveScale=std::max(loads.col(col).norm(),std::numeric_limits<double>::min());
        require((system*solution.col(col)-loads.col(col)).norm()<=2e-9*solveScale,"complete transfer stationarity residual rejected");
    }
    if(report.constraintResidual>2e-10){feasibility=Feasibility::Momentum;return std::nullopt;}
    Eigen::VectorXd closest=scale.asDiagonal()*solution.col(0).head(nq),minimum=scale.asDiagonal()*solution.col(1).head(nq);
    report.minimumKinetic=.5*minimum.dot(mass*minimum);report.unconstrainedKinetic=.5*closest.dot(mass*closest);
    const double roundoff=4096*std::numeric_limits<double>::epsilon()*std::max({report.kineticBefore,report.minimumKinetic,std::abs(surfaceChange)});
    if(report.minimumKinetic>report.kineticCeiling+roundoff){feasibility=Feasibility::Energy;return std::nullopt;}
    if(report.unconstrainedKinetic>report.kineticCeiling+roundoff){
        const Eigen::VectorXd direction=closest-minimum;const double radius=.5*direction.dot(mass*direction);
        require(radius>0,"complete transfer has no admissible nullspace radius");
        report.retainedNullspace=std::sqrt(std::max(0.,report.kineticCeiling-report.minimumKinetic)/radius);
        require(report.retainedNullspace<=1+4096*std::numeric_limits<double>::epsilon(),"complete transfer energy decomposition inconsistent");
        closest=minimum+report.retainedNullspace*direction;report.energyLimited=true;
    }
    report.kineticAfter=.5*closest.dot(mass*closest);
    const Eigen::VectorXd scaled=closest.cwiseQuotient(scale);
    report.constraintResidual=std::max(report.constraintResidual,(all*scaled-rhs).norm()/std::max({rhs.norm(),scaled.norm(),std::numeric_limits<double>::min()}));
    const double momentumScale=std::max((oldMomentum.cwiseAbs()*source.coordinates().cwiseAbs()).norm(),std::numeric_limits<double>::min());
    report.momentumResidual=(momentum*closest-desired).norm()/momentumScale;
    auto result=target.withCoordinates(closest,source.time());const auto& d=space.divergence();
    const double divScale=(d.cwiseAbs()*result.nativeMoments().cwiseAbs()).norm();
    report.divergenceResidual=(d*result.nativeMoments()).norm()/std::max(divScale,std::numeric_limits<double>::min());
    const double cross=closest.dot(projection.load);
    const double distance=2*report.kineticAfter-2*cross+projection.sourceNormSquared;
    const double distanceRoundoff=4096*std::numeric_limits<double>::epsilon()
        *(2*std::abs(report.kineticAfter)+2*std::abs(cross)+std::abs(projection.sourceNormSquared));
    require(distance>=-distanceRoundoff&&projection.sourceNormSquared>0,"complete projection distance has inconsistent material integrals");
    report.correctionNorm=std::sqrt(std::max(0.,distance));
    report.relativeProjectionError=std::sqrt((std::max(0.,distance)+distanceRoundoff)/projection.sourceNormSquared);
    require(report.constraintResidual<=2e-10&&report.momentumResidual<=2e-10&&report.divergenceResidual<=2e-10,"complete transfer physical conservation rejected");
    require(report.kineticAfter<=report.kineticCeiling+roundoff,"complete transfer exceeded its physical energy bound");
    feasibility=Feasibility::Admissible;return result;
}
}
