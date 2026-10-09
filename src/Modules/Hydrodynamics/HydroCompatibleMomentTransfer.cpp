#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentTransfer.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <Eigen/SparseQR>
#include <Eigen/SuperLUSupport>
#include <numeric>
#include <map>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

namespace Quantum::Hydrodynamics::MomentTransfer {
namespace {
using Space=HydroCompatibleMomentSpace;using Matrix=Space::Matrix;
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
Matrix momentum(const QuadtreeMesh& mesh,const Space& space,const HydroCompatibleMomentMaterial& material,
    const std::vector<int>& patches,int count){
    std::vector<Eigen::Triplet<double>> entries;
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto h=mesh.cells(mesh.cellIndexToLabel[c]).size();const auto b=Space::affineVelocityBasis({0,0},h);
        const auto& d=material.cells()[c].density;const auto& block=space.cells()[c];
        const Eigen::MatrixXd row=(d(0,0)*b+d(0,1)*(Space::affineVelocityBasis({1,0},h)-b)+d(0,2)*(Space::affineVelocityBasis({0,1},h)-b))*block.affineProjection;
        for(int k=0;k<int(block.rows.size());k++)for(int d=0;d<2;d++)if(row(d,k)!=0)entries.emplace_back(2*patches[c]+d,block.rows[k],row(d,k));
    }Matrix result(2*count,space.mass().rows());result.setFromTriplets(entries.begin(),entries.end());return result;
}
// Optional read-only decomposition of the SAME owned kinetic norm. The
// diagnostic alternatives below are not admissibility rules or fallbacks.
void diagnosticEnergy(std::ofstream& out,const char* stage,const QuadtreeMesh& mesh,
    const Space& space,const HydroCompatibleMomentMaterial& material,const Eigen::VectorXd& value){
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto x=cell.centerCoord(),h=cell.size();
        const auto& block=space.cells()[c];const auto& data=material.cells()[c];Eigen::VectorXd u(block.rows.size());
        for(int k=0;k<u.size();k++)u[k]=value[block.rows[k]];
        const Eigen::VectorXd a=block.affineProjection*u,r=u-block.affineSamples*a;
        const auto b0=Space::affineVelocityBasis({0,0},h),bx=(Space::affineVelocityBasis({1,0},h)-b0).eval(),by=(Space::affineVelocityBasis({0,1},h)-b0).eval();
        const Eigen::Vector2d velocity[3]={b0*a,bx*a,by*a};double affine=0;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)affine+=.5*data.density(i,j)*velocity[i].dot(velocity[j]);
        const double complement=.5*data.density(0,0)/data.volume(0,0)*r.dot(block.mass*r);
        const Eigen::Vector2d p=data.density(0,0)*velocity[0]+data.density(0,1)*velocity[1]+data.density(0,2)*velocity[2];
        out<<stage<<','<<c<<','<<x.x<<','<<x.y<<','<<h.x<<','<<h.y<<','<<data.liquid(0,0)/data.volume(0,0)<<','<<data.density(0,0)
            <<','<<affine<<','<<complement<<','<<p.x()<<','<<p.y()<<'\n';
    }
}
double diagnosticMinimum(const Matrix& mass,const Matrix& constraints,const Eigen::VectorXd& rhs,
    const std::vector<int>& free,const Eigen::VectorXd& prescribed,double& residual,double& stationarity){
    const int nv=int(mass.rows()),nr=int(constraints.rows()),nf=*std::max_element(free.begin(),free.end())+1;
    Eigen::VectorXd scale(nv);for(int k=0;k<nv;k++)scale[k]=1/std::sqrt(mass.coeff(k,k));
    Eigen::MatrixXd a=Eigen::MatrixXd::Zero(nr,nf),h=Eigen::MatrixXd::Zero(nf,nf);
    Eigen::VectorXd b=rhs-constraints*prescribed,load=Eigen::VectorXd::Zero(nf),lift=mass*prescribed;
    for(int i=0;i<nv;i++)if(free[i]>=0)load[free[i]]=-scale[i]*lift[i];
    for(int col=0;col<nv;col++)if(free[col]>=0){
        for(Matrix::InnerIterator it(constraints,col);it;++it)a(it.row(),free[col])=it.value()*scale[col];
        for(Matrix::InnerIterator it(mass,col);it;++it)if(free[it.row()]>=0)h(free[it.row()],free[col])=it.value()*scale[it.row()]*scale[col];}
    for(int row=0;row<nr;row++){const double norm=a.row(row).norm();if(norm>0){a.row(row)/=norm;b[row]/=norm;}}
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a.transpose());qr.setThreshold(1e-10);const int rank=int(qr.rank());
    Eigen::MatrixXd system=Eigen::MatrixXd::Zero(nf+rank,nf+rank);system.topLeftCorner(nf,nf)=h;
    Eigen::VectorXd right=Eigen::VectorXd::Zero(nf+rank);right.head(nf)=load;
    for(int k=0;k<rank;k++){const int row=qr.colsPermutation().indices()[k];system.block(nf+k,0,1,nf)=a.row(row);system.block(0,nf+k,nf,1)=a.row(row).transpose();right[nf+k]=b[row];}
    const Eigen::VectorXd solution=system.fullPivLu().solve(right);require(solution.allFinite(),"diagnostic kinetic solve failed");
    residual=(a*solution.head(nf)-b).norm()/std::max({b.norm(),solution.head(nf).norm(),std::numeric_limits<double>::min()});
    stationarity=(system*solution-right).norm()/std::max(right.norm(),std::numeric_limits<double>::min());
    require(residual<2e-10&&stationarity<2e-10,"diagnostic minimum does not satisfy its constraints");
    Eigen::VectorXd u=prescribed;for(int i=0;i<nv;i++)if(free[i]>=0)u[i]=scale[i]*solution[free[i]];
    return .5*u.dot(mass*u);
}
}
static Eigen::VectorXd velocityImpl(const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const Space& oldSpace,
    const HydroCompatibleMomentMaterial& oldMaterial,const Eigen::VectorXd& oldVelocity,double time,
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,const HydroCompatibleMomentMaterial& material,
    const InterfaceParametricProjection::MeshTransferReport& overlay,const std::vector<int>& oldCells,const std::vector<int>& cells,
    HydroCompatibleMomentFlow::TransferLevel& report,const HydroCompatibleMomentBoundary* sourceBoundary,
    const HydroCompatibleMomentBoundary* targetBoundary,bool sameTopologyEvent,double surfaceChange,Feasibility* feasibility=nullptr){
    oldMaterial.validate(oldMesh,oldDomain,oldSpace);material.validate(mesh,domain,space);
    require(oldDomain.axisymmetric()==domain.axisymmetric(),"moment remap changed physical measure");
    const int no=oldMesh.cellIndexNum(),nn=mesh.cellIndexNum(),nv=int(space.mass().rows());
    std::vector<int> parents(no+nn);std::iota(parents.begin(),parents.end(),0);
    const auto root=[&](int i){while(parents[i]!=i){parents[i]=parents[parents[i]];i=parents[i];}return i;};
    for(const auto& o:overlay.overlaps)parents[root(oldCells.at(o.source))]=root(no+cells.at(o.target));
    std::map<int,int> ids;std::vector<int> oldPatch(no),patch(nn);
    for(int k=0;k<no+nn;k++){const int r=root(k);const auto found=ids.emplace(r,int(ids.size()));if(k<no)oldPatch[k]=found.first->second;else patch[k-no]=found.first->second;}
    const int patches=int(ids.size());
    require(oldVelocity.size()==oldSpace.mass().rows()&&oldVelocity.allFinite(),"invalid moment remap source velocity");
    if(sameTopologyEvent)require(&oldMesh==&mesh && &oldDomain==&domain && &oldSpace==&space,"event projection requires unchanged native topology");
    if(sourceBoundary)sourceBoundary->validate(oldMesh,oldDomain,oldSpace);
    if(targetBoundary)targetBoundary->validate(mesh,domain,space);
    // The strictly convex projection has the exact minimizer u=0 when its
    // prior, momentum, divergence and boundary data are all exactly zero.
    // Its minimum kinetic energy is also known: no global rank analysis or
    // factorization can add information. This is not a small-velocity cutoff.
    if(oldVelocity.isZero(0)&&(!targetBoundary||targetBoundary->values().isZero(0))){
        require(std::isfinite(surfaceChange),"nonfinite event/transfer energy bound");
        report={};report.time=time;report.patches=patches;report.independentConstraints=-1;report.kineticCeiling=-surfaceChange;
        if(surfaceChange>0){
            if(feasibility){*feasibility=Feasibility::Energy;return {};}
            throw std::invalid_argument("zero-state remap has no energy-admissible state");
        }
        if(feasibility)*feasibility=Feasibility::Admissible;
        return Eigen::VectorXd::Zero(nv);
    }
    const Matrix oldMomentum=momentum(oldMesh,oldSpace,oldMaterial,oldPatch,patches),newMomentum=momentum(mesh,space,material,patch,patches);
    const Eigen::VectorXd desired=oldMomentum*oldVelocity;
    Eigen::VectorXd prior;
    if(sameTopologyEvent)prior=oldVelocity;
    else {
    const auto spline=HydroCompatibleMomentSpline::build(oldMesh,oldDomain,oldSpace,oldVelocity,time,sourceBoundary);
    Space::Partition planes;for(int c=0;c<no;c++){const auto& cell=oldMesh.cells(oldMesh.cellIndexToLabel[c]);const auto p=cell.centerCoord(),h=cell.size();
        for(double x:{p.x-.5*h.x,p.x+.5*h.x})planes[0].push_back(x);for(double y:{p.y-.5*h.y,p.y+.5*h.y})planes[1].push_back(y);}
    for(auto& p:planes){std::sort(p.begin(),p.end());p.erase(std::unique(p.begin(),p.end()),p.end());}
    prior=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=spline.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);},planes);
    }
    std::vector<bool> boundary(nv,false);Eigen::VectorXd prescribed=Eigen::VectorXd::Zero(nv);
    if(targetBoundary){for(int i=0;i<nv;i++)if(targetBoundary->fixed()[i]){boundary[i]=true;prescribed[i]=targetBoundary->values()[i];}}
    else for(const auto f:domain.faces(mesh))if(f.outwardSign&&f.area>0)for(int k=0;k<2;k++)boundary[space.faces()[f.nativeRow].first+k]=true;
    std::vector<int> free(nv,-1);int nf=0;for(int i=0;i<nv;i++)if(!boundary[i])free[i]=nf++;else prior[i]=prescribed[i];
    const auto& mass=material.kineticMass();const int nd=int(space.divergence().rows()),nr=nd+2*patches;
    Eigen::VectorXd scale(nv);for(int i=0;i<nv;i++){require(mass.coeff(i,i)>0,"moment remap mass is not positive");scale[i]=1/std::sqrt(mass.coeff(i,i));}
    std::vector<Eigen::Triplet<double>> constraints;Eigen::VectorXd rhs=Eigen::VectorXd::Zero(nr);
    rhs.head(nd)=-space.divergence()*prescribed;rhs.tail(2*patches)=desired-newMomentum*prescribed;
    for(int col=0;col<space.divergence().outerSize();col++)if(free[col]>=0)for(Matrix::InnerIterator it(space.divergence(),col);it;++it)constraints.emplace_back(it.row(),free[col],it.value()*scale[col]);
    for(int col=0;col<nv;col++)if(free[col]>=0)for(Matrix::InnerIterator it(newMomentum,col);it;++it)constraints.emplace_back(nd+it.row(),free[col],it.value()*scale[col]);
    Matrix all(nr,nf);all.setFromTriplets(constraints.begin(),constraints.end());Eigen::VectorXd norms=Eigen::VectorXd::Zero(nr);
    for(int col=0;col<all.outerSize();col++)for(Matrix::InnerIterator it(all,col);it;++it)norms[it.row()]+=it.value()*it.value();
    for(int row=0;row<nr;row++)norms[row]=norms[row]>0?1/std::sqrt(norms[row]):1;
    all=norms.asDiagonal()*all;rhs=norms.asDiagonal()*rhs;
    Eigen::SparseQR<Matrix,Eigen::COLAMDOrdering<int>> qr;qr.setPivotThreshold(1e-10);qr.compute(Matrix(all.transpose()));
    require(qr.info()==Eigen::Success,"moment remap constraint rank failed");const int rank=int(qr.rank());
    std::vector<int> selected(nr,-1);for(int k=0;k<rank;k++)selected[qr.colsPermutation().indices()[k]]=k;
    std::vector<Eigen::Triplet<double>> entries;Eigen::VectorXd load=Eigen::VectorXd::Zero(nf+rank),oldLoad=mass*(prior-prescribed);
    for(int i=0;i<nv;i++)if(free[i]>=0)load[free[i]]=scale[i]*oldLoad[i];
    for(int col=0;col<mass.outerSize();col++)if(free[col]>=0)for(Matrix::InnerIterator it(mass,col);it;++it)if(free[it.row()]>=0)
        entries.emplace_back(free[it.row()],free[col],scale[it.row()]*it.value()*scale[col]);
    for(int row=0;row<nr;row++)if(selected[row]>=0)load[nf+selected[row]]=rhs[row];
    for(int col=0;col<all.outerSize();col++)for(Matrix::InnerIterator it(all,col);it;++it)if(selected[it.row()]>=0){entries.emplace_back(nf+selected[it.row()],col,it.value());entries.emplace_back(col,nf+selected[it.row()],it.value());}
    Matrix system(nf+rank,nf+rank);system.setFromTriplets(entries.begin(),entries.end());Eigen::SuperLU<Matrix> solver;solver.compute(system);
    require(solver.info()==Eigen::Success,"moment remap constrained mass factorization failed");const Eigen::VectorXd solution=solver.solve(load);
    require(solver.info()==Eigen::Success&&solution.allFinite(),"moment remap constrained mass solve failed");
    const double residual=(all*solution.head(nf)-rhs).norm();
    const double residualScale=std::max(rhs.norm(),solution.head(nf).norm());
    report.time=time;report.patches=patches;report.independentConstraints=rank;
    report.kineticBefore=.5*oldVelocity.dot(oldMaterial.kineticMass()*oldVelocity);report.kineticCeiling=report.kineticBefore-surfaceChange;
    report.constraintResidual=residualScale>0?residual/residualScale:0;
    if(residual>2e-10*residualScale){
        if(feasibility){*feasibility=Feasibility::Momentum;return {};}
        std::ostringstream message;message.precision(17);message<<"moment remap constraints incompatible: residual="<<residual<<" rhs="<<rhs.norm()<<" rank="<<rank;throw std::invalid_argument(message.str());}
    Eigen::VectorXd result=prescribed;for(int i=0;i<nv;i++)if(free[i]>=0)result[i]=scale[i]*solution[free[i]];
    // Projection onto the SAME affine conservation space intersected with
    // the old kinetic-energy ball. The affine minimum-energy state is
    // orthogonal to its nullspace. Thus the unique constrained minimizer is
    // obtained along that nullspace radius, not by scaling physical momentum.
    Eigen::VectorXd minimumLoad=load;minimumLoad.head(nf).setZero();const Eigen::VectorXd liftedMass=mass*prescribed;
    for(int i=0;i<nv;i++)if(free[i]>=0)minimumLoad[free[i]]=-scale[i]*liftedMass[i];
    const Eigen::VectorXd minimumSolution=solver.solve(minimumLoad);
    require(solver.info()==Eigen::Success&&minimumSolution.allFinite(),"moment remap minimum-energy solve failed");
    Eigen::VectorXd minimum=prescribed;for(int i=0;i<nv;i++)if(free[i]>=0)minimum[i]=scale[i]*minimumSolution[free[i]];
    report.kineticBefore=.5*oldVelocity.dot(oldMaterial.kineticMass()*oldVelocity);
    report.minimumKinetic=.5*minimum.dot(mass*minimum);report.unconstrainedKinetic=.5*result.dot(mass*result);
    if(const char* prefix=std::getenv("QUANTUM_MOMENT_TRANSFER_DIAGNOSTICS");prefix&&report.minimumKinetic>report.kineticCeiling
        +4096*std::numeric_limits<double>::epsilon()*std::max({report.kineticBefore,report.minimumKinetic,std::abs(surfaceChange)})){
        const std::string base=std::string(prefix)+(domain.axisymmetric()?"_rz_":"_xy_")+std::to_string(no)+"_"+std::to_string(nn)+"_"+std::to_string(time);
        std::ofstream cellsOut(base+"_cells.csv"),summary(base+"_summary.csv");require(cellsOut&&summary,"cannot open transfer diagnostics");
        cellsOut<<std::setprecision(17)<<"stage,cell,x,y,hx,hy,liquid_fraction,mass,affine_energy,complement_energy,momentum_x,momentum_y\n";
        diagnosticEnergy(cellsOut,"source",oldMesh,oldSpace,oldMaterial,oldVelocity);
        diagnosticEnergy(cellsOut,"minimum",mesh,space,material,minimum);
        diagnosticEnergy(cellsOut,"closest",mesh,space,material,result);
        const double minimumConstraint=(all*minimumSolution.head(nf)-rhs).norm()/std::max({rhs.norm(),minimumSolution.head(nf).norm(),std::numeric_limits<double>::min()});
        const double stationarity=(system*minimumSolution-minimumLoad).norm()/std::max(minimumLoad.norm(),std::numeric_limits<double>::min());
        summary<<std::setprecision(17)<<"constraints,kinetic,relative_constraint_residual,relative_kkt_residual\nsource,"<<report.kineticBefore<<",,\npatch_momentum_divergence,"<<report.minimumKinetic<<','<<minimumConstraint<<','<<stationarity<<'\n';
        double diagnosticResidual=0,diagnosticStationarity=0;
        const auto recordMinimum=[&](const char* name,const Matrix& a,const Eigen::VectorXd& b){
            const double energy=diagnosticMinimum(mass,a,b,free,prescribed,diagnosticResidual,diagnosticStationarity);summary<<name<<','<<energy<<','<<diagnosticResidual<<','<<diagnosticStationarity<<'\n';};
        recordMinimum("divergence_only",space.divergence(),Eigen::VectorXd::Zero(nd));
        recordMinimum("patch_momentum_only",newMomentum,desired);
        std::vector<Eigen::Triplet<double>> globalEntries;
        for(int col=0;col<space.divergence().outerSize();col++)for(Matrix::InnerIterator it(space.divergence(),col);it;++it)globalEntries.emplace_back(it.row(),col,it.value());
        for(int col=0;col<newMomentum.outerSize();col++)for(Matrix::InnerIterator it(newMomentum,col);it;++it)globalEntries.emplace_back(nd+it.row()%2,col,it.value());
        Matrix global(nd+2,nv);global.setFromTriplets(globalEntries.begin(),globalEntries.end());Eigen::VectorXd globalRight=Eigen::VectorXd::Zero(nd+2);
        for(int k=0;k<patches;k++)globalRight.tail<2>()+=desired.segment<2>(2*k);
        recordMinimum("global_momentum_divergence",global,globalRight);
        cellsOut.flush();summary.flush();require(cellsOut&&summary,"cannot write transfer diagnostics");
    }
    const double ceiling=report.kineticBefore-surfaceChange;
    require(std::isfinite(ceiling),"nonfinite event/transfer energy bound");
    const double energyRoundoff=4096*std::numeric_limits<double>::epsilon()*std::max({report.kineticBefore,report.minimumKinetic,std::abs(surfaceChange)});
    if(report.minimumKinetic>ceiling+energyRoundoff){
        if(feasibility){*feasibility=Feasibility::Energy;return {};}
        std::ostringstream message;message.precision(17);message<<"moment remap energy/conservation infeasible: before="<<report.kineticBefore<<" surface_change="<<surfaceChange<<" minimum="<<report.minimumKinetic<<" unconstrained="<<report.unconstrainedKinetic;throw std::invalid_argument(message.str());}
    if(report.unconstrainedKinetic>ceiling+energyRoundoff){
        const Eigen::VectorXd direction=result-minimum;const double radius=.5*direction.dot(mass*direction);
        require(radius>0,"moment remap has no energy-admissible nullspace direction");
        report.retainedNullspace=std::sqrt(std::max(0.,ceiling-report.minimumKinetic)/radius);
        require(report.retainedNullspace<=1+4096*std::numeric_limits<double>::epsilon(),"moment remap kinetic decomposition inconsistent");
        result=minimum+report.retainedNullspace*direction;report.energyLimited=true;
    }
    const double divScale=(space.divergence().cwiseAbs()*result.cwiseAbs()).norm(),momentumScale=(oldMomentum.cwiseAbs()*oldVelocity.cwiseAbs()).norm();
    report.time=time;report.kineticBefore=.5*oldVelocity.dot(oldMaterial.kineticMass()*oldVelocity);report.kineticAfter=.5*result.dot(mass*result);
    report.divergenceResidual=divScale>0?(space.divergence()*result).norm()/divScale:0;
    report.momentumResidual=momentumScale>0?(newMomentum*result-desired).norm()/momentumScale:0;
    report.correctionNorm=std::sqrt(std::max(0.,(result-prior).dot(mass*(result-prior))));report.patches=patches;report.independentConstraints=rank;
    require(report.kineticAfter<=ceiling+energyRoundoff,"moment remap exceeded its declared energy bound");
    require(report.divergenceResidual<=4096*std::numeric_limits<double>::epsilon()&&report.momentumResidual<2e-10,"moment remap conservation residual rejected");
    if(feasibility)*feasibility=Feasibility::Admissible;return result;
}
Eigen::VectorXd velocity(const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const Space& oldSpace,
    const HydroCompatibleMomentMaterial& oldMaterial,const Eigen::VectorXd& oldVelocity,double time,
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,const HydroCompatibleMomentMaterial& material,
    const InterfaceParametricProjection::MeshTransferReport& overlay,const std::vector<int>& oldCells,const std::vector<int>& cells,
    HydroCompatibleMomentFlow::TransferLevel& report,const HydroCompatibleMomentBoundary* sourceBoundary,const HydroCompatibleMomentBoundary* targetBoundary){
    return velocityImpl(oldMesh,oldDomain,oldSpace,oldMaterial,oldVelocity,time,mesh,domain,space,material,overlay,oldCells,cells,report,sourceBoundary,targetBoundary,false,0);
}
Eigen::VectorXd eventVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,
    const HydroCompatibleMomentMaterial& before,const HydroCompatibleMomentMaterial& after,const Eigen::VectorXd& value,
    double time,const HydroCompatibleMomentBoundary& boundary,double surfaceChange,const std::vector<int>& cells,
    HydroCompatibleMomentFlow::TransferLevel& report,Feasibility& feasibility){
    InterfaceParametricProjection::MeshTransferReport overlay;
    for(std::size_t k=0;k<cells.size();k++)overlay.overlaps.push_back({int(k),int(k),{},0});
    return velocityImpl(mesh,domain,space,before,value,time,mesh,domain,space,after,overlay,cells,cells,report,&boundary,&boundary,true,surfaceChange,&feasibility);
}
Eigen::VectorXd pressure(const QuadtreeMesh& oldMesh,const Eigen::VectorXd& old,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceParametricProjection::MeshTransferReport& overlay,const std::vector<int>& oldCells,const std::vector<int>& cells){
    require(old.size()==3*oldMesh.cellIndexNum(),"invalid remap pressure history");
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(3*mesh.cellIndexNum());const auto& g=Geometry::Parametric::gaussEight();
    for(const auto& o:overlay.overlaps){const int from=oldCells.at(o.source),to=cells.at(o.target);const auto& a=oldMesh.cells(oldMesh.cellIndexToLabel[from]);const auto& b=mesh.cells(mesh.cellIndexToLabel[to]);
        const auto ac=a.centerCoord(),ah=a.size(),bc=b.centerCoord(),bh=b.size();const auto box=o.box;
        for(int j=0;j<8;j++)for(int i=0;i<8;i++){
            const Vector2D p{.5*(box.left+box.right)+.5*(box.right-box.left)*g.x[i],.5*(box.bottom+box.top)+.5*(box.top-box.bottom)*g.x[j]};
            const double dw=.25*(box.right-box.left)*(box.top-box.bottom)*g.w[i]*g.w[j]*(domain.axisymmetric()?2*std::acos(-1.)*p.x:1.);
            rhs.segment<3>(3*to)+=dw*Eigen::Vector3d(1,(p.x-bc.x)/bh.x,(p.y-bc.y)/bh.y)*old.segment<3>(3*from).dot(Eigen::Vector3d(1,(p.x-ac.x)/ah.x,(p.y-ac.y)/ah.y));
        }
    }
    for(int c=0;c<mesh.cellIndexNum();c++){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        Eigen::Matrix3d gram=Eigen::Matrix3d::Zero();gram.diagonal()<<volume,volume/12,volume/12;if(domain.axisymmetric())gram(0,1)=gram(1,0)=volume*cell.size().x/(12*cell.centerCoord().x);
        rhs.segment<3>(3*c)=gram.ldlt().solve(rhs.segment<3>(3*c)).eval();}
    return rhs;
}
}
