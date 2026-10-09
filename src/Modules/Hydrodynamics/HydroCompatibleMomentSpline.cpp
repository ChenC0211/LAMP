#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMassCurrent.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineBasis.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineNative.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <Eigen/SuperLUSupport>
#include <map>
#include <mutex>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>

namespace Quantum::Hydrodynamics {
namespace {
using namespace MomentSplineBasis;
void require(bool ok,const char* reason){if(!ok)throw std::invalid_argument(reason);}
// Cell-local powers of a C2 quartic, interpolating nodal values and cell
// means, with zero endpoint first derivatives. Exactly 5n constraints for
// 5n coefficients; no smoothing penalty or tunable fit parameter.
Eigen::MatrixXd quartic(int n){
    Eigen::MatrixXd a=Eigen::MatrixXd::Zero(5*n,5*n),b=Eigen::MatrixXd::Zero(5*n,2*n+1);int row=0;
    for(int i=0;i<n;i++){
        a(row,5*i)=1;b(row++,i)=1;
        for(int k=0;k<5;k++)a(row,5*i+k)=1;b(row++,i+1)=1;
        for(int k=0;k<5;k++)a(row,5*i+k)=1./(k+1);b(row++,n+1+i)=1;
    }
    for(int i=0;i+1<n;i++){
        for(int k=1;k<5;k++)a(row,5*i+k)=k;a(row++,5*(i+1)+1)=-1;
        for(int k=2;k<5;k++)a(row,5*i+k)=k*(k-1);a(row++,5*(i+1)+2)=-2;
    }
    a(row++,1)=1;for(int k=1;k<5;k++)a(row,5*(n-1)+k)=k;row++;
    require(row==5*n,"quartic streamfunction constraint count differs");
    Eigen::FullPivLU<Eigen::MatrixXd> lu(a);require(lu.rank()==5*n,"quartic streamfunction constraints lost rank");
    Eigen::MatrixXd result=lu.solve(b);require((a*result-b).norm()<2e-12*b.norm(),"quartic streamfunction interpolation residual");return result;
}
// Tensor quintic Hermite coordinates: value, first and second nodal jets
// are shared, so chi is C2. No-slip imposes value/first jet at a wall.
// At the axis chi_r=chi_rrr=0 instead; eliminate the first-cell second
// derivative using its exact quintic third-derivative identity.
Eigen::MatrixXd quintic(int n,bool axis){
    Eigen::Matrix<double,6,6> endpoint;
    for(int d=0;d<3;d++)for(int k=0;k<6;k++){endpoint(d,k)=basis(k,0,d);endpoint(3+d,k)=basis(k,1,d);}
    const Eigen::Matrix<double,6,6> hermite=endpoint.inverse();
    std::vector<int> index(3*(n+1),-1);int count=0;
    for(int k=0;k<int(index.size());k++)if(k!=1&&k!=3*n&&k!=3*n+1&&k!=(axis?2:0))index[k]=count++;
    Eigen::MatrixXd jets=Eigen::MatrixXd::Zero(index.size(),count);
    for(int k=0;k<int(index.size());k++)if(index[k]>=0)jets(k,index[k])=1;
    if(axis)jets.row(2)=(20./3)*(jets.row(3)-jets.row(0))-(8./3)*jets.row(4)+(1./3)*jets.row(5);
    Eigen::MatrixXd result(6*n,count);
    for(int i=0;i<n;i++)result.middleRows(6*i,6)=hermite*jets.middleRows(3*i,6);
    return result;
}
using Row=std::map<int,double>;
void add(Row& a,const Row& b,double scale){for(const auto& term:b)a[term.first]+=scale*term.second;}
}
struct HydroCompatibleMomentSpline::Layout {
    std::shared_ptr<const NativeMeshIdentity> identity;
    std::shared_ptr<const NativeLayout> native;
    int nx=0,ny=0,nvx=0,nvy=0,velocityCount=0,width=5,kineticDofs=0;
    bool rz=false,axis=false;
    Vector2D lo,hi,h;
    Eigen::MatrixXd x,y;
    HydroCompatibleMomentSpace::Matrix data;
    HydroCompatibleMomentSpace::Matrix polynomial,kineticSystem;
    Eigen::SuperLU<HydroCompatibleMomentSpace::Matrix> kineticSolver;
    mutable std::mutex kineticMutex;
    struct Face {InterfaceNormalFluxView::FaceChart chart;int i,j;bool horizontal,boundary;};
    std::vector<Face> faces;
    std::vector<int> cells,boundaryRows;
    std::vector<HydroCompatibleMomentBoundary::Kind> boundaryKinds;
    struct Location {int i,j;double x,y;};
    Location locate(Vector2D p) const {
        require(std::isfinite(p.x)&&std::isfinite(p.y)&&p.x>=lo.x&&p.x<=hi.x&&p.y>=lo.y&&p.y<=hi.y,"moment spline query outside physical domain");
        const double a=(p.x-lo.x)/h.x,b=(p.y-lo.y)/h.y;const int i=std::min(nx-1,int(a)),j=std::min(ny-1,int(b));return {i,j,a-i,b-j};
    }
    Eigen::VectorXd solveKinetic(const Eigen::VectorXd& rhs) const {
        std::lock_guard<std::mutex> lock(kineticMutex);
        Eigen::VectorXd result=kineticSolver.solve(rhs);
        require(kineticSolver.info()==Eigen::Success&&result.allFinite(),"RZ minimum-energy C1 solve failed");
        require((kineticSystem*result-rhs).norm()<=2e-9*rhs.norm(),"RZ minimum-energy C1 residual unresolved");return result;
    }
    void buildKinetic(const HydroCompatibleMomentSpace& space,const std::vector<Row>& moments){
        width=6;x=quintic(nx,axis);y=quintic(ny,false);const int nodalDofs=int(x.cols()*y.cols());kineticDofs=nodalDofs+4*nx*ny;
        std::vector<Eigen::Triplet<double>> energy,polynomialEntries,constraints,targets;
        int constraintCount=0;const auto& g=Geometry::Parametric::gaussEight();const GaussTwelve radialRule;
        for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
            std::vector<int> xx,yy,ids;
            for(int k=0;k<x.cols();k++)if(!x.block(6*i,k,6,1).isZero(0))xx.push_back(k);
            for(int k=0;k<y.cols();k++)if(!y.block(6*j,k,6,1).isZero(0))yy.push_back(k);
            Eigen::MatrixXd lift=Eigen::MatrixXd::Zero(40,xx.size()*yy.size()+4);int column=0;
            for(int b:yy)for(int a:xx){ids.push_back(a+int(x.cols())*b);
                for(int ar=0;ar<6;ar++)for(int az=0;az<6;az++)lift(6*ar+az,column)=x(6*i+ar,a)*y(6*j+az,b);column++;}
            for(int k=0;k<4;k++){ids.push_back(nodalDofs+4*(j*nx+i)+k);lift(36+k,column++)=1;}
            for(int ar=0;ar<6;ar++)for(int az=0;az<6;az++)for(int k=0;k<column;k++)if(lift(6*ar+az,k)!=0)
                polynomialEntries.emplace_back(6*i+ar+6*nx*(6*j+az),ids[k],lift(6*ar+az,k));
            for(int k=0;k<4;k++)polynomialEntries.emplace_back(36*nx*ny+4*(j*nx+i)+k,nodalDofs+4*(j*nx+i)+k,1);
            const double rc=lo.x+(i+.5)*h.x,s0=(lo.x+i*h.x)/rc,d=h.x/rc;
            Eigen::Matrix<double,40,40> gram=Eigen::Matrix<double,40,40>::Zero();
            Eigen::Matrix<double,4,40> integrals=Eigen::Matrix<double,4,40>::Zero();
            for(int a=0;a<12;a++)for(int b=0;b<8;b++){
                const double qr=.5+.5*radialRule.x[a],qz=.5+.5*g.x[b],r=lo.x+(i+qr)*h.x,w=.25*radialRule.w[a]*g.w[b]*h.x*h.y*CartesianMeasure::twoPi*r;
                Eigen::Matrix<double,2,40> velocity;
                for(int ar=0;ar<6;ar++)for(int az=0;az<6;az++){
                    const double value=basis(ar,qr,0)*basis(az,qz,0);
                    velocity(0,6*ar+az)=r/h.y*basis(ar,qr,0)*basis(az,qz,1);
                    velocity(1,6*ar+az)=-2*value-r/h.x*basis(ar,qr,1)*basis(az,qz,0);}
                const auto bubbles=interiorJets(qr,qz,rc,h.x,axis&&i==0);
                for(int k=0;k<4;k++){velocity(0,36+k)=r/h.y*bubbles[k][2];velocity(1,36+k)=-2*bubbles[k][0]-r/h.x*bubbles[k][1];}
                gram+=w*velocity.transpose()*velocity;
            }
            for(int panel=0;panel<4;panel++)for(int a=0;a<8;a++)for(int b=0;b<8;b++){
                const double qr=(panel+.5+.5*g.x[a])/4,qz=.5+.5*g.x[b],w=.25*g.w[a]*g.w[b]/4,X=qr-.5,Y=qz-.5,s=s0+d*qr;
                for(int ar=0;ar<6;ar++)for(int az=0;az<6;az++){
                    const double value=basis(ar,qr,0)*basis(az,qz,0),increment=(ar==0&&az==0)?d*qr*(2*s0+d*qr):s*s*value;
                    integrals(0,6*ar+az)+=w*s*s*value;
                    integrals(1,6*ar+az)+=w*increment*X/s;
                    integrals(2,6*ar+az)+=w*increment*X*X/(s*s);
                    integrals(3,6*ar+az)+=w*increment*Y/(s*s);
                }
                const auto bubbles=interiorJets(qr,qz,rc,h.x,axis&&i==0);
                for(int k=0;k<4;k++){const double bubble=bubbles[k][0];integrals(0,36+k)+=w*s*s*bubble;integrals(1,36+k)+=w*s*X*bubble;integrals(2,36+k)+=w*X*X*bubble;integrals(3,36+k)+=w*Y*bubble;}
            }
            Eigen::MatrixXd local=lift.transpose()*gram*lift;local=.5*(local+local.transpose()).eval();
            for(int a=0;a<column;a++)for(int b=0;b<column;b++)if(local(a,b)!=0)energy.emplace_back(ids[a],ids[b],local(a,b));
            const auto append=[&](const Eigen::RowVectorXd& functional,const Row& target){
                const Eigen::RowVectorXd row=functional*lift;const double scale=row.norm();require(scale>0,"RZ C1 moment constraint is empty");
                for(int k=0;k<column;k++)if(row[k]!=0)constraints.emplace_back(constraintCount,ids[k],row[k]/scale);
                for(const auto& entry:target)if(entry.second!=0)targets.emplace_back(constraintCount,entry.first,entry.second/scale);
                constraintCount++;
            };
            const auto source=[&](int a,int b)->const Row&{return moments[b*nvx+a-int(axis)];};
            Eigen::RowVectorXd value=Eigen::RowVectorXd::Zero(40);
            if(i>0&&j>0){value[0]=1;append(value,source(i,j));value.setZero();}
            if(i>0){for(int b=0;b<6;b++)value[b]=1./(b+1);append(value,source(i,ny+1+j));value.setZero();}
            if(j>0){for(int a=0;a<6;a++)value[6*a]=s0*s0/(a+1)+2*s0*d/(a+2)+d*d/(a+3);append(value,source(nx+1+i,j));}
            append(integrals.row(0),source(nx+1+i,ny+1+j));
            for(int k=0;k<3;k++){Row target;target[space.firstCellVelocity(cells[j*nx+i])+3+k]=std::sqrt(h.x*h.y)/rc;append(integrals.row(k+1),target);}
        }
        require(constraintCount==(2*nx-1)*(2*ny-1)+3*nx*ny,"RZ C1 independent constraint count differs");
        polynomial.resize(40*nx*ny,kineticDofs);polynomial.setFromTriplets(polynomialEntries.begin(),polynomialEntries.end());
        data.resize(constraintCount,velocityCount);data.setFromTriplets(targets.begin(),targets.end());
        for(const auto& v:constraints){energy.emplace_back(kineticDofs+v.row(),v.col(),v.value());energy.emplace_back(v.col(),kineticDofs+v.row(),v.value());}
        kineticSystem.resize(kineticDofs+constraintCount,kineticDofs+constraintCount);kineticSystem.setFromTriplets(energy.begin(),energy.end());
        kineticSolver.compute(kineticSystem);
        if(kineticSolver.info()!=Eigen::Success){
            const Eigen::MatrixXd dense(kineticSystem),a=dense.bottomLeftCorner(constraintCount,kineticDofs);
            Eigen::JacobiSVD<Eigen::MatrixXd> rank(a);int resolved=0;
            for(int k=0;k<rank.singularValues().size();k++)if(rank.singularValues()[k]>1e-11*rank.singularValues()[0])resolved++;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> metric(dense.topLeftCorner(kineticDofs,kineticDofs),Eigen::EigenvaluesOnly);
            std::ostringstream message;message<<"RZ C1 constrained kinetic space lost rank: dofs="<<kineticDofs<<" constraints="<<constraintCount
                <<" constraint_rank="<<resolved<<" singular_min="<<rank.singularValues()[constraintCount-1]<<" kinetic_min="<<metric.eigenvalues()[0];
            throw std::invalid_argument(message.str());
        }
    }
    Layout(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentBoundary* boundary){
        space.validate(mesh,domain);rz=domain.axisymmetric();
        if(boundary)boundary->validate(mesh,domain,space);
        const bool general=boundary&&boundary->usesGeneralTraceSpace();boundaryKinds.resize(domain.nativeFaceCount(),HydroCompatibleMomentBoundary::Kind::Wall);
        if(boundary)for(const auto f:domain.faces(mesh))boundaryKinds[f.nativeRow]=boundary->faces()[f.nativeRow].kind;
        if(!domain.uniformLeaves()||general){
            native=std::make_shared<const NativeLayout>(mesh,domain,space,boundary);lo=native->lower;hi=native->upper;
            identity=domain.meshIdentity();velocityCount=int(space.mass().rows());
            for(const auto f:domain.faces(mesh))if(f.outwardSign&&HydroCompatibleMomentBoundary::impermeable(boundaryKinds[f.nativeRow])&&space.faces()[f.nativeRow].first>=0)
                for(int k=0;k<2;k++)boundaryRows.push_back(space.faces()[f.nativeRow].first+k);
            return;
        }
        UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);nx=grid.nx();ny=grid.ny();nvx=2*nx+1;nvy=2*ny+1;
        h=grid.spacing();lo=grid.minCoord();hi=grid.maxCoord();axis=rz&&lo.x==0;require(!rz||lo.x>=0,"negative RZ spline radius");
        identity=domain.meshIdentity();velocityCount=int(space.mass().rows());if(!rz){x=quartic(nx);y=quartic(ny);}
        cells.resize(nx*ny);for(int k=0;k<nx*ny;k++)cells[k]=mesh.cells(grid.cells()[k]).index();
        std::vector<int> vertical((nx+1)*ny,-1),horizontal(nx*(ny+1),-1);faces.resize(domain.nativeFaceCount());
        for(const auto f:domain.faces(mesh)){
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const int i=int(std::lround((a.x-lo.x)/h.x)),j=int(std::lround((a.y-lo.y)/h.y));
            faces[f.nativeRow]={{a,edge.midpointCoord(),b},i,j,f.horizontal,f.outwardSign!=0};
            (f.horizontal?horizontal[j*nx+i]:vertical[j*(nx+1)+i])=f.nativeRow;
            if(f.outwardSign&&space.faces()[f.nativeRow].first>=0)for(int k=0;k<2;k++)boundaryRows.push_back(space.faces()[f.nativeRow].first+k);
        }
        // Streamfunction on vertices: mean of the vertical and horizontal
        // paths from the zero wall gauge. On null(D), both paths coincide.
        // Their linear extension gives a defined transpose off null(D).
        std::vector<Row> nodes((nx+1)*(ny+1));
        for(int j=1;j<ny;j++)for(int i=1;i<nx;i++){
            auto& node=nodes[j*(nx+1)+i];
            for(int k=0;k<j;k++)node[space.faces()[vertical[k*(nx+1)+i]].first]+=.5*h.y*(rz?lo.x+i*h.x:1.);
            for(int k=0;k<i;k++)node[space.faces()[horizontal[j*nx+k]].first]-=.5*h.x*(rz?lo.x+(k+.5)*h.x:1.);
        }
        std::vector<Row> rows(nvx*nvy);const auto at=[&](int i,int j)->Row&{return rows[j*nvx+i];};
        for(int j=0;j<=ny;j++)for(int i=0;i<=nx;i++)at(i,j)=nodes[j*(nx+1)+i];
        for(int j=0;j<ny;j++)for(int i=0;i<=nx;i++){
            auto& target=at(i,ny+1+j);add(target,nodes[j*(nx+1)+i],.5);add(target,nodes[(j+1)*(nx+1)+i],.5);
            if(i>0&&i<nx){const auto& f=space.faces()[vertical[j*(nx+1)+i]];target[f.first+1]-=f.tangentScale*(rz?lo.x+i*h.x:1.);}
        }
        for(int j=0;j<=ny;j++)for(int i=0;i<nx;i++){
            auto& target=at(nx+1+i,j);const auto& f=space.faces()[horizontal[j*nx+i]];
            const double rc=lo.x+(i+.5)*h.x,shift=rz?(f.tangentCenter-rc)/h.x:0;
            add(target,nodes[j*(nx+1)+i],.5+shift);add(target,nodes[j*(nx+1)+i+1],.5-shift);
            if(j>0&&j<ny)target[f.first+1]+=f.tangentScale*(rz?rc:1.);
        }
        for(int j=0;j<ny;j++)for(int i=0;i<nx;i++){
            auto& target=at(nx+1+i,ny+1+j);const double scale=1/(1/h.x+1/h.y);
            target[space.firstCellVelocity(cells[j*nx+i])+2]=scale*(rz?lo.x+(i+.5)*h.x:1.);
            add(target,at(i,ny+1+j),.5*scale/h.x);add(target,at(i+1,ny+1+j),.5*scale/h.x);
            add(target,at(nx+1+i,j),.5*scale/h.y);add(target,at(nx+1+i,j+1),.5*scale/h.y);
        }
        if(rz){
            const int oldWidth=nvx;nvx-=int(axis);std::vector<Row> chi(nvx*nvy);
            for(int j=0;j<nvy;j++)for(int i=int(axis);i<oldWidth;i++){
                const double r=i<=nx?lo.x+i*h.x:lo.x+(i-nx-.5)*h.x;
                add(chi[j*nvx+i-int(axis)],rows[j*oldWidth+i],1/(r*r));
            }
            buildKinetic(space,chi);return;
        }
        std::vector<Eigen::Triplet<double>> entries;for(int r=0;r<int(rows.size());r++)for(const auto& v:rows[r])if(v.second!=0)entries.emplace_back(r,v.first,v.second);
        data.resize(nvx*nvy,velocityCount);data.setFromTriplets(entries.begin(),entries.end());
    }
};
HydroCompatibleMomentSpline HydroCompatibleMomentSpline::bind(std::shared_ptr<const Layout> layout,const HydroCompatibleMomentSpace& space,const Eigen::VectorXd& u,double time){
    require(u.size()==layout->velocityCount&&u.allFinite()&&std::isfinite(time),"invalid moment spline velocity/time");
    for(int row:layout->boundaryRows)require(u[row]==0,"moment C1 velocity requires the owned zero normal wall trace");
    const Eigen::VectorXd divergence=space.divergence()*u,scale=space.divergence().cwiseAbs()*u.cwiseAbs();
    require(divergence.stableNorm()<=1e-10*scale.stableNorm(),"moment C1 streamfunction requires a solenoidal state");
    HydroCompatibleMomentSpline result;result.layout_=std::move(layout);result.time_=time;const auto& l=*result.layout_;
    if(l.native)result.coefficients_=l.native->coefficients(u);
    else if(l.rz){
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(l.kineticSystem.rows());rhs.tail(l.data.rows())=l.data*u;
        const Eigen::VectorXd solution=l.solveKinetic(rhs),polynomial=l.polynomial*solution.head(l.kineticDofs);
        result.coefficients_=Eigen::Map<const Eigen::MatrixXd>(polynomial.data(),6*l.nx,6*l.ny);
        result.interior_=polynomial.tail(4*l.nx*l.ny);
    }else {const Eigen::VectorXd values=l.data*u;const Eigen::Map<const Eigen::MatrixXd> moments(values.data(),l.nvx,l.nvy);
        result.coefficients_=l.x*moments*l.y.transpose();}
    require(result.coefficients_.allFinite(),"nonfinite moment streamfunction");
    return result;
}
HydroCompatibleMomentSpline HydroCompatibleMomentSpline::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const Eigen::VectorXd& u,double time,const HydroCompatibleMomentBoundary* boundary){
    return bind(std::make_shared<const Layout>(mesh,domain,space,boundary),space,u,time);
}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::initializeNative(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const Eigen::VectorXd& u,double time,const HydroCompatibleMomentBoundary* boundary){
    space.validate(mesh,domain);if(boundary)boundary->validate(mesh,domain,space);
    require(u.size()==space.mass().rows()&&u.allFinite(),"invalid complete initial native moments");
    // With exactly zero data, the minimum physical kinetic energy is already
    // attained by the zero field. No rank reveal or Riesz solve is needed.
    const bool zero=u.isZero(0.);
    auto layout=std::make_shared<const NativeLayout>(mesh,domain,space,boundary,!zero);
    if(zero){const Eigen::VectorXd q=Eigen::VectorXd::Zero(layout->dofs);return bindNative(std::move(layout),q,time);}
    const Eigen::VectorXd q=layout->fieldCoordinates(layout->coefficients(u));return bindNative(std::move(layout),q,time);
}
HydroCompatibleMomentSpline HydroCompatibleMomentSpline::withVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const Eigen::VectorXd& u,double time) const{
    space.validate(mesh,domain);require(layout_&&layout_->identity==domain.meshIdentity(),"moment spline geometry identity changed");return bind(layout_,space,u,time);
}
HydroCompatibleMomentSpline::RefinementSource HydroCompatibleMomentSpline::prepareNativeRefinement(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentBoundary* boundary) const{
    require(layout_&&layout_->native,"exact native refinement requires the complete Q8 source representation");space.validate(mesh,domain);
    if(boundary)boundary->validate(mesh,domain,space);const bool general=boundary&&boundary->usesGeneralTraceSpace();
    require(general==layout_->native->generalTrace&&domain.axisymmetric()==layout_->rz,"native refinement changed boundary/measure kind");
    for(const auto f:domain.faces(mesh)){
        const auto kind=boundary?boundary->faces()[f.nativeRow].kind:HydroCompatibleMomentBoundary::Kind::Wall;
        if(f.outwardSign&&f.area>0){
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            double covered=0;for(const auto& old:layout_->native->faces)if(old.boundary&&old.horizontal==f.horizontal){
                if((f.horizontal?old.chart.begin.y:old.chart.begin.x)!=(f.horizontal?a.y:a.x))continue;
                const double lo=std::max(f.horizontal?a.x:a.y,f.horizontal?old.chart.begin.x:old.chart.begin.y),hi=std::min(f.horizontal?b.x:b.y,f.horizontal?old.chart.end.x:old.chart.end.y);
                if(hi>lo){require(old.wall==(kind==HydroCompatibleMomentBoundary::Kind::Wall)&&old.impermeable==HydroCompatibleMomentBoundary::impermeable(kind),"native refinement changed a physical wall/port partition");covered+=hi-lo;}}
            const double length=f.horizontal?b.x-a.x:b.y-a.y;require(std::abs(covered-length)<64*std::numeric_limits<double>::epsilon()*length,"native refinement boundary is not covered by source");
        }
    }
    RefinementSource result;result.layout_=std::make_shared<const NativeLayout>(mesh,domain,space,boundary);result.time_=time_;
    result.coefficients_=result.layout_->restrictCoefficients(*layout_->native,coefficients_);
    result.moments_=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=result.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);});
    return result;
}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentSpline::RefinementSource::firstJet(Vector2D p) const{
    require(bool(layout_),"moved native refinement source");return layout_->firstJet(coefficients_,p);
}
InterfaceNormalFluxView HydroCompatibleMomentSpline::RefinementSource::normalFlux() const{
    require(bool(layout_),"moved native refinement source");return layout_->normalFlux(coefficients_);
}
HydroCompatibleMomentSpline::NativeMaterialForms HydroCompatibleMomentSpline::nativeMaterialForms(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentMaterial& material,double tolerance) const{
    material.validate(mesh,domain,space);require(layout_&&layout_->native&&layout_->identity==domain.meshIdentity(),"native material forms require the same complete Q8 layout");
    require(material.geometryTime()==time_,"native material forms changed the physical time");return layout_->native->materialForms(material,coefficients_,tolerance);
}
HydroCompatibleMomentSpline::NativeMaterialForms HydroCompatibleMomentSpline::RefinementSource::nativeMaterialForms(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentMaterial& material,double tolerance) const{
    material.validate(mesh,domain,space);require(layout_&&layout_->identity==domain.meshIdentity(),"refined material forms changed the native layout");
    require(material.geometryTime()==time_,"refined material forms changed the physical time");return layout_->materialForms(material,coefficients_,tolerance);
}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::bindNative(std::shared_ptr<const NativeLayout> layout,const Eigen::VectorXd& coordinates,double time){
    require(layout&&coordinates.size()==layout->dofs&&coordinates.allFinite()&&std::isfinite(time),"invalid complete native velocity coordinates/time");
    require((layout->kinematicConstraints*coordinates).norm()<=2e-9*coordinates.norm(),"complete native velocity violates wall or axis conditions");
    NativeVelocity result;result.layout_=std::move(layout);result.coordinates_=coordinates;result.time_=time;
    const Eigen::VectorXd flat=result.layout_->polynomial*coordinates;result.coefficients_=Eigen::Map<const Eigen::MatrixXd>(flat.data(),81,result.layout_->cells.size());
    result.moments_=result.layout_->momentMap*coordinates;return result;
}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::nativeVelocity() const{
    require(layout_&&layout_->native,"complete velocity requires the native Q8 representation");return bindNative(layout_->native,layout_->native->fieldCoordinates(coefficients_),time_);
}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::RefinementSource::nativeVelocity() const{
    require(bool(layout_),"moved complete refinement source");return bindNative(layout_,layout_->fieldCoordinates(coefficients_),time_);
}
const HydroCompatibleMomentSpace::Matrix& HydroCompatibleMomentSpline::NativeVelocity::momentMap() const{
    require(bool(layout_),"moved complete native velocity");return layout_->momentMap;
}
const HydroCompatibleMomentSpace::Matrix& HydroCompatibleMomentSpline::NativeVelocity::kinematicConstraints() const{
    require(bool(layout_),"moved complete native velocity");return layout_->kinematicConstraints;
}
const std::vector<std::array<int,9>>& HydroCompatibleMomentSpline::NativeVelocity::cellInteriorCoordinates() const{
    require(bool(layout_),"moved complete native velocity");return layout_->cellInteriorCoordinates;
}
void HydroCompatibleMomentSpline::NativeVelocity::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);require(layout_&&layout_->identity==domain.meshIdentity()&&layout_->momentMap.rows()==space.mass().rows(),"complete native velocity changed its source layout");
}
HydroCompatibleMomentSpace::Matrix HydroCompatibleMomentSpline::NativeVelocity::normalBoundaryConstraints(const HydroCompatibleMomentBoundary& boundary) const{
    require(bool(layout_),"moved complete native velocity");return layout_->normalBoundaryConstraints(boundary);
}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::NativeVelocity::withCoordinates(const Eigen::VectorXd& q,double time) const{return bindNative(layout_,q,time);}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentSpline::NativeVelocity::firstJet(Vector2D p) const{
    require(bool(layout_),"moved complete native velocity");return layout_->firstJet(coefficients_,p);
}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentSpline::NativeVelocity::transportJet(Vector2D p) const{
    require(bool(layout_),"moved complete native velocity");return layout_->transportJet(coefficients_,p);
}
Vector2D HydroCompatibleMomentSpline::NativeVelocity::transportLower() const{
    require(bool(layout_),"moved complete native velocity");const double infinity=std::numeric_limits<double>::infinity();
    return {layout_->outletSpans[0].empty()?layout_->lower.x:-infinity,layout_->outletSpans[2].empty()?layout_->lower.y:-infinity};
}
Vector2D HydroCompatibleMomentSpline::NativeVelocity::transportUpper() const{
    require(bool(layout_),"moved complete native velocity");const double infinity=std::numeric_limits<double>::infinity();
    return {layout_->outletSpans[1].empty()?layout_->upper.x:infinity,layout_->outletSpans[3].empty()?layout_->upper.y:infinity};
}
InterfaceNormalFluxView HydroCompatibleMomentSpline::NativeVelocity::normalFlux() const{
    require(bool(layout_),"moved complete native velocity");return layout_->normalFlux(coefficients_,true);
}
Vector2D HydroCompatibleMomentSpline::NativeVelocity::lower() const{require(bool(layout_),"moved complete native velocity");return layout_->lower;}
Vector2D HydroCompatibleMomentSpline::NativeVelocity::upper() const{require(bool(layout_),"moved complete native velocity");return layout_->upper;}
bool HydroCompatibleMomentSpline::NativeVelocity::sameLayout(const NativeVelocity& other) const{return layout_&&layout_==other.layout_;}
bool HydroCompatibleMomentSpline::NativeVelocity::openBoundary() const{require(bool(layout_),"moved complete native velocity");return layout_->permeable;}
bool HydroCompatibleMomentSpline::NativeVelocity::usesGeneralTraceSpace() const{require(bool(layout_),"moved complete native velocity");return layout_->generalTrace;}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::NativeVelocity::refine(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentBoundary* boundary) const{
    require(bool(layout_),"moved complete native velocity");space.validate(mesh,domain);if(boundary)boundary->validate(mesh,domain,space);
    const bool general=boundary&&boundary->usesGeneralTraceSpace();require(general==layout_->generalTrace&&domain.axisymmetric()==layout_->rz,"native refinement changed boundary/measure kind");
    NativeVelocity target;target.layout_=std::make_shared<const NativeLayout>(mesh,domain,space,boundary,false);return target.restrictFrom(*this);
}
HydroCompatibleMomentSpline::NativeVelocity HydroCompatibleMomentSpline::NativeVelocity::restrictFrom(const NativeVelocity& source) const{
    require(layout_&&source.layout_&&layout_->generalTrace==source.layout_->generalTrace&&layout_->permeable==source.layout_->permeable&&layout_->rz==source.layout_->rz
        &&layout_->lower.x==source.layout_->lower.x&&layout_->lower.y==source.layout_->lower.y&&layout_->upper.x==source.layout_->upper.x&&layout_->upper.y==source.layout_->upper.y,
        "native restriction changed physical domain");
    for(const auto& face:layout_->faces)if(face.boundary){double covered=0;
        const double a=face.horizontal?face.chart.begin.x:face.chart.begin.y,b=face.horizontal?face.chart.end.x:face.chart.end.y;
        for(const auto& old:source.layout_->faces)if(old.boundary&&old.horizontal==face.horizontal){
            if((face.horizontal?face.chart.begin.y:face.chart.begin.x)!=(old.horizontal?old.chart.begin.y:old.chart.begin.x))continue;
            const double lo=std::max(a,old.horizontal?old.chart.begin.x:old.chart.begin.y),hi=std::min(b,old.horizontal?old.chart.end.x:old.chart.end.y);
            if(hi>lo){require(face.wall==old.wall&&face.impermeable==old.impermeable,"native restriction changed wall/port partition");covered+=hi-lo;}}
        require(std::abs(covered-(b-a))<=64*std::numeric_limits<double>::epsilon()*(b-a),"native restriction boundary not covered");
    }
    const auto coefficients=layout_->restrictCoefficients(*source.layout_,source.coefficients_);
    return bindNative(layout_,layout_->fieldCoordinates(coefficients),source.time_);
}
HydroCompatibleMomentSpline::NativeCellBalance HydroCompatibleMomentSpline::NativeVelocity::cellBalance(const NativeMaterialForms& forms) const{
    require(layout_&&forms.solenoidalCoordinates==coordinates_.size()&&forms.cells.size()==layout_->cells.size(),"material balance belongs to another complete layout");
    NativeCellBalance result;result.momentum=Eigen::MatrixXd::Zero(2,layout_->slotByCell.size());result.kinetic=Eigen::VectorXd::Zero(layout_->slotByCell.size());
    for(const auto& f:forms.cells){const int slot=layout_->slotByCell.at(f.index);const auto u=coefficients_.col(slot);
        result.momentum.col(f.index)=f.momentum.leftCols(81)*u;result.kinetic[f.index]=.5*u.dot(f.mass.topLeftCorner(81,81)*u);}
    return result;
}
HydroCompatibleMomentSpace::Matrix HydroCompatibleMomentSpline::NativeVelocity::cellMomentumMap(const NativeMaterialForms& forms) const{
    require(layout_&&forms.solenoidalCoordinates==coordinates_.size()&&forms.cells.size()==layout_->cells.size(),"momentum map belongs to another complete layout");
    std::vector<Eigen::Triplet<double>> entries;
    for(const auto& f:forms.cells){const int slot=layout_->slotByCell.at(f.index);
        for(int d=0;d<2;d++)for(int k=0;k<81;k++)if(f.momentum(d,k)!=0)entries.emplace_back(2*f.index+d,81*slot+k,f.momentum(d,k));}
    HydroCompatibleMomentSpace::Matrix local(2*layout_->slotByCell.size(),layout_->polynomial.rows());local.setFromTriplets(entries.begin(),entries.end());
    return local*layout_->polynomial;
}
HydroCompatibleMomentSpline::NativeProjectionForms HydroCompatibleMomentSpline::NativeVelocity::projectionForms(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const HydroCompatibleMomentMaterial& material,const NativeVelocity& source,const HydroCompatibleMomentComplement* complement,double tolerance) const{
    validate(mesh,domain,space);material.validate(mesh,domain,space);
    if(complement)complement->validate(mesh,domain,space);
    require(source.layout_&&source.time_==time_&&material.geometryTime()==time_,"projection changed physical time");
    require(layout_->rz==source.layout_->rz&&layout_->generalTrace==source.layout_->generalTrace&&layout_->permeable==source.layout_->permeable
        &&lower().x==source.lower().x&&lower().y==source.lower().y&&upper().x==source.upper().x&&upper().y==source.upper().y,"projection changed physical domain");
    NativeProjectionForms result;
    if(sameLayout(source)){
        result.material=layout_->materialForms(material,coefficients_,tolerance,complement);
        result.load=result.material.mass.topLeftCorner(coordinates_.size(),coordinates_.size())*source.coordinates_;result.sourceNormSquared=source.coordinates_.dot(result.load);
    }else{
        NativeLayout::ProjectionRequest request{*source.layout_,source.coefficients_};
        result.material=layout_->materialForms(material,coefficients_,tolerance,complement,&request);
        result.load=request.load.head(coordinates_.size());result.sourceNormSquared=request.sourceNormSquared;
    }
    return result;
}
HydroCompatibleMomentSpline::NativeFaceProducts HydroCompatibleMomentSpline::NativeVelocity::faceProducts(int row) const{
    require(bool(layout_),"moved complete native velocity");return layout_->faceProducts(coefficients_,row);
}
HydroCompatibleMomentSpline::NativeVolumeTransport HydroCompatibleMomentSpline::NativeVelocity::nativeVolumeTransport(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const HydroCompatibleMomentMassCurrent& current,const HydroCompatibleMomentComplement* complement) const{
    validate(mesh,domain,space);current.validate(mesh,domain,space);if(complement)complement->validate(mesh,domain,space);
    require(current.time()==time_,"volume momentum uses a different midpoint time");return layout_->volumeTransport(current,complement);
}
HydroCompatibleMomentSpline::NativeMomentumTransport HydroCompatibleMomentSpline::NativeVelocity::nativeMomentumTransport(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,
    const HydroCompatibleMomentMassCurrent& current,const HydroCompatibleMomentBoundary& boundary,const HydroCompatibleMomentComplement* complement) const{
    validate(mesh,domain,space);current.validate(mesh,domain,space);boundary.validate(mesh,domain,space);if(complement)complement->validate(mesh,domain,space);
    require(current.time()==time_,"complete momentum uses a different midpoint time");return layout_->momentumTransport(current,boundary,complement,coordinates_);
}
Eigen::VectorXd HydroCompatibleMomentSpline::NativeVelocity::transposePointForces(const std::vector<PointForce>& forces) const{
    require(bool(layout_),"moved complete native velocity");return layout_->nativeTranspose(forces);
}
HydroCompatibleMomentSpace::Matrix HydroCompatibleMomentSpline::NativeVelocity::pointWorkGram(const std::vector<PointForce>& forces) const{
    require(bool(layout_),"moved complete native velocity");return layout_->pointWorkGram(forces);
}
HydroCompatibleMomentSpline::NativeMaterialForms HydroCompatibleMomentSpline::NativeVelocity::nativeMaterialForms(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentMaterial& material,double tolerance) const{
    material.validate(mesh,domain,space);require(layout_&&layout_->identity==domain.meshIdentity(),"complete native material forms changed the layout");
    require(material.geometryTime()==time_,"complete native material forms changed the physical time");return layout_->materialForms(material,coefficients_,tolerance);
}
HydroCompatibleMomentSpline::NativeMaterialForms HydroCompatibleMomentSpline::NativeVelocity::mixedMaterialForms(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentMaterial& material,
    const HydroCompatibleMomentComplement& complement,double tolerance) const{
    material.validate(mesh,domain,space);complement.validate(mesh,domain,space);
    require(layout_&&layout_->identity==domain.meshIdentity(),"mixed native material forms changed the layout");
    require(material.geometryTime()==time_,"mixed native material forms changed the physical time");
    return layout_->materialForms(material,coefficients_,tolerance,&complement);
}
Vector2D HydroCompatibleMomentSpline::lower() const{return layout_->lo;}
Vector2D HydroCompatibleMomentSpline::upper() const{return layout_->hi;}
bool HydroCompatibleMomentSpline::sameLayout(const HydroCompatibleMomentSpline& other) const{return layout_&&other.layout_&&layout_->identity==other.layout_->identity&&layout_->rz==other.layout_->rz&&layout_->boundaryKinds==other.layout_->boundaryKinds;}
bool HydroCompatibleMomentSpline::openBoundary() const{return layout_&&layout_->native&&layout_->native->permeable;}
bool HydroCompatibleMomentSpline::usesGeneralTraceSpace() const{return layout_&&layout_->native&&layout_->native->generalTrace;}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentSpline::firstJet(Vector2D p) const{
    const auto& l=*layout_;if(l.native)return l.native->firstJet(coefficients_,p);const auto q=l.locate(p);
    // All jets share these powers; flow-map quadrature queries this field
    // repeatedly. Do not evaluate pow for every coefficient and derivative.
    double px[3][6]{},py[3][6]{};px[0][0]=py[0][0]=1;
    for(int k=1;k<l.width;k++){px[0][k]=px[0][k-1]*q.x;py[0][k]=py[0][k-1]*q.y;}
    for(int k=1;k<l.width;k++){px[1][k]=k*px[0][k-1]/l.h.x;py[1][k]=k*py[0][k-1]/l.h.y;}
    for(int k=2;k<l.width;k++){px[2][k]=k*(k-1)*px[0][k-2]/(l.h.x*l.h.x);py[2][k]=k*(k-1)*py[0][k-2]/(l.h.y*l.h.y);}
    const auto derivative=[&](int dx,int dy){double out=0;for(int a=dx;a<l.width;a++)for(int b=dy;b<l.width;b++)out+=coefficients_(l.width*q.i+a,l.width*q.j+b)*px[dx][a]*py[dy][b];return out;};
    if(l.rz){
        std::array<double,6> chi{derivative(0,0),derivative(1,0),derivative(0,1),derivative(2,0),derivative(1,1),derivative(0,2)};
        const int cell=q.j*l.nx+q.i;
        const auto bubbles=interiorJets(q.x,q.y,l.lo.x+(q.i+.5)*l.h.x,l.h.x,l.axis&&q.i==0);
        const double scaling[6]={1,1/l.h.x,1/l.h.y,1/(l.h.x*l.h.x),1/(l.h.x*l.h.y),1/(l.h.y*l.h.y)};
        for(int b=0;b<4;b++)for(int k=0;k<6;k++)chi[k]+=interior_[4*cell+b]*bubbles[b][k]*scaling[k];
        // Psi=r^2 chi gives a polynomial physical velocity, including r=0.
        // No division by r, radial clipping or axis epsilon is needed.
        FirstJet out;out.nativeCell=l.cells[cell];out.value={p.x*chi[2],-2*chi[0]-p.x*chi[1]};
        out.dx={chi[2]+p.x*chi[4],-3*chi[1]-p.x*chi[3]};out.dy={p.x*chi[5],-2*chi[2]-p.x*chi[4]};
        if(p.x==0){out.value.x=0;out.dy.x=0;out.dx.y=0;}
        else if(p.x==l.lo.x||p.x==l.hi.x){out.value={0,0};out.dy={0,0};out.dx.x=0;}
        if(p.y==l.lo.y||p.y==l.hi.y){out.value={0,0};out.dx={0,0};out.dy.y=0;}
        return out;
    }
    FirstJet out;out.nativeCell=l.cells[q.j*l.nx+q.i];out.value={derivative(0,1),-derivative(1,0)};
    const double mixed=derivative(1,1);out.dx={mixed,-derivative(2,0)};out.dy={derivative(0,2),-mixed};
    // Essential no-slip traces are known exactly. Their tangential derivatives
    // are zero; no clipping is applied at an interior point.
    if(p.x==l.lo.x||p.x==l.hi.x){out.value={0,0};out.dy={0,0};out.dx.x=0;}
    if(p.y==l.lo.y||p.y==l.hi.y){out.value={0,0};out.dx={0,0};out.dy.y=0;}
    return out;
}
InterfaceNormalFluxView HydroCompatibleMomentSpline::normalFlux() const{
    if(layout_->native)return layout_->native->normalFlux(coefficients_);
    const auto layout=layout_;auto polynomials=std::make_shared<std::vector<std::array<Quantum::Geometry::UnitPolynomialTerm,2>>>(layout->faces.size());
    for(int r=0;r<int(layout->faces.size());r++){
        const auto& f=layout->faces[r];if(f.boundary)continue;
        // Differentiate the SAME streamfunction along the face.
        // Its half-parameter derivative already includes the physical ds.
        const int i=f.horizontal?f.i:f.i-1,j=f.horizontal?f.j-1:f.j;
        const int width=layout->width;double primitive[8]{};for(int k=0;k<width;k++)for(int normal=0;normal<width;normal++)
            primitive[k]+=(f.horizontal?-1:1)*coefficients_(width*i+(f.horizontal?k:normal),width*j+(f.horizontal?normal:k));
        int degree=width-1;
        if(layout->rz){
            if(f.horizontal){double chi[6];for(int k=0;k<width;k++){chi[k]=primitive[k];primitive[k]=0;}
                const double r0=layout->lo.x+i*layout->h.x,h=layout->h.x;
                for(int k=0;k<width;k++){primitive[k]+=r0*r0*chi[k];primitive[k+1]+=2*r0*h*chi[k];primitive[k+2]+=h*h*chi[k];}degree=width+1;
            }else for(double& c:primitive)c*=f.chart.begin.x*f.chart.begin.x;
            for(double& c:primitive)c*=CartesianMeasure::twoPi;
        }
        for(int half=0;half<2;half++){
            auto& term=(*polynomials)[r][half];term.factor=.5;
            for(int k=1;k<=degree;k++)for(int a=0;a<k;a++){
                double choose=1;for(int s=1;s<=a;s++)choose*=double(k-s)/s;
                term.coefficients[a]+=k*primitive[k]*choose*std::pow(.5*half,k-1-a)*std::pow(.5,a);
            }
        }
    }
    InterfaceNormalFluxView::Queries q;q.retainedMomentDegree=2;q.chart=[layout](int r){return layout->faces.at(r).chart;};
    q.polynomial=[polynomials](int r,int half){require(half==0||half==1,"invalid moment face half");return std::vector<Quantum::Geometry::UnitPolynomialTerm>{polynomials->at(r)[half]};};
    q.integral=[polynomial=q.polynomial](int r,int h,double a,double b,double c,double d){return Quantum::Geometry::integrateUnitPolynomial(polynomial(r,h),a,b,c,d);};
    q.density=[polynomial=q.polynomial](int r,int h,bool end,long double offset){return Quantum::Geometry::evaluateUnitPolynomial(polynomial(r,h),double(end?1-offset:offset));};
    q.phaseIntegral=[](int,int,const LevelSetMappedSegment&)->long double{throw std::invalid_argument("moment spline exports parametric phase flux, not legacy mapped-segment geometry");};
    return InterfaceNormalFluxView(std::move(q));
}
Eigen::VectorXd HydroCompatibleMomentSpline::transposePointForces(const std::vector<PointForce>& forces) const{
    if(layout_->native)return layout_->native->transpose(coefficients_,forces);
    const auto& l=*layout_;Eigen::MatrixXd coefficients=Eigen::MatrixXd::Zero(l.width*l.nx,l.width*l.ny);
    if(l.rz){
        Eigen::VectorXd interior=Eigen::VectorXd::Zero(4*l.nx*l.ny);
        for(auto f:forces){
            require(std::isfinite(f.measure.x)&&std::isfinite(f.measure.y)&&std::isfinite(f.dx.x)&&std::isfinite(f.dx.y)&&std::isfinite(f.dy.x)&&std::isfinite(f.dy.y),"nonfinite RZ moment work");
            const auto p=l.locate(f.point);
            if(f.point.x==0){f.measure.x=0;f.dy.x=0;f.dx.y=0;}
            else if(f.point.x==l.lo.x||f.point.x==l.hi.x){f.measure={0,0};f.dy={0,0};f.dx.x=0;}
            if(f.point.y==l.lo.y||f.point.y==l.hi.y){f.measure={0,0};f.dx={0,0};f.dy.y=0;}
            const double r=f.point.x;const std::array<double,6> work{-2*f.measure.y,(-r*f.measure.y-3*f.dx.y)/l.h.x,
                (r*f.measure.x+f.dx.x-2*f.dy.y)/l.h.y,-r*f.dx.y/(l.h.x*l.h.x),r*(f.dx.x-f.dy.y)/(l.h.x*l.h.y),r*f.dy.x/(l.h.y*l.h.y)};
            for(int a=0;a<6;a++)for(int b=0;b<6;b++)coefficients(6*p.i+a,6*p.j+b)+=
                work[0]*basis(a,p.x,0)*basis(b,p.y,0)+work[1]*basis(a,p.x,1)*basis(b,p.y,0)+work[2]*basis(a,p.x,0)*basis(b,p.y,1)
                +work[3]*basis(a,p.x,2)*basis(b,p.y,0)+work[4]*basis(a,p.x,1)*basis(b,p.y,1)+work[5]*basis(a,p.x,0)*basis(b,p.y,2);
            const auto bubbles=interiorJets(p.x,p.y,l.lo.x+(p.i+.5)*l.h.x,l.h.x,l.axis&&p.i==0);
            for(int b=0;b<4;b++)for(int k=0;k<6;k++)interior[4*(p.j*l.nx+p.i)+b]+=work[k]*bubbles[b][k];
        }
        Eigen::VectorXd flat(40*l.nx*l.ny);flat.head(coefficients.size())=Eigen::Map<const Eigen::VectorXd>(coefficients.data(),coefficients.size());flat.tail(interior.size())=interior;
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(l.kineticSystem.rows());rhs.head(l.kineticDofs)=l.polynomial.transpose()*flat;
        const Eigen::VectorXd solution=l.solveKinetic(rhs);return l.data.transpose()*solution.tail(l.data.rows());
    }
    for(const auto& f:forces){require(std::isfinite(f.measure.x)&&std::isfinite(f.measure.y),"nonfinite moment point force");const auto p=l.locate(f.point);
        // Essential wall velocities do no work, including rounding-free trace
        // evaluation. Interior work is the exact transpose of firstJet.
        require(std::isfinite(f.dx.x)&&std::isfinite(f.dx.y)&&std::isfinite(f.dy.x)&&std::isfinite(f.dy.y),"nonfinite moment derivative work");
        const bool xWall=f.point.x==l.lo.x||f.point.x==l.hi.x,yWall=f.point.y==l.lo.y||f.point.y==l.hi.y;
        for(int a=0;a<5;a++)for(int b=0;b<5;b++){
            double value=0;
            if(!xWall&&!yWall)value+=f.measure.x*basis(a,p.x,0)*basis(b,p.y,1)/l.h.y-f.measure.y*basis(a,p.x,1)*basis(b,p.y,0)/l.h.x
                +(f.dx.x-f.dy.y)*basis(a,p.x,1)*basis(b,p.y,1)/(l.h.x*l.h.y);
            if(!yWall)value-=f.dx.y*basis(a,p.x,2)*basis(b,p.y,0)/(l.h.x*l.h.x);
            if(!xWall)value+=f.dy.x*basis(a,p.x,0)*basis(b,p.y,2)/(l.h.y*l.h.y);
            coefficients(5*p.i+a,5*p.j+b)+=value;
        }
    }
    const Eigen::MatrixXd moments=l.x.transpose()*coefficients*l.y;const Eigen::Map<const Eigen::VectorXd> flat(moments.data(),moments.size());return l.data.transpose()*flat;
}
}
