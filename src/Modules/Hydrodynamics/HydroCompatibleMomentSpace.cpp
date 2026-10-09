#include "Modules/Hydrodynamics/HydroCompatibleMomentSpace.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* reason){if(!ok)throw std::invalid_argument(reason);}
constexpr double gx[4]={-.8611363115940525752,-.3399810435848562648,.3399810435848562648,.8611363115940525752};
constexpr double gw[4]={.3478548451374538574,.6521451548625461426,.6521451548625461426,.3478548451374538574};
double weight(Vector2D p,bool rz){return rz?CartesianMeasure::twoPi*p.x:1.;}
struct Basis {
    Eigen::Vector3d pressure;
    Eigen::Matrix<double,5,1> potential;
    Eigen::Matrix<double,2,6> velocity;
};
Basis basis(Vector2D p,Vector2D center,Vector2D h){
    const double x=(p.x-center.x)/h.x,y=(p.y-center.y)/h.y,l=std::sqrt(h.x*h.y);
    Basis b;b.pressure<<1,x,y;b.potential<<l*x,l*y,l*x*x/2,l*x*y,l*y*y/2;
    b.velocity=HydroCompatibleMomentSpace::affineVelocityBasis({x,y},h);
    return b;
}
template<class F>void faceQuadrature(const QuadtreeMesh& mesh,const NativeLeafFace& face,F&& f){
    const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
    const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
    const double length=face.horizontal?b.x-a.x:b.y-a.y;
    for(int k=0;k<4;k++){const double t=.5*(1+gx[k]);f(Vector2D{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)},.5*length*gw[k]);}
}
template<class F>void cellQuadrature(Vector2D center,Vector2D h,F&& f){
    for(int j=0;j<4;j++)for(int i=0;i<4;i++)f(Vector2D{center.x+.5*h.x*gx[i],center.y+.5*h.y*gx[j]},.25*h.x*h.y*gw[i]*gw[j]);
}
}
Eigen::Matrix<double,2,6> HydroCompatibleMomentSpace::affineVelocityBasis(Vector2D local,Vector2D h){
    const double l=std::sqrt(h.x*h.y),x=local.x,y=local.y;
    Eigen::Matrix<double,2,6> result;
    result<<l/h.x,0,l*x/h.x,l*y/h.x,0,-y,
        0,l/h.y,0,l*x/h.y,l*y/h.y,x;
    return result;
}
HydroCompatibleMomentSpace HydroCompatibleMomentSpace::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain){
    domain.validateMeshIdentity(mesh);
    require(domain.cellIndexByRow().size()==std::size_t(mesh.cellIndexNum()),"moment space requires the complete native domain");
    HydroCompatibleMomentSpace result;result.identity_=domain.meshIdentity();result.axisymmetric_=domain.axisymmetric();
    const bool rz=result.axisymmetric_;const int count=mesh.cellIndexNum();
    result.faces_.resize(domain.nativeFaceCount());result.cells_.resize(count);
    std::vector<std::vector<int>> incident(count);int velocityCount=0;
    for(const auto face:domain.faces(mesh)){
        auto& saved=result.faces_[face.nativeRow];saved.area=face.area;
        if(face.area==0){require(rz&&!face.horizontal&&face.center.x==0,"unexpected zero-area moment face");continue;}
        saved.first=velocityCount;velocityCount+=2;
        long double first=0,second=0;
        // Center differences avoid subtracting large squared coordinates.
        const double anchor=face.horizontal?face.center.x:face.center.y;
        faceQuadrature(mesh,face,[&](Vector2D p,double ds){const double t=(face.horizontal?p.x:p.y)-anchor;
            first+=ds*weight(p,rz)*t;second+=ds*weight(p,rz)*t*t;});
        const double shift=static_cast<double>(first/face.area);
        saved.tangentCenter=anchor+shift;saved.tangentScale=std::sqrt(static_cast<double>(second/face.area)-shift*shift);
        require(std::isfinite(saved.tangentScale)&&saved.tangentScale>0,"invalid face moment scale");
        if(face.negativeCell>=0)incident[face.negativeCell].push_back(face.nativeRow);
        if(face.positiveCell>=0)incident[face.positiveCell].push_back(face.nativeRow);
    }
    const int cellDofs=result.cellVelocityCount();
    result.cellVelocityBegin_=velocityCount;velocityCount+=cellDofs*count;
    std::vector<Eigen::Triplet<double>> massEntries,divergenceEntries;
    result.report_.minimumScaledEigenvalue=std::numeric_limits<double>::infinity();
    constexpr double tolerance=16384*std::numeric_limits<double>::epsilon();
    for(int c=0;c<count;c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[c]],flat=h.x*h.y;
        const int nf=int(incident[c].size()),first=2*nf,n=first+cellDofs;
        auto& block=result.cells_[c];block.rows.reserve(n);
        Eigen::MatrixXd samples=Eigen::MatrixXd::Zero(n,6),loads=Eigen::MatrixXd::Zero(n,6),div=Eigen::MatrixXd::Zero(3,n);
        Eigen::VectorXd baseline=Eigen::VectorXd::Constant(n,volume/(2*nf));
        Eigen::Matrix<double,3,5> potentialMoments=Eigen::Matrix<double,3,5>::Zero();
        for(int i=0;i<nf;i++){
            const auto face=domain.face(mesh,incident[c][i]);const auto& info=result.faces_[face.nativeRow];
            const double sign=face.negativeCell==c?1:-1;
            block.rows.push_back(info.first);block.rows.push_back(info.first+1);
            faceQuadrature(mesh,face,[&](Vector2D p,double ds){
                const auto b=basis(p,center,h);const double dw=ds*weight(p,rz);
                const double shape[2]={1,((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale};
                for(int m=0;m<2;m++){
                    samples.row(2*i+m)+=dw*shape[m]/face.area*b.velocity.row(face.horizontal?1:0);
                    loads.block(2*i+m,0,1,5)+=sign*dw*shape[m]*b.potential.transpose();
                    div.col(2*i+m)+=sign*dw*shape[m]*b.pressure;
                }
            });
            // The centered moment has identically zero integral. Preserve
            // this definition exactly, rather than a quadrature roundoff.
            div(0,2*i)=sign*face.area;div(0,2*i+1)=0;
        }
        for(int k=0;k<cellDofs;k++)block.rows.push_back(result.firstCellVelocity(c)+k);
        baseline[first]=baseline[first+1]=volume;baseline[first+2]=6*volume;
        if(rz)baseline.segment(first+3,3).setConstant(12*volume);
        div(1,first)=-volume/h.x;div(2,first+1)=-volume/h.y;
        cellQuadrature(center,h,[&](Vector2D p,double da){const auto b=basis(p,center,h);const double dv=da*weight(p,rz);
            samples.row(first)+=dv/volume*b.velocity.row(0);samples.row(first+1)+=dv/volume*b.velocity.row(1);
            samples.row(first+2)+=dv/volume*b.velocity.col(5).transpose()*b.velocity;
            potentialMoments+=da*b.pressure*b.potential.transpose();
        });
        if(rz){
            samples.middleRows(first+3,3)=rzExtraMoments(center,h,[&](Vector2D p)->Eigen::MatrixXd{return basis(p,center,h).velocity;});
            block.flatVelocityMoments=result.rzFlatMoments(mesh,domain,c,incident[c],div);
        }
        // div(w*u) is P1 in the flat (x,y)/(r,z) coordinates. Recover its
        // three coefficients using the flat polynomial Gram, without 1/r.
        Eigen::Vector3d inverseGram;inverseGram<<1/flat,12/flat,12/flat;
        loads.leftCols(5)-=div.transpose()*inverseGram.asDiagonal()*potentialMoments;
        loads(first+2,5)=volume; // the retained rotational cell moment
        const Eigen::VectorXd root=baseline.cwiseSqrt();
        const Eigen::MatrixXd u=root.asDiagonal()*samples,v=root.cwiseInverse().asDiagonal()*loads;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(u,Eigen::ComputeThinU|Eigen::ComputeThinV);
        require(svd.singularValues()[5]>tolerance*svd.singularValues()[0],"affine velocity moments lost rank");
        const Eigen::MatrixXd q=svd.matrixU(),z=v*svd.matrixV()*svd.singularValues().cwiseInverse().asDiagonal();
        const Eigen::MatrixXd a=q.transpose()*z;
        const double skew=(a-a.transpose()).norm()/a.norm();
        require(skew<tolerance,"incompatible pressure/velocity moments: no symmetric Hodge exists");
        const Eigen::MatrixXd symmetric=.5*(a+a.transpose());
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> aeig(symmetric);
        require(aeig.info()==Eigen::Success&&aeig.eigenvalues()[0]>tolerance*aeig.eigenvalues()[5],"nonpositive moment consistency Gram");
        const Eigen::MatrixXd scaled=z*symmetric.ldlt().solve(z.transpose())+Eigen::MatrixXd::Identity(n,n)-q*q.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> heig(scaled);
        require(heig.info()==Eigen::Success&&heig.eigenvalues()[0]>0,"nonpositive moment Hodge");
        const double defect=(scaled*u-v).norm()/v.norm();require(defect<tolerance,"moment Hodge consistency failed");
        block.mass=root.asDiagonal()*scaled*root.asDiagonal();
        const Eigen::MatrixXd gram=samples.transpose()*loads;
        block.affineProjection=(.5*(gram+gram.transpose())).ldlt().solve(loads.transpose());
        block.affineSamples=samples;
        result.report_.maximumConsistencyDefect=std::max(result.report_.maximumConsistencyDefect,defect);
        result.report_.maximumSymmetryDefect=std::max(result.report_.maximumSymmetryDefect,skew);
        result.report_.minimumScaledEigenvalue=std::min(result.report_.minimumScaledEigenvalue,heig.eigenvalues()[0]);
        result.report_.maximumScaledEigenvalue=std::max(result.report_.maximumScaledEigenvalue,heig.eigenvalues()[n-1]);
        result.report_.maximumLocalVelocityDofs=std::max(result.report_.maximumLocalVelocityDofs,n);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++)if(block.mass(i,j)!=0)massEntries.emplace_back(block.rows[i],block.rows[j],block.mass(i,j));
            for(int k=0;k<3;k++)if(div(k,i)!=0)divergenceEntries.emplace_back(3*c+k,block.rows[i],div(k,i));
        }
    }
    result.mass_.resize(velocityCount,velocityCount);result.mass_.setFromTriplets(massEntries.begin(),massEntries.end());
    result.divergence_.resize(3*count,velocityCount);result.divergence_.setFromTriplets(divergenceEntries.begin(),divergenceEntries.end());
    result.pressureLoad_=-result.divergence_.transpose();
    if(rz){
        // The RZ streamfunction moments can be small for a physical velocity
        // with substantial energy. A diagonal penalty in their raw coordinates
        // therefore does not define its kinetic scale. Use the Gram of the
        // SAME complete, moment-preserving minimum-L2 local reconstruction.
        // D, its negative transpose and the affine projection stay unchanged.
        const auto reconstruction=HydroCompatibleMomentTransport::build(mesh,domain,result);
        massEntries.clear();result.report_.minimumScaledEigenvalue=std::numeric_limits<double>::infinity();result.report_.maximumScaledEigenvalue=0;
        for(int c=0;c<count;c++){
            auto& block=result.cells_[c];const auto& mass=reconstruction.cells()[c].kineticMass;
            const Eigen::MatrixXd loads=block.mass*block.affineSamples;
            const double defect=(mass*block.affineSamples-loads).norm()/loads.norm();
            require(defect<2e-9,"physical RZ Hodge lost affine work consistency");
            const Eigen::VectorXd scale=mass.diagonal().cwiseSqrt();
            require((scale.array()>0).all(),"physical RZ Hodge has a zero-energy degree of freedom");
            const Eigen::MatrixXd scaled=scale.cwiseInverse().asDiagonal()*mass*scale.cwiseInverse().asDiagonal();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(scaled);
            require(eigen.info()==Eigen::Success&&eigen.eigenvalues()[0]>0,"physical RZ Hodge is not positive");
            result.report_.maximumConsistencyDefect=std::max(result.report_.maximumConsistencyDefect,defect);
            result.report_.minimumScaledEigenvalue=std::min(result.report_.minimumScaledEigenvalue,eigen.eigenvalues()[0]);
            result.report_.maximumScaledEigenvalue=std::max(result.report_.maximumScaledEigenvalue,eigen.eigenvalues()[mass.rows()-1]);
            block.mass=mass;
            for(int i=0;i<mass.rows();i++)for(int j=0;j<mass.cols();j++)if(mass(i,j)!=0)massEntries.emplace_back(block.rows[i],block.rows[j],mass(i,j));
        }
        result.mass_.setZero();result.mass_.setFromTriplets(massEntries.begin(),massEntries.end());
    }
    return result;
}
void HydroCompatibleMomentSpace::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const {
    domain.validateMeshIdentity(mesh);
    require(identity_&&identity_==domain.meshIdentity()&&axisymmetric_==domain.axisymmetric(),"moment metric/domain source differs");
}
HydroCompatibleMomentSpace::DivergenceComplement HydroCompatibleMomentSpace::divergenceComplement(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const std::vector<bool>& fixed) const{
    validate(mesh,domain);const int n=mesh.cellIndexNum(),nv=int(mass_.rows()),np=int(divergence_.rows());
    require(int(fixed.size())>=nv&&np==3*n,"invalid divergence-complement velocity/pressure layout");
    struct Link{int neighbor,face;};std::vector<std::vector<Link>> adjacency(n);
    DivergenceComplement result;result.rootCell=0;
    for(const auto f:domain.faces(mesh))if(f.area>0){
        const int row=faces_[f.nativeRow].first;require(row>=0,"positive native face has no moment");
        if(f.outwardSign){
            if(!fixed[row]&&result.outletFace<0){result.outletFace=f.nativeRow;result.rootCell=f.negativeCell>=0?f.negativeCell:f.positiveCell;}
        }else{
            require(!fixed[row],"divergence complement cannot cross a prescribed interior face");
            adjacency[f.negativeCell].push_back({f.positiveCell,f.nativeRow});adjacency[f.positiveCell].push_back({f.negativeCell,f.nativeRow});
        }
    }
    const bool closed=result.outletFace<0;const int omitted=3*result.rootCell,nc=np-(closed?1:0);
    std::vector<int> parent(n,-1),parentFace(n,-1),order{result.rootCell};parent[result.rootCell]=result.rootCell;
    for(std::size_t i=0;i<order.size();i++)for(const auto link:adjacency[order[i]])if(parent[link.neighbor]<0){
        parent[link.neighbor]=order[i];parentFace[link.neighbor]=link.face;order.push_back(link.neighbor);
    }
    require(int(order.size())==n,"divergence complement requires a connected native domain");
    std::vector<Eigen::Triplet<double>> range,values;
    for(int row=0,col=0;row<np;row++){
        if(closed&&row==omitted)continue;
        range.emplace_back(row,col,1);
        if(row%3==0){
            const int cell=row/3;
            // A spanning-tree edge is a local coordinate for the same
            // divergence range. Root-to-cell path sums give an equivalent
            // basis but make the physical mass matrix globally coupled.
            if(cell!=result.rootCell){
                range.emplace_back(3*parent[cell],col,-1);
                const auto face=domain.face(mesh,parentFace[cell]);const double sign=face.negativeCell==cell?1.:-1.;
                values.emplace_back(faces_[face.nativeRow].first,col,sign/face.area);
            }
            else if(!closed){const auto face=domain.face(mesh,result.outletFace);values.emplace_back(faces_[face.nativeRow].first,col,face.outwardSign/face.area);}
        }++col;
    }
    result.divergence.resize(np,nc);result.divergence.setFromTriplets(range.begin(),range.end());
    Matrix normal(nv,nc);normal.setFromTriplets(values.begin(),values.end());
    // The P1 divergence identities own the two cell velocity means. Correct
    // those linear tests locally after the graph supplies mean face fluxes.
    const Matrix residual=result.divergence-divergence_*normal;
    for(int col=0;col<residual.outerSize();col++)for(Matrix::InnerIterator it(residual,col);it;++it)if(it.row()%3){
        const int cell=int(it.row())/3,component=int(it.row())%3-1,row=firstCellVelocity(cell)+component;
        require(!fixed[row]&&divergence_.coeff(it.row(),row)!=0,"divergence complement lost its owned cell mean");
        values.emplace_back(row,col,it.value()/divergence_.coeff(it.row(),row));
    }
    result.moments.resize(nv,nc);result.moments.setFromTriplets(values.begin(),values.end());
    result.commutingResidual=(divergence_*result.moments-result.divergence).norm()/result.divergence.norm();
    require(result.commutingResidual<2e-12,"native divergence complement does not commute");
    for(int col=0;col<result.moments.outerSize();col++)for(Matrix::InnerIterator it(result.moments,col);it;++it)
        require(!fixed[it.row()]||it.value()==0,"divergence complement changes prescribed boundary data");
    result.pressureLoad=-Matrix(result.divergence.transpose());return result;
}
}
