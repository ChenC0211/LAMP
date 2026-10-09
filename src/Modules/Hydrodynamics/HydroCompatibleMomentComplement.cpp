#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineBasis.h"
#include <algorithm>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
std::array<double,2> bubble(double t){return {t*t*(1-t)*(1-t),2*t*(1-t)*(1-2*t)};}
std::array<double,2> end(double t){return {(1-t)*(1-t)*(1+2*t),-6*t*(1-t)};}
double weight(Vector2D p,bool rz){return rz?CartesianMeasure::twoPi*p.x:1.;}
void canonical(std::vector<double>& a){std::sort(a.begin(),a.end());a.erase(std::unique(a.begin(),a.end()),a.end());}
}
struct HydroCompatibleMomentComplement::Impl {
    struct Patch {Vector2D lower,upper,size;int component=0,side=-1;double scale=900;};
    std::shared_ptr<const NativeMeshIdentity> identity;
    CartesianLeafLocator locator;
    std::vector<Cell> cells;
    std::vector<std::vector<Patch>> patches;
    Matrix basisMap,moments,divergence,pressure;
    int nativeSize=0,root=0,outlet=-1;bool rz=false;double commuting=0;
    explicit Impl(const QuadtreeMesh& mesh):locator(mesh){}
    Eigen::MatrixXd jets(int cell,Vector2D p) const{
        const auto& info=cells.at(cell);const auto& list=patches.at(cell);Eigen::MatrixXd out=Eigen::MatrixXd::Zero(6,list.size());
        require(std::isfinite(p.x)&&std::isfinite(p.y)&&p.x>=info.partition[0].front()&&p.x<=info.partition[0].back()
            &&p.y>=info.partition[1].front()&&p.y<=info.partition[1].back(),"complement query outside its native cell");
        for(int k=0;k<int(list.size());k++){
            const auto& patch=list[k];
            if(p.x<patch.lower.x||p.x>patch.upper.x||p.y<patch.lower.y||p.y>patch.upper.y)continue;
            const double x=p.x==patch.upper.x?1.:(p.x-patch.lower.x)/patch.size.x,
                y=p.y==patch.upper.y?1.:(p.y-patch.lower.y)/patch.size.y;
            auto a=bubble(x),b=bubble(y);
            if(patch.side>=0){const double s=patch.component==0?x:y,t=patch.side?1-s:s;auto n=end(t);if(patch.side)n[1]=-n[1];if(patch.component==0)a=n;else b=n;}
            out(patch.component,k)=patch.scale*a[0]*b[0];
            out(2+patch.component,k)=patch.scale*a[1]*b[0]/patch.size.x;
            out(4+patch.component,k)=patch.scale*a[0]*b[1]/patch.size.y;
        }return out;
    }
};
HydroCompatibleMomentComplement HydroCompatibleMomentComplement::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const std::vector<bool>& fixed){
    space.validate(mesh,domain);const auto algebra=space.divergenceComplement(mesh,domain,fixed);
    auto p=std::make_shared<Impl>(mesh);p->identity=domain.meshIdentity();p->rz=domain.axisymmetric();p->nativeSize=int(space.mass().rows());p->root=algebra.rootCell;p->outlet=algebra.outletFace;
    const int n=mesh.cellIndexNum();p->cells.resize(n);p->patches.resize(n);
    std::vector<int> faceBasis(domain.nativeFaceCount(),-1),nativeFace(p->nativeSize,-1);
    for(const auto f:domain.faces(mesh))if(f.area>0)nativeFace[space.faces()[f.nativeRow].first]=f.nativeRow;
    std::vector<bool> used(domain.nativeFaceCount(),false);
    for(int col=0;col<algebra.moments.outerSize();col++)for(Matrix::InnerIterator it(algebra.moments,col);it;++it)
        if(it.value()!=0&&nativeFace[it.row()]>=0)used[nativeFace[it.row()]]=true;
    for(int cell=0;cell<n;cell++){
        const auto& c=mesh.cells(mesh.cellIndexToLabel[cell]);auto& info=p->cells[cell];
        info.lower=mesh.verts(c.cellToVert[0]).coord();const auto upper=mesh.verts(c.cellToVert[2]).coord();info.size=upper-info.lower;
        // Use the mesh endpoints shared by the solenoidal field and its
        // traces; recovering them from center +/- half-size changes identity.
        info.partition={std::vector<double>{info.lower.x,upper.x},std::vector<double>{info.lower.y,upper.y}};
    }
    int count=0;
    for(const auto f:domain.faces(mesh))if(used[f.nativeRow]){
        const int id=count++;faceBasis[f.nativeRow]=id;const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        for(int cell:{f.negativeCell,f.positiveCell})if(cell>=0){
            auto& info=p->cells[cell];Impl::Patch patch;patch.lower=info.lower;patch.upper={info.partition[0][1],info.partition[1][1]};
            patch.size=info.size;patch.component=f.horizontal?1:0;patch.scale=30;
            patch.side=cell==f.negativeCell?1:0;
            if(f.horizontal){patch.lower.x=a.x;patch.upper.x=b.x;patch.size.x=b.x-a.x;info.partition[0].push_back(a.x);info.partition[0].push_back(b.x);}
            else{patch.lower.y=a.y;patch.upper.y=b.y;patch.size.y=b.y-a.y;info.partition[1].push_back(a.y);info.partition[1].push_back(b.y);}
            info.basis.push_back(id);p->patches[cell].push_back(patch);
        }
    }
    std::vector<int> cellBasis(2*n);
    for(int cell=0;cell<n;cell++){
        auto& info=p->cells[cell];const Vector2D upper{info.partition[0][1],info.partition[1][1]};
        for(int d=0;d<2;d++){cellBasis[2*cell+d]=count;info.basis.push_back(count++);p->patches[cell].push_back({info.lower,upper,info.size,d,-1,900});}
        for(auto& partition:info.partition)canonical(partition);
    }
    std::vector<Eigen::Triplet<double>> moments;const auto& g=Geometry::Parametric::gaussEight();
    for(const auto f:domain.faces(mesh))if(faceBasis[f.nativeRow]>=0){
        const int cell=f.negativeCell>=0?f.negativeCell:f.positiveCell;const auto& info=p->cells[cell];const int id=faceBasis[f.nativeRow];
        const int local=int(std::find(info.basis.begin(),info.basis.end(),id)-info.basis.begin());
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const double length=f.horizontal?b.x-a.x:b.y-a.y;const auto& face=space.faces()[f.nativeRow];double mean=0,linear=0;
        for(int q=0;q<8;q++){
            const auto point=a+(b-a)*(.5+.5*g.x[q]);const double v=p->jets(cell,point)(f.horizontal?1:0,local),w=.5*length*g.w[q]*weight(point,p->rz)/face.area;
            mean+=w*v;linear+=w*v*((f.horizontal?point.x:point.y)-face.tangentCenter)/face.tangentScale;
        }
        moments.emplace_back(face.first,id,mean);moments.emplace_back(face.first+1,id,linear);
    }
    for(int cell=0;cell<n;cell++){
        const auto& info=p->cells[cell];const auto center=info.lower+info.size*.5;const double volume=domain.cellVolumes()[domain.cellRowByIndex()[cell]];
        Eigen::MatrixXd local=Eigen::MatrixXd::Zero(space.cellVelocityCount(),info.basis.size());
        const auto& xs=info.partition[0];const auto& ys=info.partition[1];
        for(std::size_t a=1;a<xs.size();a++)for(std::size_t b=1;b<ys.size();b++)for(int i=0;i<8;i++)for(int j=0;j<8;j++){
            const Vector2D point{xs[a-1]+.5*(xs[a]-xs[a-1])*(1+g.x[i]),ys[b-1]+.5*(ys[b]-ys[b-1])*(1+g.x[j])};
            const auto v=p->jets(cell,point);const double w=.25*(xs[a]-xs[a-1])*(ys[b]-ys[b-1])*g.w[i]*g.w[j]*weight(point,p->rz)/volume;
            local.row(0)+=w*v.row(0);local.row(1)+=w*v.row(1);local.row(2)+=w*(-(point.y-center.y)/info.size.y*v.row(0)+(point.x-center.x)/info.size.x*v.row(1));
        }
        if(p->rz)local.bottomRows(3)=space.rzExtraMoments(center,info.size,[&](Vector2D point)->Eigen::MatrixXd{return p->jets(cell,point).topRows(2);},info.partition);
        for(int row=0;row<local.rows();row++)for(int k=0;k<local.cols();k++)if(local(row,k)!=0)moments.emplace_back(space.firstCellVelocity(cell)+row,info.basis[k],local(row,k));
    }
    Matrix samples(p->nativeSize,count);samples.setFromTriplets(moments.begin(),moments.end());const Matrix d=space.divergence()*samples;
    std::vector<Eigen::Triplet<double>> coefficients;
    for(int col=0;col<algebra.moments.outerSize();col++)for(Matrix::InnerIterator it(algebra.moments,col);it;++it)
        if(it.value()!=0&&nativeFace[it.row()]>=0){const int face=nativeFace[it.row()],basis=faceBasis[face];coefficients.emplace_back(basis,col,it.value()/samples.coeff(it.row(),basis));}
    Matrix faceMap(count,algebra.divergence.cols());faceMap.setFromTriplets(coefficients.begin(),coefficients.end());
    const Matrix remaining=algebra.divergence-d*faceMap;
    for(int col=0;col<remaining.outerSize();col++)for(Matrix::InnerIterator it(remaining,col);it;++it)if(it.row()%3){
        const int basis=cellBasis[2*(it.row()/3)+it.row()%3-1];require(d.coeff(it.row(),basis)!=0,"smooth complement lost its cell divergence moment");
        coefficients.emplace_back(basis,col,it.value()/d.coeff(it.row(),basis));
    }
    p->basisMap.resize(count,algebra.divergence.cols());p->basisMap.setFromTriplets(coefficients.begin(),coefficients.end());p->moments=samples*p->basisMap;
    p->divergence=algebra.divergence;p->pressure=-Matrix(p->divergence.transpose());
    p->commuting=(space.divergence()*p->moments-p->divergence).norm()/p->divergence.norm();require(p->commuting<2e-12,"smooth physical complement does not reproduce native divergence");
    for(int col=0;col<p->moments.outerSize();col++)for(Matrix::InnerIterator it(p->moments,col);it;++it)require(!fixed[it.row()]||std::abs(it.value())<2e-12,"smooth complement changes prescribed boundary moments");
    HydroCompatibleMomentComplement result;result.impl_=std::move(p);return result;
}
void HydroCompatibleMomentComplement::validate(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentSpace& space) const{
    space.validate(mesh,domain);require(impl_&&impl_->identity==domain.meshIdentity()&&impl_->nativeSize==space.mass().rows(),"physical complement belongs to another native space");
}
const std::vector<HydroCompatibleMomentComplement::Cell>& HydroCompatibleMomentComplement::cells() const{require(bool(impl_),"moved physical complement");return impl_->cells;}
Eigen::MatrixXd HydroCompatibleMomentComplement::cellJets(int cell,Vector2D point) const{require(bool(impl_),"moved physical complement");return impl_->jets(cell,point);}
HydroCompatibleMomentSpline::NativeFaceTrace HydroCompatibleMomentComplement::faceTrace(int cell,
    const InterfaceNormalFluxView::FaceChart& face,bool horizontal) const{
    require(bool(impl_),"moved physical complement");using namespace MomentSplineBasis::Bernstein;
    const auto& patches=impl_->patches.at(cell);HydroCompatibleMomentSpline::NativeFaceTrace out;
    out.x=Eigen::MatrixXd::Zero(11,patches.size());out.y=out.x;
    const double lo=horizontal?face.begin.x:face.begin.y,hi=horizontal?face.end.x:face.end.y;
    require(hi>lo,"invalid physical complement face chart");
    for(int k=0;k<int(patches.size());k++){const auto& p=patches[k];
        if(p.side<0||p.component!=(horizontal?1:0))continue;
        const double normal=horizontal?face.begin.y:face.begin.x;
        const double fixed=horizontal?(p.side?p.upper.y:p.lower.y):(p.side?p.upper.x:p.lower.x);
        if(normal!=fixed)continue;
        const double begin=horizontal?p.lower.x:p.lower.y,finish=horizontal?p.upper.x:p.upper.y,size=horizontal?p.size.x:p.size.y;
        if(hi<=begin||lo>=finish)continue;
        require(lo>=begin&&hi<=finish,"native face crosses a physical complement support edge");
        const long double a=(lo-begin)/size,b=hi==finish?1.:(hi-begin)/size;
        Polynomial q{p.scale};q=linear(q,a,b);q=linear(q,a,b);q=linear(q,1-a,1-b);q=linear(q,1-a,1-b);
        if(impl_->rz&&horizontal&&lo>0)q=linear(q,1,hi/lo);
        q=elevate(std::move(q),10);for(int j=0;j<=10;j++)(horizontal?out.y:out.x)(j,k)=double(q[j]);
    }return out;
}
const HydroCompatibleMomentComplement::Matrix& HydroCompatibleMomentComplement::basisMap() const{return impl_->basisMap;}
const HydroCompatibleMomentComplement::Matrix& HydroCompatibleMomentComplement::nativeMoments() const{return impl_->moments;}
const HydroCompatibleMomentComplement::Matrix& HydroCompatibleMomentComplement::divergence() const{return impl_->divergence;}
const HydroCompatibleMomentComplement::Matrix& HydroCompatibleMomentComplement::pressureLoad() const{return impl_->pressure;}
double HydroCompatibleMomentComplement::commutingResidual() const{return impl_->commuting;}
int HydroCompatibleMomentComplement::rootCell() const{return impl_->root;}
int HydroCompatibleMomentComplement::outletFace() const{return impl_->outlet;}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentComplement::firstJet(Vector2D point,const Eigen::VectorXd& eta) const{
    require(impl_&&eta.size()==impl_->basisMap.cols()&&eta.allFinite(),"invalid physical complementary coordinates");
    const auto incident=impl_->locator.incident(point);require(incident.count>0,"physical complement query outside domain");const int cell=incident.indices[0];
    const Eigen::VectorXd coefficients=impl_->basisMap*eta;const auto& rows=impl_->cells[cell].basis;Eigen::VectorXd local(rows.size());for(int k=0;k<int(rows.size());k++)local[k]=coefficients[rows[k]];
    const Eigen::VectorXd v=impl_->jets(cell,point)*local;return {{v[0],v[1]},{v[2],v[3]},{v[4],v[5]},cell};
}
Eigen::VectorXd HydroCompatibleMomentComplement::transposePointForces(const std::vector<HydroCompatibleMomentSpline::PointForce>& forces) const{
    require(bool(impl_),"moved physical complement");Eigen::VectorXd local=Eigen::VectorXd::Zero(impl_->basisMap.rows());
    for(const auto& f:forces){const auto incident=impl_->locator.incident(f.point);require(incident.count>0,"complement work outside domain");
        const int cell=incident.indices[0];Eigen::Matrix<double,6,1> work;work<<f.measure.x,f.measure.y,f.dx.x,f.dx.y,f.dy.x,f.dy.y;
        require(work.allFinite(),"nonfinite complement work");const Eigen::VectorXd values=impl_->jets(cell,f.point).transpose()*work;
        const auto& rows=impl_->cells[cell].basis;for(int k=0;k<int(rows.size());k++)local[rows[k]]+=values[k];
    }return impl_->basisMap.transpose()*local;
}
}
