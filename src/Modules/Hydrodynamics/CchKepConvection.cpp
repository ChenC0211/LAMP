#include "Modules/Hydrodynamics/CchKepConvection.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/SVD>
#include <Eigen/SparseCholesky>
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
CchKineticHodge::Matrix CchKepConvection::modalBulkTransportMatrix(const QuadtreeMesh& mesh,
    const CchKineticHodge& hodge,const CchMaterialFlux& flux)
{
    hodge.validateFor(mesh,hodge.geometryIdentity(),hodge.timeLayer());
    flux.validateFor(mesh,hodge.materialMass().geometry(),hodge.geometryIdentity(),flux.advectorIdentity(),hodge.timeLayer());
    const auto& source=hodge.reconstruction();
    if(!source.sameIdentity(flux.phaseRate().spatialAdvector().reconstruction())
        ||flux.liquidDensity()!=hodge.materialMass().liquidDensity()||flux.gasDensity()!=hodge.materialMass().gasDensity())
        throw std::invalid_argument("[CchKepConvection] different modal source/material owner");
    const auto& regular=source.regularGeometry();const int nu=hodge.layout().velocityCount();
    const int n=nu+static_cast<int>(source.tangentialModes().size());
    std::vector<std::optional<Detail::CchRegularCellStencil>> stencils(hodge.layout().pressureCount());
    std::map<std::pair<int,int>,long double> entries;
    for(const auto& q:hodge.strainQuadrature()){
        const auto f=flux.physicalFlux(q.pressure,q.point,q.liquid);
        if(f.x==0&&f.y==0)continue;
        if(!stencils[q.pressure])stencils[q.pressure]=regular.cellStencil(mesh,q.pressure);
        struct Row{int column;Vector2D value,transport;};std::vector<Row> rows;
        for(const auto& w:regular.cellJetWeights(*stencils[q.pressure],q.point))
            rows.push_back({w.velocityDof,{w.jet[0],w.jet[1]},
                {f.x*w.jet[2]+f.y*w.jet[3],f.x*w.jet[4]+f.y*w.jet[5]}});
        const int mode=source.tangentialModeForPressure(q.pressure);
        if(mode>=0){const auto j=source.tangentialModeJet(mode,q.point,q.liquid);
            rows.push_back({nu+mode,j.value,{f.x*j.derivative.xx+f.y*j.derivative.xy,
                f.x*j.derivative.yx+f.y*j.derivative.yy}});}
        for(const auto& a:rows)for(const auto& b:rows)
            entries[{a.column,b.column}]+=q.physicalWeight*(static_cast<long double>(a.value.x)*b.transport.x
                +static_cast<long double>(a.value.y)*b.transport.y);
    }
    std::vector<Eigen::Triplet<double,Eigen::Index>> triplets;
    for(const auto& e:entries){const double value=static_cast<double>(e.second);
        if(!std::isfinite(value))throw std::runtime_error("[CchKepConvection] nonfinite modal bulk transport entry");
        if(value!=0)triplets.emplace_back(e.first.first,e.first.second,value);}
    CchKineticHodge::Matrix result(n,n);result.setFromTriplets(triplets.begin(),triplets.end());result.makeCompressed();return result;
}
namespace
{
constexpr long double tolerance=256.L*std::numeric_limits<double>::epsilon();
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[CchKepConvection] ")+why);}
long double finite(long double x,const char* why){require(std::isfinite(x),why);return x;}
long double product(long double a,long double b)
{
    const auto v=finite(a*b,"nonfinite convection product");
    require(a==0 || b==0 || v!=0,"unrepresentable convection product");return v;
}
long double identity(long double difference,long double scale,const char* why)
{
    finite(difference,why);finite(scale,why);require(scale>=0 && (scale>0 || difference==0),why);
    const auto error=scale>0?std::abs(difference)/scale:0;require(error<=tolerance,why);return error;
}
int component(const CchDofLayout& l,int u){return l.normalAxis(u)==CchDofLayout::NormalAxis::First?0:1;}
bool internal(const CchMomentumSegment& g){return g.negativeVelocity>=0 && g.positiveVelocity>=0;}
int owner(const CchMomentumSegment& g){return g.negativeVelocity>=0?g.negativeVelocity:g.positiveVelocity;}
Vector2D center(const CchMomentumSegment& g,bool rz)
{
    auto p=g.begin+(g.end-g.begin)*.5;
    if(rz && g.normalAxis==1){
        const long double a=g.begin.x,h=static_cast<long double>(g.end.x)-a,d=a+h/2;
        require(d>0,"invalid radial segment centroid");
        p.x=static_cast<double>(a+h*(a/(2*d)+h/(3*d)));
    }
    require(std::isfinite(p.x) && std::isfinite(p.y),"nonfinite segment centroid");return p;
}
void insert(std::vector<int>& values,int v)
{if(std::find(values.begin(),values.end(),v)==values.end())values.push_back(v);}
std::vector<CchConvectionWeight> reconstruct(const CchMomentumSegment& g,Vector2D point,
    const CchMomentumGeometry& geometry,const std::vector<NativeLeafFace>& faces,
    const std::vector<std::vector<int>>& cellFaces,const std::vector<std::vector<int>>& neighbors)
{
    const auto& l=geometry.layout();std::vector<int> roots,patch,support;
    for(int u:{g.negativeVelocity,g.positiveVelocity})for(const auto& c:geometry.volumes()[u].overlaps)insert(roots,c.pressure);
    patch=roots;for(int p:roots)for(int n:neighbors[p])insert(patch,n);
    for(int p:patch)for(int u:cellFaces[p])if(l.normalAxis(u)==l.normalAxis(g.negativeVelocity))insert(support,u);
    std::sort(support.begin(),support.end());require(!support.empty(),"empty common trace support");
    const auto& a=geometry.volumes()[g.negativeVelocity];const auto& b=geometry.volumes()[g.positiveVelocity];
    const Vector2D width{std::max(a.upper.x-a.lower.x,b.upper.x-b.lower.x),std::max(a.upper.y-a.lower.y,b.upper.y-b.lower.y)};
    require(width.x>0 && width.y>0,"invalid common trace scale");
    Eigen::MatrixXd constraints(3,support.size());Eigen::VectorXd w=Eigen::VectorXd::Zero(support.size());
    for(Eigen::Index j=0;j<constraints.cols();++j){const int u=support[static_cast<std::size_t>(j)];
        constraints(0,j)=1;constraints(1,j)=(faces[u].center.x-point.x)/width.x;constraints(2,j)=(faces[u].center.y-point.y)/width.y;
        if(u==g.negativeVelocity || u==g.positiveVelocity)w[j]=.5;
    }
    require(constraints.allFinite(),"nonfinite common trace constraints");
    Eigen::Vector3d target;target<<1,0,0;const Eigen::Vector3d residual=target-constraints*w;
    // Preserve ordinary center averaging exactly where it already reproduces
    // the physical common point. Else solve only this identifiable functional.
    if(residual.cwiseAbs().maxCoeff()>8*std::numeric_limits<double>::epsilon()){
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraints,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        w+=svd.solve(residual);
    }
    require(w.allFinite(),"nonfinite common trace weights");
    const Eigen::Vector3d difference=constraints*w-target;
    for(int k=0;k<3;++k){long double gross=std::abs(target[k]);for(Eigen::Index j=0;j<w.size();++j)gross+=std::abs(constraints(k,j)*w[j]);
        identity(difference[k],gross,"unidentifiable affine common trace");}
    std::vector<CchConvectionWeight> result;
    for(Eigen::Index j=0;j<w.size();++j)if(w[j]!=0)result.push_back({support[static_cast<std::size_t>(j)],w[j]});
    return result;
}
// Shared sparse native pairs, with block constraints on the FINAL operator.
// Local 1/6-row SVD blocks provide a minimum-change patch preconditioner;
// the transient multipliers are discarded after preparing the shared pairs.
struct MomentEntry {std::size_t edge;std::array<double,6> a{};};
struct MomentBlock
{
    int offset=0,count=0;std::vector<MomentEntry> entries;Eigen::MatrixXd whitening;
};
std::array<double,6> polynomial(Vector2D a,Vector2D b,double h)
{
    const double x=(b.x-a.x)/h,y=(b.y-a.y)/h;
    return {1,x,y,x*x,x*y,y*y};
}
void accumulate(double value,double& sum,double& compensation)
{const double term=value-compensation,next=sum+term;compensation=(next-sum)-term;sum=next;}
long double dot(const std::vector<double>& a,const std::vector<double>& b)
{double value=0,compensation=0;for(std::size_t i=0;i<a.size();++i)accumulate(a[i]*b[i],value,compensation);return finite(value,"nonfinite joint moment dot product");}
double maximum(const std::vector<double>& a)
{double value=0;for(double x:a)value=std::max(value,std::abs(x));return value;}
struct NativeTopology
{
    bool paired=false;Vector2D period{};
    std::vector<int> partner,representative,count;
    std::vector<std::vector<int>> roots;
};
Vector2D displacement(Vector2D a,Vector2D b,const NativeTopology& topology)
{
    auto d=b-a;
    if(topology.period.x>0)d.x-=topology.period.x*std::round(d.x/topology.period.x);
    if(topology.period.y>0)d.y-=topology.period.y*std::round(d.y/topology.period.y);
    return d;
}
std::array<double,6> polynomial(Vector2D a,Vector2D b,double h,const NativeTopology& topology)
{return polynomial({},displacement(a,b,topology),h);}
NativeTopology topologyFor(const CchDualMaterial& m,const std::vector<NativeLeafFace>& faces,
    std::vector<std::vector<int>>& neighbors,const std::vector<CchConvectionPeriodicPair>* pairs)
{
    NativeTopology t;t.paired=pairs!=nullptr;const int n=m.layout().velocityCount();const auto& segments=m.geometry().segments();
    t.partner.assign(segments.size(),-1);t.representative.resize(n);std::iota(t.representative.begin(),t.representative.end(),0);
    const auto root=[&](int i){while(t.representative[i]!=i)i=t.representative[i];return i;};
    if(pairs)for(const auto& p:*pairs){const int g=p.firstSegment,h=p.secondSegment;
        require(g>=0 && h>=0 && g!=h && static_cast<std::size_t>(g)<segments.size() && static_cast<std::size_t>(h)<segments.size(),"invalid periodic segment indices");
        const auto& a=segments[g];const auto& b=segments[h];require(!internal(a) && !internal(b) && t.partner[g]<0 && t.partner[h]<0,"repeated or internal periodic segment");
        const int i=owner(a),j=owner(b);require(component(m.layout(),i)==component(m.layout(),j) && a.normalAxis==b.normalAxis && (a.negativeVelocity>=0)!=(b.negativeVelocity>=0),"incompatible periodic segment orientation");
        require(std::isfinite(p.translation.x) && std::isfinite(p.translation.y),"nonfinite periodic translation");
        const double length=a.normalAxis==0?std::abs(p.translation.x):std::abs(p.translation.y);
        require(length>0 && (a.normalAxis==0?p.translation.y==0:p.translation.x==0),"periodic translation must be normal to segments");
        require(!m.layout().domain().axisymmetric() || a.normalAxis==1,"radial periodicity changes RZ measure");
        double& declared=a.normalAxis==0?t.period.x:t.period.y;if(declared>0)identity(length-declared,length+declared,"different periodic lengths");else declared=length;
        const auto begin=b.begin+p.translation,end=b.end+p.translation;
        identity(std::hypot(begin.x-a.begin.x,begin.y-a.begin.y),length,"periodic segment begin mismatch");
        identity(std::hypot(end.x-a.end.x,end.y-a.end.y),length,"periodic segment end mismatch");
        identity(a.area-b.area,a.area+b.area,"periodic segment measure mismatch");
        const auto fa=(a.negativeVelocity>=0?1:-1)*m.faceRates()[g],fb=(b.negativeVelocity>=0?1:-1)*m.faceRates()[h];
        identity(fa+fb,std::abs(fa)+std::abs(fb),"periodic material rates mismatch");t.partner[g]=h;t.partner[h]=g;
        if(component(m.layout(),i)==a.normalAxis){
            const auto d=displacement(faces[i].center,faces[j].center,t);identity(std::hypot(d.x,d.y),length,"periodic normal samples do not coincide");
            const int ri=root(i),rj=root(j);t.representative[std::max(ri,rj)]=std::min(ri,rj);
            // These actual primal boundary faces close pressure adjacency.
            for(const auto& x:m.geometry().volumes()[i].overlaps)for(const auto& y:m.geometry().volumes()[j].overlaps){insert(neighbors[x.pressure],y.pressure);insert(neighbors[y.pressure],x.pressure);}
        }
    }
    if(pairs)for(std::size_t g=0;g<segments.size();++g)if(!internal(segments[g]) && m.faceRates()[g]!=0)
        require(t.partner[g]>=0,"unpaired nonzero material boundary rate");
    t.roots.resize(n);t.count.resize(n);
    for(int i=0;i<n;++i){const int r=root(i);t.representative[i]=r;++t.count[r];for(const auto& p:m.geometry().volumes()[i].overlaps)insert(t.roots[r],p.pressure);}
    return t;
}
std::vector<CchConvectionCoupling> constrainOpaqueMoments(const CchDualMaterial& m,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,const std::vector<Vector2D>& points,
    const std::vector<std::vector<CchConvectionWeight>>& weights,CchConvectionMomentReport& report)
{
    const auto& l=m.layout();const int n=l.velocityCount();const auto& segments=m.geometry().segments();
    long double rateScale=0;for(auto f:m.faceRates())rateScale=std::max(rateScale,std::abs(f));
    if(rateScale==0)return {};
    using Pair=std::pair<int,int>;std::map<Pair,long double> original;
    const auto add=[&](int i,int j,long double value){if(i==j)return;if(i>j){std::swap(i,j);value=-value;}original[{i,j}]+=value;};
    std::vector<double> rowScale(n),h(n);std::vector<bool> exterior(n);
    std::vector<std::array<double,6>> rhs(n);
    for(int i=0;i<n;++i){const auto& box=m.geometry().volumes()[i];h[i]=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);require(h[i]>0,"invalid native moment scale");}
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const long double f=m.faceRates()[g]/rateScale;
        for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0)rowScale[i]+=static_cast<double>(std::abs(f));
        if(internal(s))for(const auto& w:weights[g]){add(s.negativeVelocity,w.velocity,.5L*f*w.weight);add(s.positiveVelocity,w.velocity,-.5L*f*w.weight);}
        else exterior[owner(s)]=true;
    }
    for(const auto& edge:original){rowScale[edge.first.first]+=static_cast<double>(std::abs(edge.second));rowScale[edge.first.second]+=static_cast<double>(std::abs(edge.second));}
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const double f=static_cast<double>(m.faceRates()[g]/rateScale);
        for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0 && !exterior[i]){const auto p=polynomial(faces[i].center,points[g],h[i]);const double sign=i==s.negativeVelocity?1:-1;
            rhs[i][1]+=sign*f*p[1];rhs[i][2]+=sign*f*p[2];}
    }
    for(const auto& edge:original){const int i=edge.first.first,j=edge.first.second;const double k=static_cast<double>(edge.second);
        for(int row:{i,j})if(!exterior[row]){const auto p=polynomial(faces[row].center,faces[row==i?j:i].center,h[row]);const double sign=row==i?1:-1;
            for(int c=1;c<6;++c)rhs[row][c]-=sign*k*p[c];}}
    // All pairs in the same geometry-only one-ring patch, shared globally.
    std::map<Pair,std::size_t> ids;
    for(int u=0;u<n;++u){std::vector<int> roots,patch,support;
        for(const auto& p:m.geometry().volumes()[u].overlaps)insert(roots,p.pressure);
        patch=roots;for(int p:roots)for(int q:neighbors[p])insert(patch,q);
        for(int p:patch)for(int v:cellFaces[p])if(component(l,v)==component(l,u))insert(support,v);
        std::sort(support.begin(),support.end());
        for(std::size_t a=0;a<support.size();++a)for(std::size_t b=a+1;b<support.size();++b){const int i=support[a],j=support[b];
            if(rowScale[i]>0 && rowScale[j]>0)ids.emplace(Pair{i,j},0);}
    }
    std::vector<Pair> pairs;std::vector<double> pairScale;
    for(auto& entry:ids){entry.second=pairs.size();pairs.push_back(entry.first);pairScale.push_back(std::sqrt(rowScale[entry.first.first]*rowScale[entry.first.second]));}
    std::vector<MomentBlock> blocks(n);std::vector<double> target;
    for(int i=0;i<n;++i){auto& b=blocks[i];b.offset=static_cast<int>(target.size());b.count=rowScale[i]>0?(exterior[i]?1:6):0;
        if(b.count==6)++report.interiorRows;for(int c=0;c<b.count;++c)target.push_back(rhs[i][c]/rowScale[i]);}
    report.constraints=static_cast<int>(target.size());
    for(std::size_t id=0;id<pairs.size();++id){const auto [i,j]=pairs[id];
        for(int row:{i,j}){MomentEntry e{id,polynomial(faces[row].center,faces[row==i?j:i].center,h[row])};const double scale=(row==i?1:-1)*pairScale[id]/rowScale[row];
            for(double& a:e.a)a*=scale;blocks[row].entries.push_back(e);}}
    const auto originalBlocks=blocks;const auto originalTarget=target;
    for(auto& b:blocks)if(b.count){Eigen::MatrixXd a(b.count,b.entries.size());for(Eigen::Index j=0;j<a.cols();++j)for(int c=0;c<b.count;++c)a(c,j)=b.entries[static_cast<std::size_t>(j)].a[c];
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        require(svd.rank()==b.count,"unidentifiable native moment patch");
        const Eigen::MatrixXd whitening=svd.singularValues().cwiseInverse().asDiagonal()*svd.matrixU().transpose();b.whitening=whitening;
        const Eigen::MatrixXd whitened=whitening*a;Eigen::VectorXd rhsBlock(b.count);for(int c=0;c<b.count;++c)rhsBlock[c]=target[b.offset+c];
        const Eigen::VectorXd whitenedTarget=whitening*rhsBlock;
        require(whitened.allFinite() && whitenedTarget.allFinite(),"nonfinite native moment preconditioner");
        for(int c=0;c<b.count;++c){target[b.offset+c]=whitenedTarget[c];for(std::size_t j=0;j<b.entries.size();++j)b.entries[j].a[c]=whitened(c,static_cast<Eigen::Index>(j));}}
    const auto transpose=[&](const std::vector<double>& y){std::vector<double> x(pairs.size()),compensation(x.size());for(const auto& b:blocks)for(const auto& e:b.entries){double value=0,local=0;for(int c=0;c<b.count;++c)accumulate(e.a[c]*y[b.offset+c],value,local);accumulate(value,x[e.edge],compensation[e.edge]);}return x;};
    const auto multiply=[&](const std::vector<double>& x){std::vector<double> y(target.size());for(const auto& b:blocks)for(int c=0;c<b.count;++c){double compensation=0;for(const auto& e:b.entries)accumulate(e.a[c]*x[e.edge],y[b.offset+c],compensation);}return y;};
    const auto precondition=[](const std::vector<double>& r){return r;};
    const auto residual=[&](const std::vector<double>& x){auto r=originalTarget;for(const auto& b:originalBlocks)for(int c=0;c<b.count;++c){double actual=0,compensation=0;for(const auto& e:b.entries)accumulate(e.a[c]*x[e.edge],actual,compensation);r[b.offset+c]-=actual;}return r;};
    const auto residualScale=[&](const std::vector<double>& x){std::vector<double> scale(originalTarget.size(),1);for(const auto& b:originalBlocks)for(int c=0;c<b.count;++c)for(const auto& e:b.entries)scale[b.offset+c]+=std::abs(e.a[c]*x[e.edge]);return scale;};
    const auto whiten=[&](const std::vector<double>& r){std::vector<double> y(r.size());for(const auto& b:blocks)for(int c=0;c<b.count;++c){double compensation=0;for(int k=0;k<b.count;++k)accumulate(b.whitening(c,k)*r[b.offset+k],y[b.offset+c],compensation);}return y;};
    std::vector<double> solution(pairs.size()),solutionCompensation(solution.size()),r=target,z=precondition(r),direction=z;
    long double rz=dot(r,z);const double stopping=8*std::numeric_limits<double>::epsilon()*std::max(1.,maximum(target));
    int iterations=0;const auto maximumIterations=std::max<std::size_t>(10000,4*target.size());
    while(static_cast<std::size_t>(iterations)<maximumIterations){
        if(maximum(r)<=stopping){const auto actual=residual(solution),scale=residualScale(solution);double error=0;for(std::size_t i=0;i<actual.size();++i)error=std::max(error,std::abs(actual[i])/scale[i]);if(error<=8*std::numeric_limits<double>::epsilon())break;
            r=whiten(actual);direction=precondition(r);rz=dot(r,direction);}
        const auto step=transpose(direction),a=multiply(step);const auto denominator=dot(direction,a);
        require(denominator>0 && rz>0,"incompatible native spatial/skew moments");const double alpha=static_cast<double>(rz/denominator);require(std::isfinite(alpha),"nonfinite native moment step");
        // Accumulate the minimum-norm pair solution directly: forming A^T
        // from a large final multiplier loses the small resolved moments.
        for(std::size_t i=0;i<solution.size();++i)accumulate(alpha*step[i],solution[i],solutionCompensation[i]);
        for(std::size_t i=0;i<r.size();++i)r[i]-=alpha*a[i];++iterations;
        // Retain conjugacy across the long-wavelength patch modes. Restart
        // only for explicit residual refinement, never on an arbitrary timer.
        const bool restart=maximum(r)<=stopping;
        if(restart)r=whiten(residual(solution));
        z=precondition(r);const auto next=dot(r,z);const double beta=restart?0:static_cast<double>(next/rz);
        for(std::size_t i=0;i<direction.size();++i)direction[i]=z[i]+beta*direction[i];rz=next;
    }
    report.iterations=iterations;
    const auto scale=residualScale(solution);
    for(const auto& b:originalBlocks)for(int c=0;c<b.count;++c){double actual=0,compensation=0;for(const auto& e:b.entries)accumulate(e.a[c]*solution[e.edge],actual,compensation);const auto error=std::abs(actual-originalTarget[b.offset+c]);
        const auto normalized=error/scale[b.offset+c];
        if(c==0)report.constantError=std::max(report.constantError,static_cast<long double>(normalized));else if(c<3)report.affineError=std::max(report.affineError,static_cast<long double>(normalized));else report.quadraticError=std::max(report.quadraticError,static_cast<long double>(normalized));}
    if(std::max({report.constantError,report.affineError,report.quadraticError})>tolerance)
        throw std::invalid_argument("[CchKepConvection] native moment solve did not converge: iterations="+std::to_string(iterations)+
            " scaled residual="+std::to_string(std::max({report.constantError,report.affineError,report.quadraticError})/tolerance)+" tolerances");
    std::vector<CchConvectionCoupling> result;
    for(std::size_t id=0;id<pairs.size();++id)if(solution[id]!=0){const auto value=product(rateScale,static_cast<long double>(pairScale[id])*solution[id]);result.push_back({pairs[id].first,pairs[id].second,value});}
    return result;
}
// A closed final skew operator: exact C/P1, Taylor P2 slabs and bounded
// coefficient norm. ADMM variables are temporary shared-pair algebra;
// no prepared failure can select the opaque operator.
// A shared-pair convex solve. Row normalization is numerical conditioning,
// never a change of physical moment or coefficient budgets.
using Sparse=Eigen::SparseMatrix<double>;
using Triplet=Eigen::Triplet<double>;
Eigen::VectorXd boundedProjection(const Sparse& raw,const Eigen::VectorXd& lower,
    const Eigen::VectorXd& upper,const Eigen::VectorXd& seed,CchConvectionMomentReport& report)
{
    Eigen::VectorXd scale=Eigen::VectorXd::Zero(raw.rows());
    for(int j=0;j<raw.outerSize();++j)for(Sparse::InnerIterator e(raw,j);e;++e)scale[e.row()]+=e.value()*e.value();
    for(int i=0;i<scale.size();++i){require(scale[i]>0,"unidentifiable bounded moment row");scale[i]=std::sqrt((lower[i]==upper[i]?1000.:1.)/scale[i]);}
    const Sparse a=scale.asDiagonal()*raw;const Eigen::VectorXd lo=scale.cwiseProduct(lower),hi=scale.cwiseProduct(upper);
    const auto clip=[&](const Eigen::VectorXd& x)->Eigen::VectorXd{return x.cwiseMax(lo).cwiseMin(hi);};
    const Sparse transpose=a.transpose();const double rho=10;
    Sparse normal=a*transpose;for(int i=0;i<normal.rows();++i)normal.coeffRef(i,i)+=1/rho;normal.makeCompressed();
    Eigen::SimplicialLDLT<Sparse> factor;factor.compute(normal);require(factor.info()==Eigen::Success,"bounded moment factorization failed");
    const Eigen::VectorXd origin=a*seed;Eigen::VectorXd x=seed,z=clip(origin),y=Eigen::VectorXd::Zero(a.rows());
    bool converged=false;
    for(int iteration=0;iteration<20000;++iteration){
        const Eigen::VectorXd multiplier=factor.solve(z-y-origin);require(factor.info()==Eigen::Success && multiplier.allFinite(),"nonfinite bounded moment solve");
        x=seed+transpose*multiplier;const Eigen::VectorXd value=a*x,previous=z,relaxed=1.6*value-.6*previous;
        z=clip(relaxed+y);y+=relaxed-z;++report.iterations;
        if(iteration%25==0){const double primal=(value-z).lpNorm<Eigen::Infinity>(),dual=(rho*transpose*(z-previous)).lpNorm<Eigen::Infinity>();
            if(std::max(primal,dual)<=1e-10){converged=true;break;}}
    }
    require(converged,"bounded moment ADMM did not converge");
    // Polish the identified active constraints jointly. Regularization only
    // stabilizes a rank-deficient factorization; residual refinement removes it.
    const Eigen::VectorXd value=a*x;std::vector<int> active(a.rows(),-1);std::vector<double> targets;
    for(int i=0;i<a.rows();++i){double target=0;bool use=false;
        if(lo[i]==hi[i]){use=true;target=lo[i];}
        else if(std::isfinite(lo[i]) && std::abs(value[i]-lo[i])<1e-8){use=true;target=lo[i];}
        else if(std::isfinite(hi[i]) && std::abs(value[i]-hi[i])<1e-8){use=true;target=hi[i];}
        if(use){active[i]=static_cast<int>(targets.size());targets.push_back(target);}}
    std::vector<Triplet> entries;for(int j=0;j<a.outerSize();++j)for(Sparse::InnerIterator e(a,j);e;++e)if(active[e.row()]>=0)entries.emplace_back(active[e.row()],j,e.value());
    Sparse hard(targets.size(),a.cols());hard.setFromTriplets(entries.begin(),entries.end());const Sparse ht=hard.transpose();
    Sparse schur=hard*ht;for(int i=0;i<schur.rows();++i)schur.coeffRef(i,i)+=1e-9;schur.makeCompressed();
    Eigen::SimplicialLDLT<Sparse> polish;polish.compute(schur);require(polish.info()==Eigen::Success,"bounded active moment factorization failed");
    const Eigen::Map<const Eigen::VectorXd> target(targets.data(),targets.size());
    for(int pass=0;pass<32;++pass){const Eigen::VectorXd residual=target-hard*x;if(residual.lpNorm<Eigen::Infinity>()<=8*std::numeric_limits<double>::epsilon())break;
        const Eigen::VectorXd dual=polish.solve(residual);require(dual.allFinite(),"nonfinite bounded active refinement");x+=ht*dual;++report.projectionCalls;}
    require(x.allFinite(),"nonfinite final bounded coefficients");return x;
}
std::vector<CchConvectionCoupling> constrainBoundedMoments(const CchDualMaterial& m,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,const std::vector<Vector2D>&,
    const std::vector<std::vector<CchConvectionWeight>>& weights,const NativeTopology& topology,CchConvectionMomentReport& report)
{
    const int n=m.layout().velocityCount();const auto& segments=m.geometry().segments();const auto& actual=m.material().midpointRate().spatialAdvector();
    long double rateScale=0;for(auto f:m.faceRates())rateScale=std::max(rateScale,std::abs(f));if(rateScale==0)return {};
    using Pair=std::pair<int,int>;std::map<Pair,long double> original,central;
    const auto add=[&](std::map<Pair,long double>& map,int i,int j,long double k){i=topology.representative[i];j=topology.representative[j];if(i==j)return;if(i>j){std::swap(i,j);k=-k;}map[{i,j}]+=k;};
    std::vector<double> h(n),gross(n);std::vector<std::array<double,3>> target(n);std::vector<int> row(n,-1);
    Vector2D lower=faces.front().center,upper=lower;int count=0;
    for(int i=0;i<n;++i){const int r=topology.representative[i];const auto& box=m.geometry().volumes()[i];h[r]=std::max(h[r],std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y));
        lower.x=std::min(lower.x,faces[i].center.x);lower.y=std::min(lower.y,faces[i].center.y);upper.x=std::max(upper.x,faces[i].center.x);upper.y=std::max(upper.y,faces[i].center.y);
        const auto b=actual.nativeVelocity(i);const auto mass=m.midpointMass().masses()[i]/rateScale;target[r][1]+=static_cast<double>(mass*b.x);target[r][2]+=static_cast<double>(mass*b.y);}
    const double length=std::max({upper.x-lower.x,upper.y-lower.y,topology.period.x,topology.period.y});require(length>0,"invalid bounded domain scale");
    for(std::size_t g=0;g<segments.size();++g){const auto& s=segments[g];const long double f=m.faceRates()[g]/rateScale;
        for(int i:{s.negativeVelocity,s.positiveVelocity})if(i>=0){const int r=topology.representative[i];gross[r]+=static_cast<double>(std::abs(f));target[r][0]+=.5*(i==s.negativeVelocity?1:-1)*static_cast<double>(f);}
        if(internal(s)){add(central,s.negativeVelocity,s.positiveVelocity,.5L*f);for(const auto& w:weights[g]){add(original,s.negativeVelocity,w.velocity,.5L*f*w.weight);add(original,s.positiveVelocity,w.velocity,-.5L*f*w.weight);}}
        else if(static_cast<int>(g)<topology.partner[g]){const auto k=.5L*(s.negativeVelocity>=0?1:-1)*f;const int j=owner(segments[topology.partner[g]]);add(original,owner(s),j,k);add(central,owner(s),j,k);}}
    for(int i=0;i<n;++i)if(topology.representative[i]==i){require(gross[i]>0 && h[i]>0,"zero-scale bounded native row");row[i]=count++;target[i][1]/=h[i];target[i][2]/=h[i];}
    std::map<Pair,int> ids;
    for(int u=0;u<n;++u)if(topology.representative[u]==u){const auto& roots=topology.roots[u];auto patch=roots;std::vector<int> support;
        for(int p:roots)for(int q:neighbors[p])insert(patch,q);for(int p:patch)for(int v:cellFaces[p])if(component(m.layout(),v)==component(m.layout(),u))insert(support,topology.representative[v]);
        std::sort(support.begin(),support.end());for(std::size_t a=0;a<support.size();++a)for(std::size_t b=a+1;b<support.size();++b)ids.emplace(Pair{support[a],support[b]},0);}
    std::vector<Pair> pairs;std::vector<double> pairScale;std::vector<std::vector<std::pair<int,double>>> incident(n);std::vector<Triplet> entries;
    for(auto& e:ids){const int k=static_cast<int>(pairs.size());e.second=k;pairs.push_back(e.first);const auto [i,j]=e.first;const double scale=std::sqrt(gross[i]*gross[j]);pairScale.push_back(scale);
        for(int r:{i,j}){const auto p=polynomial(faces[r].center,faces[r==i?j:i].center,h[r],topology);const double sign=r==i?1:-1;incident[r].push_back({k,scale/gross[r]});
            report.maximumSupportRatio=std::max(report.maximumSupportRatio,static_cast<long double>(std::hypot(p[1],p[2])));
            for(int c=0;c<6;++c)if(p[c]!=0)entries.emplace_back(6*row[r]+c,k,sign*scale/gross[r]*p[c]);}}
    for(const auto& e:original)require(e.second==0 || ids.count(e.first)==1,"original pair outside declared one-ring");
    for(int i=0;i<n;++i)if(row[i]>=0){Eigen::MatrixXd local(3,incident[i].size());
        for(std::size_t j=0;j<incident[i].size();++j){const auto& edge=pairs[incident[i][j].first];const int other=edge.first==i?edge.second:edge.first;const auto p=polynomial(faces[i].center,faces[other].center,h[i],topology);for(int c=0;c<3;++c)local(c,j)=p[c];}
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(local);svd.setThreshold(1e-12);require(svd.rank()==3,"unidentifiable bounded affine patch");}
    Eigen::VectorXd seed(pairs.size());for(int j=0;j<seed.size();++j)seed[j]=static_cast<double>(central[pairs[j]])/pairScale[j];
    std::vector<double> lo,hi;for(int i=0;i<n;++i)if(row[i]>=0)for(int c=0;c<6;++c){lo.push_back(c<3?target[i][c]/gross[i]:-h[i]/length);hi.push_back(c<3?target[i][c]/gross[i]:h[i]/length);}
    Sparse matrix(lo.size(),pairs.size());matrix.setFromTriplets(entries.begin(),entries.end());Eigen::VectorXd solution=seed;
    // L1 balls are imposed by their supporting half-spaces. Add only violated
    // facets, always solving the same final-pair QP; no alternate operator.
    bool bounded=false;
    for(int pass=0;pass<32;++pass){
        const Eigen::Map<const Eigen::VectorXd> lowerBound(lo.data(),lo.size()),upperBound(hi.data(),hi.size());
        const Eigen::VectorXd value=matrix*solution;bool feasible=true;
        for(int i=0;i<value.size();++i)if(value[i]<lo[i]-tolerance || value[i]>hi[i]+tolerance)feasible=false;
        if(!feasible)solution=boundedProjection(matrix,lowerBound,upperBound,seed,report);
        bounded=true;for(int i=0;i<n;++i)if(row[i]>=0){double norm=0;for(const auto& e:incident[i])norm+=e.second*std::abs(solution[e.first]);
            if(norm>2+tolerance){bounded=false;const int r=static_cast<int>(lo.size());lo.push_back(-std::numeric_limits<double>::infinity());hi.push_back(2);
                for(const auto& e:incident[i])entries.emplace_back(r,e.first,std::copysign(e.second,solution[e.first]));}}
        if(bounded)break;matrix.resize(lo.size(),pairs.size());matrix.setFromTriplets(entries.begin(),entries.end());}
    require(bounded,"bounded coefficient facets did not converge");
    const Eigen::VectorXd value=matrix*solution;report.boundedRows=report.interiorRows=count;report.constraints=static_cast<int>(lo.size());
    for(int i=0;i<n;++i)if(row[i]>=0){double norm=0;for(const auto& e:incident[i])norm+=e.second*std::abs(solution[e.first]);report.coefficientRatio=std::max(report.coefficientRatio,static_cast<long double>(norm/2));
        for(int c=0;c<6;++c){const int r=6*row[i]+c;if(c<3){const auto error=std::abs(value[r]-lo[r]);if(c==0)report.constantError=std::max(report.constantError,static_cast<long double>(error));else report.affineError=std::max(report.affineError,static_cast<long double>(error));}
            else{report.quadraticBudgetRatio=std::max(report.quadraticBudgetRatio,static_cast<long double>(std::abs(value[r])/hi[r]));report.quadraticBoundError=std::max(report.quadraticBoundError,static_cast<long double>(std::max(0.,std::abs(value[r])-hi[r])));}}}
    require(std::max({report.constantError,report.affineError,report.quadraticBoundError})<=tolerance,"final bounded native moments failed");
    std::vector<CchConvectionCoupling> result;for(int j=0;j<solution.size();++j){const auto correction=static_cast<long double>(solution[j])*pairScale[j]-original[pairs[j]];if(correction!=0)result.push_back({pairs[j].first,pairs[j].second,product(rateScale,correction)});}return result;
}

}
struct DissipationPatch
{
    std::vector<int> support;std::vector<double> squareRootMassWeights;
    Eigen::MatrixXd image;
};
struct DissipationConnection
{
    int first=-1,second=-1,firstSegment=-1,secondSegment=-1;
};
namespace
{
void validateDissipation(const CchConvectionDissipationConfig& c)
{
    require(std::isfinite(c.strength) && c.strength>=0 && c.strength<=1,
        "dissipation strength must be finite and in [0,1]");
    require(std::isfinite(c.activation) && std::isfinite(c.saturation) &&
        c.activation>=0 && c.activation<c.saturation,
        "dissipation thresholds must be finite with 0 <= activation < saturation");
}
void prepareDissipation(const CchDualMaterial& m,const NativeTopology& topology,
    const std::vector<NativeLeafFace>& faces,const std::vector<std::vector<int>>& cellFaces,
    const std::vector<std::vector<int>>& neighbors,std::vector<DissipationPatch>& patches,
    std::vector<DissipationConnection>& connections)
{
    const int n=m.layout().velocityCount();std::vector<long double> mass(n);
    std::vector<std::vector<int>> supports(n);
    for(int i=0;i<n;++i)mass[topology.representative[i]]+=m.midpointMass().masses()[i];
    for(int i=0;i<n;++i)if(topology.representative[i]==i){
        auto cells=topology.roots[i];
        for(int c:topology.roots[i])for(int j:neighbors[c])insert(cells,j);
        for(int c:cells)for(int j:cellFaces[c])if(component(m.layout(),i)==component(m.layout(),j))
            insert(supports[i],topology.representative[j]);
    }
    const auto& segments=m.geometry().segments();
    for(std::size_t g=0;g<segments.size();++g){
        const auto& s=segments[g];int i=-1,j=-1,other=-1;
        if(internal(s)){i=s.negativeVelocity;j=s.positiveVelocity;}
        else if(topology.partner[g]>static_cast<int>(g)){other=topology.partner[g];i=owner(s);j=owner(segments[other]);}
        else continue;
        const int ri=topology.representative[i],rj=topology.representative[j];
        if(ri==rj)continue; // Identical periodic samples have exactly zero jump.
        DissipationPatch p;p.support=supports[ri];for(int k:supports[rj])insert(p.support,k);
        std::sort(p.support.begin(),p.support.end());require(!p.support.empty(),"empty dissipation sensor support");
        const auto& box=m.geometry().volumes()[i];const double h=std::max(box.upper.x-box.lower.x,box.upper.y-box.lower.y);
        const auto point=faces[i].center+displacement(faces[i].center,faces[j].center,topology)*.5;
        long double largest=0,total=0;for(int j:p.support)largest=std::max(largest,mass[j]);
        require(std::isfinite(largest) && largest>0 && h>0,"invalid dissipation sensor mass or scale");
        for(int j:p.support)total+=mass[j]/largest;
        Eigen::MatrixXd a(p.support.size(),6);
        for(std::size_t k=0;k<p.support.size();++k){
            const double w=static_cast<double>((mass[p.support[k]]/largest)/total);
            require(w>0 && std::isfinite(w),"unrepresentable dissipation sensor mass weight");p.squareRootMassWeights.push_back(std::sqrt(w));
            const auto q=polynomial(point,faces[p.support[k]].center,h,topology);
            for(int c=0;c<6;++c)a(k,c)=p.squareRootMassWeights.back()*q[c];
        }
        // A sampling projection for the sensor only. No constraint or change
        // to the final convection operator, material flux, or actual B/J.
        require(a.allFinite(),"nonfinite dissipation sensor geometry");
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);
        p.image=svd.matrixU().leftCols(svd.rank());require(p.image.allFinite(),"nonfinite dissipation sensor projection");
        patches.push_back(std::move(p));connections.push_back({i,j,static_cast<int>(g),other});
    }
}
long double sensorValue(const DissipationPatch& p,const std::vector<double>& u,std::vector<long double>& y)
{
    const int root=p.support.front();long double scale=0;
    for(int j:p.support)scale=std::max(scale,std::abs(.5L*u[j]-.5L*u[root]));
    if(scale==0)return 0;
    std::array<long double,6> coefficient{};long double mean=0;y.resize(p.support.size());
    for(std::size_t k=0;k<p.support.size();++k){
        y[k]=p.squareRootMassWeights[k]*(.5L*u[p.support[k]]-.5L*u[root])/scale;
        mean+=p.squareRootMassWeights[k]*y[k];
        for(int c=0;c<p.image.cols();++c)coefficient[c]+=p.image(k,c)*y[k];
    }
    long double variance=0,residual=0;
    for(std::size_t k=0;k<p.support.size();++k){
        const auto deviation=y[k]-p.squareRootMassWeights[k]*mean;variance+=deviation*deviation;
        auto r=y[k];for(int c=0;c<p.image.cols();++c)r-=p.image(k,c)*coefficient[c];residual+=r*r;
    }
    return variance>0?finite(std::sqrt(residual/variance),"nonfinite Favre dissipation sensor"):0;
}
}
struct CchConvectionBoundary::Data
{
    Data(CchDualMaterial m,std::vector<double> t):material(std::move(m)),traces(std::move(t)){}
    CchDualMaterial material;std::vector<double> traces;
};
CchConvectionBoundary::CchConvectionBoundary(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchConvectionBoundary CchConvectionBoundary::bind(const QuadtreeMesh& mesh,CchDualMaterial m,
    const std::vector<double>& traces,CchMaterialSources sources,Interface::HermiteGalsTimeLayer clock)
{
    m.validateFor(mesh,sources,clock);const auto& segments=m.geometry().segments();
    require(traces.size()==segments.size(),"boundary trace dimensions do not match shared segments");
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];require(std::isfinite(traces[id]),"nonfinite boundary trace");
        if(internal(g))require(traces[id]==0,"nonzero unused internal boundary trace");
        else if(m.layout().domain().axisymmetric() && g.normalAxis==0 && g.begin.x==0 && component(m.layout(),owner(g))==0)
            require(traces[id]==0,"nonzero radial axis boundary trace");
    }
    return CchConvectionBoundary(std::make_shared<Data>(std::move(m),traces));
}
const CchDualMaterial& CchConvectionBoundary::material() const noexcept{return data_->material;}
const std::vector<double>& CchConvectionBoundary::traces() const noexcept{return data_->traces;}
void CchConvectionBoundary::validateFor(const QuadtreeMesh& mesh,const CchDualMaterial& expected,
    CchMaterialSources sources,Interface::HermiteGalsTimeLayer clock) const
{
    require(&material().faceRates()==&expected.faceRates(),"different prepared material owner");
    material().validateFor(mesh,sources,clock);
}
struct CchKepConvection::Data
{
    explicit Data(CchDualMaterial m):material(std::move(m)){}
    CchDualMaterial material;std::vector<Vector2D> points;std::vector<std::vector<CchConvectionWeight>> weights;
    std::vector<CchConvectionCoupling> couplings;CchConvectionMomentReport momentReport;
    std::vector<CchConvectionPeriodicPair> periodicPairs;NativeTopology topology;std::vector<long double> outflow;
    CchConvectionDissipationConfig dissipation;
    std::vector<DissipationPatch> sensorPatches;std::vector<DissipationConnection> dissipationConnections;
};
CchKepConvection::CchKepConvection(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchKepConvection CchKepConvection::build(const QuadtreeMesh& mesh,CchDualMaterial m)
{return buildWithPorts(mesh,std::move(m),nullptr);}
CchKepConvection CchKepConvection::build(const QuadtreeMesh& mesh,CchDualMaterial m,const std::vector<CchConvectionPeriodicPair>& pairs)
{return buildWithPorts(mesh,std::move(m),&pairs);}
CchKepConvection CchKepConvection::buildControlled(const QuadtreeMesh& mesh,CchDualMaterial m,CchConvectionDissipationConfig config)
{validateDissipation(config);return buildWithPorts(mesh,std::move(m),nullptr,&config);}
CchKepConvection CchKepConvection::buildControlled(const QuadtreeMesh& mesh,CchDualMaterial m,const std::vector<CchConvectionPeriodicPair>& pairs,CchConvectionDissipationConfig config)
{validateDissipation(config);return buildWithPorts(mesh,std::move(m),&pairs,&config);}
CchKepConvection CchKepConvection::buildWithPorts(const QuadtreeMesh& mesh,CchDualMaterial m,const std::vector<CchConvectionPeriodicPair>* pairs,const CchConvectionDissipationConfig* config)
{
    m.validateFor(mesh,m.sources(),m.timeLayer());auto d=std::make_shared<Data>(std::move(m));const auto& l=d->material.layout();
    std::vector<NativeLeafFace> faces;std::vector<std::vector<int>> cellFaces(l.pressureCount()),neighbors(l.pressureCount());
    for(int u=0;u<l.velocityCount();++u){const auto f=l.domain().face(mesh,l.nativeFace(u));faces.push_back(f);
        const int a=f.negativeCell<0?-1:l.pressureDof(f.negativeCell),b=f.positiveCell<0?-1:l.pressureDof(f.positiveCell);
        if(a>=0)cellFaces[a].push_back(u);if(b>=0)cellFaces[b].push_back(u);
        if(a>=0 && b>=0){insert(neighbors[a],b);insert(neighbors[b],a);}
    }
    for(const auto& g:d->material.geometry().segments()){
        const auto p=center(g,l.domain().axisymmetric());d->points.push_back(p);
        d->weights.push_back(internal(g)?reconstruct(g,p,d->material.geometry(),faces,cellFaces,neighbors):std::vector<CchConvectionWeight>{});
    }
    d->topology=topologyFor(d->material,faces,neighbors,pairs);if(pairs)d->periodicPairs=*pairs;
    d->couplings=d->topology.paired?constrainBoundedMoments(d->material,faces,cellFaces,neighbors,d->points,d->weights,d->topology,d->momentReport):
        constrainOpaqueMoments(d->material,faces,cellFaces,neighbors,d->points,d->weights,d->momentReport);
    d->outflow.resize(l.velocityCount());for(std::size_t g=0;g<d->material.faceRates().size();++g){const auto& s=d->material.geometry().segments()[g];const auto f=d->material.faceRates()[g];if(s.negativeVelocity>=0)d->outflow[s.negativeVelocity]+=f;if(s.positiveVelocity>=0)d->outflow[s.positiveVelocity]-=f;}
    if(config && config->strength>0){d->dissipation=*config;prepareDissipation(d->material,d->topology,faces,cellFaces,neighbors,d->sensorPatches,d->dissipationConnections);}
    return CchKepConvection(std::move(d));
}
const CchDualMaterial& CchKepConvection::material() const noexcept{return data_->material;}
const CchDofLayout& CchKepConvection::layout() const noexcept{return material().layout();}
const std::vector<Vector2D>& CchKepConvection::points() const noexcept{return data_->points;}
const std::vector<std::vector<CchConvectionWeight>>& CchKepConvection::weights() const noexcept{return data_->weights;}
const std::vector<CchConvectionCoupling>& CchKepConvection::couplings() const noexcept{return data_->couplings;}
const CchConvectionMomentReport& CchKepConvection::momentReport() const noexcept{return data_->momentReport;}
const std::vector<CchConvectionPeriodicPair>& CchKepConvection::periodicPairs() const noexcept{return data_->periodicPairs;}
const std::vector<int>& CchKepConvection::representatives() const noexcept{return data_->topology.representative;}
struct CchKepConvectionResult::Data
{
    Data(CchDualMaterial m,CchConvectionBoundary b):material(std::move(m)),boundary(std::move(b)){}
    CchDualMaterial material;CchConvectionBoundary boundary;std::vector<long double> center,rates,loads,dissipation;CchKepConvectionReport report;
    std::vector<CchConvectionDissipationSample> dissipationSamples;
};
CchKepConvectionResult::CchKepConvectionResult(std::shared_ptr<const Data> d):data_(std::move(d)){}
CchKepConvectionResult CchKepConvection::apply(const QuadtreeMesh& mesh,const CchNativeValues& values,
    const CchConvectionBoundary& boundary,CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{
    material().validateFor(mesh,source,clock);values.validateFor(layout(),mesh);boundary.validateFor(mesh,material(),source,clock);
    const auto impulse=material().consistencyImpulse(mesh,values,source,clock);const auto& u=values.normalVelocity();
    if(data_->topology.paired){
        for(int i=0;i<layout().velocityCount();++i){const int r=representatives()[i];if(i!=r)identity(u[i]-static_cast<long double>(u[r]),std::abs(static_cast<long double>(u[i]))+std::abs(static_cast<long double>(u[r])),"different periodic alias velocities");}
        for(const auto& p:periodicPairs()){const int i=owner(material().geometry().segments()[p.firstSegment]),j=owner(material().geometry().segments()[p.secondSegment]);const auto trace=.5L*u[i]+.5L*u[j];
            for(int g:{p.firstSegment,p.secondSegment})identity(boundary.traces()[g]-trace,std::abs(trace)+std::abs(static_cast<long double>(boundary.traces()[g])),"different periodic port trace");}
    }
    auto d=std::make_shared<CchKepConvectionResult::Data>(material(),boundary);d->rates.resize(u.size());d->loads.resize(u.size());
    std::array<long double,2> momentumScale{},powerScale{};const auto& segments=material().geometry().segments();
    for(std::size_t id=0;id<segments.size();++id){const auto& g=segments[id];const auto f=material().faceRates()[id];long double trace=0;
        if(internal(g))for(const auto& w:weights()[id])trace+=product(w.weight,u[w.velocity]);else trace=boundary.traces()[id];
        const auto raw=product(f,finite(trace,"nonfinite common velocity trace"));d->center.push_back(raw);
        if(internal(g)){
            const int i=g.negativeVelocity,j=g.positiveVelocity,a=component(layout(),i);const auto halfJump=.5L*u[i]-.5L*u[j];
            const auto left=product(f,.5L*trace+.5L*u[i]),right=product(f,.5L*trace+.5L*u[j]);d->rates[i]+=left;d->rates[j]-=right;
            momentumScale[a]+=std::abs(left)+std::abs(right);powerScale[a]+=std::abs(product(u[i],left))+std::abs(product(u[j],right));
            const auto factor=product(f,halfJump);for(const auto& w:weights()[id]){const auto correction=product(factor,w.weight);d->rates[w.velocity]-=correction;
                momentumScale[a]+=std::abs(correction);powerScale[a]+=std::abs(product(u[w.velocity],correction));}
        }else{
            const int i=owner(g),a=component(layout(),i);const long double sign=g.negativeVelocity>=0?1:-1;
            const auto flux=sign*raw;d->rates[i]+=flux;d->report.components[a].boundaryMomentumRate+=flux;
            // Multiply by the mass rate before the second velocity, avoiding
            // an unrepresentable u*u intermediate when the physical power fits.
            const auto work=product(u[i],flux-product(sign*f,.5L*u[i]));d->report.components[a].boundaryKineticPower+=work;
            momentumScale[a]+=std::abs(flux);powerScale[a]+=std::abs(work);
        }
    }
    for(const auto& e:couplings()){
        const int i=e.first,j=e.second,a=component(layout(),i);const auto left=product(e.coefficient,u[j]),right=product(e.coefficient,u[i]);
        d->rates[i]+=left;d->rates[j]-=right;momentumScale[a]+=std::abs(left)+std::abs(right);
        powerScale[a]+=std::abs(product(u[i],left))+std::abs(product(u[j],right));
    }
    d->dissipation.resize(u.size());
    if(!data_->sensorPatches.empty()){
        std::vector<long double> scratch;
        const auto& config=data_->dissipation;
        for(std::size_t k=0;k<data_->dissipationConnections.size();++k){const auto& e=data_->dissipationConnections[k];
            const int i=e.first,j=e.second,a=component(layout(),i);
            CchConvectionDissipationSample sample;sample.first=i;sample.second=j;sample.firstSegment=e.firstSegment;sample.secondSegment=e.secondSegment;
            sample.sensor=sensorValue(data_->sensorPatches[k],u,scratch);
            const auto ramp=sample.sensor<=config.activation?0.L:sample.sensor>=config.saturation?1.L:
                (sample.sensor-config.activation)/(config.saturation-config.activation);
            sample.chi=config.strength*ramp*ramp;sample.halfJump=.5L*u[i]-.5L*u[j];
            sample.momentumRate=product(product(sample.chi,std::abs(material().faceRates()[e.firstSegment])),sample.halfJump);
            sample.power=finite(2*product(sample.halfJump,sample.momentumRate),"nonfinite dissipation power");
            require(sample.power>=0,"negative dissipation power");
            d->dissipation[i]+=sample.momentumRate;d->dissipation[j]-=sample.momentumRate;
            d->rates[i]+=sample.momentumRate;d->rates[j]-=sample.momentumRate;
            auto& r=d->report.components[a];r.dissipationPower+=sample.power;r.maximumSensor=std::max(r.maximumSensor,sample.sensor);
            ++r.dissipationConnections;if(sample.chi>0 && material().faceRates()[e.firstSegment]!=0)++r.activeDissipationConnections;
            momentumScale[a]+=2*std::abs(sample.momentumRate);
            powerScale[a]+=std::abs(product(u[i],sample.momentumRate))+std::abs(product(u[j],sample.momentumRate))+sample.power;
            d->dissipationSamples.push_back(sample);
        }
    }
    std::vector<long double> acceleration(u.size()),groupMass(u.size()),groupDissipation(u.size());
    if(data_->topology.paired){
        for(int i=0;i<layout().velocityCount();++i){const int r=representatives()[i];if(data_->topology.count[r]>1){acceleration[r]+=-d->rates[i]+product(u[i],data_->outflow[i]);groupMass[r]+=material().midpointMass().masses()[i];groupDissipation[r]+=d->dissipation[i];}}
        for(int r=0;r<layout().velocityCount();++r)if(data_->topology.count[r]>1)acceleration[r]=finite(acceleration[r]/finite(groupMass[r],"nonfinite periodic class mass"),"nonfinite periodic acceleration");
        for(int i=0;i<layout().velocityCount();++i){const int r=representatives()[i];if(data_->topology.count[r]<=1)continue;const int a=component(layout(),i);
            const auto massLoad=product(material().midpointMass().masses()[i],acceleration[r]),outward=product(u[i],data_->outflow[i]);d->rates[i]=finite(outward-massLoad,"nonfinite periodic native convection");
            d->dissipation[i]=product(material().midpointMass().masses()[i]/groupMass[r],groupDissipation[r]);
            momentumScale[a]+=std::abs(outward)+std::abs(massLoad);powerScale[a]+=std::abs(product(u[i],outward))+std::abs(product(u[i],massLoad));}
    }
    for(int i=0;i<layout().velocityCount();++i){const int a=component(layout(),i);auto& report=d->report.components[a];
        const auto consistency=finite(impulse[i]/material().deltaTime(),"nonfinite consistency rate");
        d->loads[i]=finite(-finite(d->rates[i],"nonfinite compatible convection")+consistency,"nonfinite material momentum load");
        const int r=representatives()[i];if(data_->topology.paired && data_->topology.count[r]>1)
            d->loads[i]=finite(product(material().midpointMass().masses()[i],acceleration[r])+product(u[i],material().cells()[i].massChange/material().deltaTime()),"nonfinite periodic material load");
        const auto power=product(u[i],d->rates[i]),sourcePower=product(.5L*u[i],consistency);
        const auto massPower=product(.5L*u[i],product(u[i],material().cells()[i].massChange/material().deltaTime()));
        const auto kinetic=finite(product(u[i],d->loads[i])-massPower,"nonfinite material kinetic rate");
        report.momentumRate+=d->rates[i];report.convectivePower+=power;report.materialKineticPower+=sourcePower;report.kineticRate+=kinetic;
        momentumScale[a]+=std::abs(d->rates[i]);powerScale[a]+=std::abs(power)+std::abs(2*sourcePower)+std::abs(massPower)+std::abs(kinetic);
    }
    for(int a=0;a<2;++a){auto& r=d->report.components[a];
        r.momentumIdentity=identity(r.momentumRate-r.boundaryMomentumRate,momentumScale[a],"compatible momentum balance failed");
        r.kineticIdentity=identity(r.kineticRate-r.materialKineticPower+r.boundaryKineticPower+r.dissipationPower,powerScale[a],"nonconservative kinetic identity failed");
    }
    return CchKepConvectionResult(std::move(d));
}
CchKepConvectionResult CchKepConvection::applyPaired(const QuadtreeMesh& mesh,const CchNativeValues& values,CchMaterialSources source,Interface::HermiteGalsTimeLayer clock) const
{
    require(data_->topology.paired,"paired apply requires declared port topology");material().validateFor(mesh,source,clock);values.validateFor(layout(),mesh);
    std::vector<double> traces(material().geometry().segments().size());for(const auto& p:periodicPairs()){
        const auto& s=material().geometry().segments();const double trace=.5*values.normalVelocity()[owner(s[p.firstSegment])]+.5*values.normalVelocity()[owner(s[p.secondSegment])];traces[p.firstSegment]=traces[p.secondSegment]=trace;}
    return apply(mesh,values,CchConvectionBoundary::bind(mesh,material(),traces,source,clock),source,clock);
}
const CchDualMaterial& CchKepConvectionResult::material() const noexcept{return data_->material;}
const CchConvectionBoundary& CchKepConvectionResult::boundary() const noexcept{return data_->boundary;}
const std::vector<long double>& CchKepConvectionResult::centerMomentumRates() const noexcept{return data_->center;}
const std::vector<long double>& CchKepConvectionResult::momentumRates() const noexcept{return data_->rates;}
const std::vector<long double>& CchKepConvectionResult::momentumLoads() const noexcept{return data_->loads;}
const std::vector<long double>& CchKepConvectionResult::dissipationRates() const noexcept{return data_->dissipation;}
const std::vector<CchConvectionDissipationSample>& CchKepConvectionResult::dissipationSamples() const noexcept{return data_->dissipationSamples;}
const CchKepConvectionReport& CchKepConvectionResult::report() const noexcept{return data_->report;}
}
