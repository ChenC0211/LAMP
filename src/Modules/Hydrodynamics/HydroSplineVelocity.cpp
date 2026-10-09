#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Field/VEdgeField.h"
#include "Core/Field/HEdgeField.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace {
void require(bool ok,const char* reason){if(!ok)throw std::invalid_argument(reason);}
using Coefficients=std::array<long double,4>;
// Each polynomial uses the local unit coordinate, including for non-dyadic
// native coordinates. Natural endpoint second derivatives are zero.
std::vector<Coefficients> natural(const std::vector<double>& x,const std::vector<long double>& values) {
    const std::size_t n=x.size();require(n>=3 && values.size()==n,"spline requires at least three actual nodes");
    std::vector<long double> diagonal(n,1),rhs(n,0),upper(n,0),second(n,0);
    for(std::size_t i=1;i+1<n;++i) {
        const long double a=static_cast<long double>(x[i])-x[i-1],b=static_cast<long double>(x[i+1])-x[i];
        require(a>0 && b>0,"unordered native spline coordinates");
        diagonal[i]=2*(a+b);upper[i]=b;
        rhs[i]=6*((values[i+1]-values[i])/b-(values[i]-values[i-1])/a);
        const long double factor=a/diagonal[i-1];diagonal[i]-=factor*upper[i-1];rhs[i]-=factor*rhs[i-1];
    }
    for(std::size_t i=n-1;i-->0;)second[i]=(rhs[i]-upper[i]*second[i+1])/diagonal[i];
    std::vector<Coefficients> out(n-1);
    for(std::size_t i=0;i+1<n;++i) {
        const long double h=static_cast<long double>(x[i+1])-x[i];
        out[i]={values[i],values[i+1]-values[i]-h*h*(2*second[i]+second[i+1])/6,
            h*h*second[i]/2,h*h*(second[i+1]-second[i])/6};
    }
    return out;
}
struct TensorPrimitive {
    std::vector<std::array<long double,16>> cells;
    TensorPrimitive(const std::vector<double>& x,const std::vector<double>& y,const std::vector<long double>& values) {
        const std::size_t nx=x.size()-1,ny=y.size()-1;
        std::vector<std::vector<Coefficients>> rows(y.size());
        for(std::size_t j=0;j<y.size();++j) {
            std::vector<long double> row(x.size());
            for(std::size_t i=0;i<x.size();++i)row[i]=values[j*x.size()+i];
            rows[j]=natural(x,row);
        }
        cells.resize(nx*ny);
        for(std::size_t i=0;i<nx;++i)for(int a=0;a<4;++a) {
            std::vector<long double> column(y.size());for(std::size_t j=0;j<y.size();++j)column[j]=rows[j][i][a];
            const auto splines=natural(y,column);
            for(std::size_t j=0;j<ny;++j)for(int b=0;b<4;++b)cells[j*nx+i][4*a+b]=splines[j][b];
        }
    }
    long double evaluate(std::size_t cell,long double x,long double y,int dx,int dy) const {
        const auto& c=cells[cell];long double result=0;
        for(int a=3;a>=dx;--a) {
            long double row=0;for(int b=3;b>=dy;--b) {
                int factor=1;for(int k=0;k<dx;++k)factor*=a-k;for(int k=0;k<dy;++k)factor*=b-k;
                row=row*y+c[4*a+b]*factor;
            }
            result=result*x+row;
        }
        return result;
    }
};
struct SplineLayout {
    std::vector<double> x,y;
    std::vector<int> nativeCells;
    std::vector<int> verticalRows,horizontalRows,cellSlots;
    std::vector<double> faceAreas;
    SplineLayout(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) {
        domain.validateMeshIdentity(mesh);
        require(!domain.axisymmetric() && domain.uniformLeaves() && domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
            "C1 spline velocity currently requires a complete uniform XY domain");
        const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);const int nx=grid.nx(),ny=grid.ny();
        require(nx>=2 && ny>=2,"spline velocity requires at least two cells in each direction");
        for(int i=0;i<=nx;++i)x.push_back(mesh.verts(grid.vertices()[i]).coord().x);
        for(int j=0;j<=ny;++j)y.push_back(mesh.verts(grid.vertices()[j*(nx+1)]).coord().y);
        for(int j=0;j<=ny;++j)for(int i=0;i<=nx;++i) {
            const auto p=mesh.verts(grid.vertices()[j*(nx+1)+i]).coord();
            require(p.x==x[i] && p.y==y[j],"spline velocity requires a native tensor-product grid");
        }
        for(int label:grid.cells())nativeCells.push_back(mesh.cells(label).index());
        cellSlots.resize(nativeCells.size());for(std::size_t slot=0;slot<nativeCells.size();++slot)cellSlots[nativeCells[slot]]=static_cast<int>(slot);
        for(int label:grid.vertical())verticalRows.push_back(mesh.vEdges(label).index());
        for(int label:grid.horizontal())horizontalRows.push_back(mesh.vEdgeIndexNum()+mesh.hEdges(label).index());
        faceAreas.resize(domain.nativeFaceCount());for(const auto face:domain.faces(mesh))faceAreas[face.nativeRow]=face.area;
    }
};
// Coefficients of the original linear natural() operator applied to unit
// nodal data. Transposition reuses that exact construction, including endpoint
// conditions and metric; there is no sampled/least-squares velocity basis.
using NaturalBasis=std::vector<std::vector<Coefficients>>;
NaturalBasis naturalBasis(const std::vector<double>& nodes) {
    NaturalBasis basis;std::vector<long double> unit(nodes.size(),0);
    for(std::size_t node=0;node<nodes.size();++node){unit[node]=1;basis.push_back(natural(nodes,unit));unit[node]=0;}return basis;
}
std::vector<long double> tensorTranspose(const NaturalBasis& bx,const NaturalBasis& by,
    const std::vector<std::array<long double,16>>& moments) {
    const std::size_t nx=bx.size()-1,ny=by.size()-1;
    std::vector<long double> rows((ny+1)*nx*4,0),result((nx+1)*(ny+1),0);
    for(std::size_t node=0;node<=ny;++node)for(std::size_t i=0;i<nx;++i)for(int a=0;a<4;++a)
        for(std::size_t j=0;j<ny;++j)for(int b=0;b<4;++b)rows[(node*nx+i)*4+a]+=by[node][j][b]*moments[j*nx+i][4*a+b];
    for(std::size_t j=0;j<=ny;++j)for(std::size_t node=0;node<=nx;++node)
        for(std::size_t i=0;i<nx;++i)for(int a=0;a<4;++a)result[j*(nx+1)+node]+=bx[node][i][a]*rows[(j*nx+i)*4+a];
    return result;
}
}

struct HydroSplineVelocity::Impl : SplineLayout {
    std::unique_ptr<TensorPrimitive> px,py;
    std::vector<HydroFittedFaceFlux> faces;
    Impl(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const VEdgeField<double>& ux,const HEdgeField<double>& uy)
        :SplineLayout(mesh,domain) {
        require(ux.value && uy.value && ux.size()==mesh.vEdgeIndexNum() && uy.size()==mesh.hEdgeIndexNum(),"spline native field size mismatch");
        const int nx=static_cast<int>(x.size())-1,ny=static_cast<int>(y.size())-1;
        std::vector<long double> fx((nx+1)*(ny+1),0),fy(fx.size(),0);
        // Interpolation of the primitive and differentiation commute:
        // ux = d_y Sx Sy Fx, uy = d_x Sx Sy Fy. Each primitive increment is
        // the original FV flux. A zero discrete divergence gives d_xy(Fx+Fy)=0.
        for(int i=0;i<=nx;++i)for(int j=0;j<ny;++j) {
            const int row=verticalRows[j*(nx+1)+i];const double u=ux.value[row];
            require(std::isfinite(u) && (i!=0 && i!=nx || u==0),"spline velocity requires finite closed-box normal traces");
            const long double q=static_cast<long double>(faceAreas[row])*u;
            fx[(j+1)*(nx+1)+i]=fx[j*(nx+1)+i]+q;
        }
        for(int j=0;j<=ny;++j)for(int i=0;i<nx;++i) {
            const int row=horizontalRows[j*nx+i];const double u=uy.value[row-mesh.vEdgeIndexNum()];
            require(std::isfinite(u) && (j!=0 && j!=ny || u==0),"spline velocity requires finite closed-box normal traces");
            const long double q=static_cast<long double>(faceAreas[row])*u;
            fy[j*(nx+1)+i+1]=fy[j*(nx+1)+i]+q;
        }
        px=std::make_unique<TensorPrimitive>(x,y,fx);py=std::make_unique<TensorPrimitive>(x,y,fy);
        faces.reserve(domain.nativeFaceCount());
        for(const auto face:domain.faces(mesh)) {
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const HydroFittedFaceFlux::Chart chart{mesh.verts(edge.edgeToVert[0]).coord(),edge.midpointCoord(),mesh.verts(edge.edgeToVert[1]).coord()};
            const long double q=static_cast<long double>(face.area)*(face.horizontal?uy.value[edge.index()]:ux.value[edge.index()]);
            require(face.nativeRow==static_cast<int>(faces.size()),"spline flux lost native row ordering");
            // The spline normal trace is quadratic on this face. Endpoint
            // values and the native moment uniquely determine that polynomial.
            faces.push_back(HydroFittedFaceFlux::fromEndpointVelocities(chart,face.horizontal,q,
                firstJet(chart.begin).value,firstJet(chart.end).value,false));
        }
    }
    HydroFittedVelocitySampler::FirstJet firstJet(Vector2D p) const {
        require(std::isfinite(p.x) && std::isfinite(p.y) && p.x>=x.front() && p.x<=x.back() && p.y>=y.front() && p.y<=y.back(),
            "spline velocity query left the native closed domain");
        const auto interval=[](const auto& nodes,double q){return static_cast<std::size_t>(std::min<std::ptrdiff_t>(nodes.size()-2,
            std::upper_bound(nodes.begin(),nodes.end(),q)-nodes.begin()-1));};
        const std::size_t i=interval(x,p.x),j=interval(y,p.y),cell=j*(x.size()-1)+i;
        const long double hx=static_cast<long double>(x[i+1])-x[i],hy=static_cast<long double>(y[j+1])-y[j];
        const long double a=(static_cast<long double>(p.x)-x[i])/hx,b=(static_cast<long double>(p.y)-y[j])/hy;
        HydroFittedVelocitySampler::FirstJet out;out.nativeCell=nativeCells[cell];
        out.value={static_cast<double>(px->evaluate(cell,a,b,0,1)/hy),static_cast<double>(py->evaluate(cell,a,b,1,0)/hx)};
        out.dx={static_cast<double>(px->evaluate(cell,a,b,1,1)/(hx*hy)),static_cast<double>(py->evaluate(cell,a,b,2,0)/(hx*hx))};
        out.dy={static_cast<double>(px->evaluate(cell,a,b,0,2)/(hy*hy)),static_cast<double>(py->evaluate(cell,a,b,1,1)/(hx*hy))};
        // Essential traces of this closed spline space are exactly zero.
        // Evaluating a terminal power polynomial at one otherwise introduces
        // cancellation roundoff and fictitious bidirectional wall transport.
        // Use the prescribed trace for BOTH velocity and its tangential jet;
        // the normal derivative is still evaluated from the same polynomial.
        if(p.x==x.front() || p.x==x.back()){out.value.x=0;out.dy.x=0;}
        if(p.y==y.front() || p.y==y.back()){out.value.y=0;out.dx.y=0;}
        return out;
    }
};
struct HydroSplineVelocityAdjoint::Impl : SplineLayout {
    NaturalBasis bx,by;
    std::vector<std::array<long double,16>> mx,my;
    Impl(const QuadtreeMesh& mesh,const NativeLeafDomain& domain):SplineLayout(mesh,domain),bx(naturalBasis(x)),by(naturalBasis(y)),
        mx(nativeCells.size()),my(nativeCells.size()){}
};
HydroSplineVelocityAdjoint::HydroSplineVelocityAdjoint(const QuadtreeMesh& mesh,const NativeLeafDomain& domain)
    :impl_(std::make_unique<Impl>(mesh,domain)){}
HydroSplineVelocityAdjoint::~HydroSplineVelocityAdjoint()=default;
HydroSplineVelocityAdjoint::Moment HydroSplineVelocityAdjoint::moment(int nativeCell,Vector2D local,Vector2D force) const {
    require(std::isfinite(local.x) && std::isfinite(local.y) && local.x>=0 && local.x<=1 && local.y>=0 && local.y<=1
        && std::isfinite(force.x) && std::isfinite(force.y),"invalid spline load quadrature");
    const auto& p=*impl_;const int slot=p.cellSlots.at(nativeCell),nx=static_cast<int>(p.x.size())-1,i=slot%nx,j=slot/nx;
    const long double hx=static_cast<long double>(p.x[i+1])-p.x[i],hy=static_cast<long double>(p.y[j+1])-p.y[j];
    std::array<long double,4> ax{{1,local.x,static_cast<long double>(local.x)*local.x,static_cast<long double>(local.x)*local.x*local.x}},
        ay{{1,local.y,static_cast<long double>(local.y)*local.y,static_cast<long double>(local.y)*local.y*local.y}};
    Moment result{};
    for(int a=0;a<4;++a)for(int b=0;b<4;++b){if(b)result[4*a+b]=force.x*ax[a]*b*ay[b-1]/hy;
        if(a)result[16+4*a+b]=force.y*a*ax[a-1]*ay[b]/hx;}return result;
}
void HydroSplineVelocityAdjoint::accumulate(int nativeCell,const Moment& values){
    const auto slot=impl_->cellSlots.at(nativeCell);
    for(int k=0;k<16;++k){require(std::isfinite(values[k]) && std::isfinite(values[16+k]),"nonfinite spline load moment");
        impl_->mx[slot][k]+=values[k];impl_->my[slot][k]+=values[16+k];}
}
std::vector<double> HydroSplineVelocityAdjoint::integratedForce() const {
    const auto& p=*impl_;const int nx=static_cast<int>(p.x.size())-1,ny=static_cast<int>(p.y.size())-1;
    const auto fx=tensorTranspose(p.bx,p.by,p.mx),fy=tensorTranspose(p.bx,p.by,p.my);std::vector<double> result(p.faceAreas.size(),0);
    for(int i=0;i<=nx;++i){long double sum=0;for(int j=ny;j-->0;){sum+=fx[(j+1)*(nx+1)+i];
        const int row=p.verticalRows[j*(nx+1)+i];result[row]=static_cast<double>(p.faceAreas[row]*sum);}}
    for(int j=0;j<=ny;++j){long double sum=0;for(int i=nx;i-->0;){sum+=fy[j*(nx+1)+i+1];
        const int row=p.horizontalRows[j*nx+i];result[row]=static_cast<double>(p.faceAreas[row]*sum);}}
    for(double value:result)require(std::isfinite(value),"spline native load overflow");return result;
}
HydroSplineVelocity::HydroSplineVelocity(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const VEdgeField<double>& ux,const HEdgeField<double>& uy):impl_(std::make_unique<Impl>(mesh,domain,ux,uy)){}
HydroSplineVelocity::~HydroSplineVelocity()=default;
HydroSplineVelocity::HydroSplineVelocity(HydroSplineVelocity&&) noexcept=default;
HydroSplineVelocity& HydroSplineVelocity::operator=(HydroSplineVelocity&&) noexcept=default;
HydroFittedVelocitySampler::FirstJet HydroSplineVelocity::firstJet(Vector2D p) const{return impl_->firstJet(p);}
const HydroFittedFaceFlux& HydroSplineVelocity::faceFlux(int row) const{return impl_->faces.at(row);}
Vector2D HydroSplineVelocity::lower() const{return {impl_->x.front(),impl_->y.front()};}
Vector2D HydroSplineVelocity::upper() const{return {impl_->x.back(),impl_->y.back()};}
bool HydroSplineVelocity::sameLayout(const HydroSplineVelocity& other) const {
    if(impl_->x!=other.impl_->x || impl_->y!=other.impl_->y || impl_->nativeCells!=other.impl_->nativeCells || impl_->faces.size()!=other.impl_->faces.size())return false;
    for(std::size_t i=0;i<impl_->faces.size();++i){const auto& a=impl_->faces[i].chart();const auto& b=other.impl_->faces[i].chart();
        if(a.begin.x!=b.begin.x || a.begin.y!=b.begin.y || a.midpoint.x!=b.midpoint.x || a.midpoint.y!=b.midpoint.y || a.end.x!=b.end.x || a.end.y!=b.end.y)return false;}
    return true;
}
