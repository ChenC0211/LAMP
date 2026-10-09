#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineNative.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineBasis.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include <Eigen/QR>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <chrono>
#include <cstdlib>
#include <iostream>

namespace Quantum::Hydrodynamics {
namespace {
using namespace MomentSplineBasis;
using Row=std::map<int,double>;
constexpr int width=9,terms=width*width;
using Jets=Eigen::Matrix<double,6,terms>;
void require(bool ok,const char* reason){if(!ok)throw std::invalid_argument(reason);}
void add(Row& a,const Row& b,double scale){for(const auto& v:b)a[v.first]+=scale*v.second;}
Eigen::Matrix<double,6,6> hermite(){
    Eigen::Matrix<double,6,6> a;
    a<<1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,.5,0,0,0,
       -10,-6,-1.5,10,-4,.5, 15,8,1.5,-15,7,-1, -6,-3,-.5,6,-3,.5;return a;
}
Eigen::Matrix<double,width,3> edgeInterior(){
    Eigen::Matrix<double,width,3> b=Eigen::Matrix<double,width,3>::Zero();
    b(3,0)=1;b(4,0)=-3;b(5,0)=3;b(6,0)=-1;
    for(int k=1;k<3;k++)for(int a=0;a<width;a++){b(a,k)-=.5*b(a,k-1);if(a+1<width)b(a+1,k)+=b(a,k-1);}return b;
}
Eigen::Matrix<double,width,width> coordinates(){
    Eigen::Matrix<double,width,width> a=Eigen::Matrix<double,width,width>::Zero();
    a(0,0)=1;a(1,1)=1;a(2,2)=2;
    for(int k=0;k<width;k++){a(3,k)=1;a(4,k)=k;a(5,k)=k*(k-1);}
    a(8,8)=-1;a(7,7)=-1;a(7,8)=-4;a(6,6)=-1;a(6,7)=-3.5;a(6,8)=-7.75;return a;
}
Eigen::Matrix<long double,width,width> restriction(double offset,double scale){
    using M=Eigen::Matrix<long double,width,width>;
    M powers=M::Zero();powers.topLeftCorner<6,6>()=hermite().cast<long double>();powers.rightCols<3>()=edgeInterior().cast<long double>();
    M change=M::Zero();
    for(int k=0;k<width;k++){long double choose=1;for(int j=0;j<=k;j++){
        if(j)choose=choose*(k-j+1)/j;
        change(j,k)=choose*std::pow(static_cast<long double>(offset),k-j)*std::pow(static_cast<long double>(scale),j);
    }}
    // Known endpoint/bubble coordinate maps, not a fitted Vandermonde or
    // numerical inverse. Dyadic child intervals have exact binary entries.
    return coordinates().cast<long double>()*change*powers;
}
Eigen::Matrix<double,3,width> oneDimensionalJets(double t){
    using Jet=std::array<double,3>;Eigen::Matrix<double,3,width> out;
    const auto multiply=[](Jet a,Jet b)->Jet{return {a[0]*b[0],a[1]*b[0]+a[0]*b[1],a[2]*b[0]+2*a[1]*b[1]+a[0]*b[2]};};
    // Only powers 0..3 occur in this factored endpoint/interior basis.
    const auto power=[](double x,int n)->Jet{
        switch(n){case 0:return {1,0,0};case 1:return {x,1,0};case 2:return {x*x,2*x,2};default:return {x*x*x,3*x*x,6*x};}
    };
    for(int end=0;end<2;end++)for(int d=0;d<3;d++){
        const double x=end?1-t:t;auto v=power(1-x,3);v[1]=-v[1];
        const Jet q=d==0?Jet{1+3*x+6*x*x,3+12*x,12}:d==1?Jet{1+3*x,3,0}:Jet{.5,0,0};
        const auto basis=multiply(multiply(power(x,d),v),q);
        for(int k=0;k<3;k++)out(k,3*end+d)=basis[k]*((end&&(d+k)%2)?-1:1);
    }
    auto right=power(1-t,3);right[1]=-right[1];const auto base=multiply(power(t,3),right);
    for(int k=0;k<3;k++){const auto value=multiply(base,power(t-.5,k));for(int d=0;d<3;d++)out(d,6+k)=value[d];}return out;
}
// On an axis cell the physical streamfunction has an exact r^2 factor.
// Divide its endpoint-jet basis analytically, before evaluation, so neither
// 1/r nor cancellation of almost equal values enters the regular limit.
Eigen::Matrix<double,3,width> axisQuotientJets(double t){
    using Jet=std::array<double,3>;
    const auto multiply=[](Jet a,Jet b)->Jet{return {a[0]*b[0],a[1]*b[0]+a[0]*b[1],a[2]*b[0]+2*a[1]*b[1]+a[0]*b[2]};};
    const Jet x{t,1,0},q{1-t,-1,0},one{1,0,0};const auto q2=multiply(q,q),q3=multiply(q2,q);
    std::array<Jet,width> v{};v[2]={.5*q3[0],.5*q3[1],.5*q3[2]};
    v[3]=multiply(x,{1+3*q[0]+6*q2[0],3*q[1]+6*q2[1],6*q2[2]});
    v[4]=multiply(multiply(x,q),{-1-3*q[0],-3*q[1],0});
    v[5]=multiply(x,{.5*q2[0],.5*q2[1],.5*q2[2]});
    Jet factor=one;for(int k=0;k<3;k++){v[6+k]=multiply(multiply(x,q3),factor);factor=multiply(factor,{t-.5,1,0});}
    Eigen::Matrix<double,3,width> out;for(int k=0;k<width;k++)for(int d=0;d<3;d++)out(d,k)=v[k][d];return out;
}
using AxisJets=Eigen::Matrix<double,3,width>;
AxisJets exteriorJets(double t){
    const int end=t<0?0:1;const double d=t-end;AxisJets out=AxisJets::Zero();
    out(0,3*end)=1;out(0,3*end+1)=d;out(0,3*end+2)=.5*d*d;
    out(1,3*end+1)=1;out(1,3*end+2)=d;out(2,3*end+2)=1;return out;
}
std::array<AxisJets,2> tensorJets(Vector2D lower,Vector2D h,Vector2D p,bool axisQuotient,std::array<bool,2> exterior={false,false}){
    const double x=(p.x-lower.x)/h.x,y=(p.y-lower.y)/h.y;
    auto px=exterior[0]?exteriorJets(x):(axisQuotient?axisQuotientJets(x):oneDimensionalJets(x)),py=exterior[1]?exteriorJets(y):oneDimensionalJets(y);
    if(axisQuotient)px/=h.x*h.x;
    px.row(1)/=h.x;px.row(2)/=h.x*h.x;py.row(1)/=h.y;py.row(2)/=h.y*h.y;
    return {px,py};
}
Jets scalarJets(Vector2D lower,Vector2D h,Vector2D p,bool axisQuotient=false){
    const auto axes=tensorJets(lower,h,p,axisQuotient);
    Jets s;const int dx[6]={0,1,0,2,1,0},dy[6]={0,0,1,0,1,2};
    for(int k=0;k<6;k++)for(int a=0;a<width;a++)for(int b=0;b<width;b++)s(k,width*a+b)=axes[0](dx[k],a)*axes[1](dy[k],b);
    return s;
}
// Rows are physical velocity and its two derivatives. All evaluation,
// moments, kinetic energy and virtual work consume this one polynomial.
template<int Columns>
Eigen::Matrix<double,6,Columns> velocityJets(const Eigen::Matrix<double,6,Columns>& s,double r,bool rz,bool physicalPsi){
    Eigen::Matrix<double,6,Columns> v;
    if(physicalPsi){v.row(0)=s.row(2)/r;v.row(1)=-s.row(1)/r;
        v.row(2)=s.row(4)/r-s.row(2)/(r*r);v.row(3)=-s.row(3)/r+s.row(1)/(r*r);
        v.row(4)=s.row(5)/r;v.row(5)=-s.row(4)/r;}
    else if(rz){v.row(0)=r*s.row(2);v.row(1)=-2*s.row(0)-r*s.row(1);
        v.row(2)=s.row(2)+r*s.row(4);v.row(3)=-3*s.row(1)-r*s.row(3);
        v.row(4)=r*s.row(5);v.row(5)=-2*s.row(2)-r*s.row(4);}
    else {v.row(0)=s.row(2);v.row(1)=-s.row(1);v.row(2)=s.row(4);v.row(3)=-s.row(3);v.row(4)=s.row(5);v.row(5)=-s.row(4);}
    return v;
}
Jets physicalJets(Vector2D lower,Vector2D h,Vector2D p,bool rz,bool generalTrace){
    const auto s=scalarJets(lower,h,p,rz&&generalTrace&&lower.x==0);
    return velocityJets(s,p.x,rz,rz&&generalTrace&&lower.x!=0);
}
Eigen::Matrix<double,6,1> physicalValue(Vector2D lower,Vector2D h,Vector2D p,bool rz,bool generalTrace,const double* coefficients,std::array<bool,2> exterior={false,false}){
    const bool quotient=rz&&generalTrace&&lower.x==0&&!exterior[0];
    const auto axes=tensorJets(lower,h,p,quotient,exterior);
    // Contract the very same tensor-product polynomial before constructing
    // its physical jet. Avoid forming six 81-entry basis rows at every
    // characteristic query. Keep the factored Hermite basis and analytic
    // axis quotient; raw power coefficients would lose endpoint accuracy.
    const Eigen::Map<const Eigen::Matrix<double,width,width,Eigen::RowMajor>> patch(coefficients);
    const AxisJets contracted=axes[0]*patch;
    Eigen::Matrix<double,6,1> s;const int dx[6]={0,1,0,2,1,0},dy[6]={0,0,1,0,1,2};
    for(int k=0;k<6;k++)s[k]=contracted.row(dx[k]).dot(axes[1].row(dy[k]));
    return velocityJets(s,p.x,rz,rz&&generalTrace&&!quotient);
}
double measure(double r,bool rz){return rz?CartesianMeasure::twoPi*r:1.;}
void wallJet(Eigen::Matrix<double,6,1>& v,Vector2D p,Vector2D lo,Vector2D hi,bool rz,bool generalTrace){
    if(rz&&p.x==0){v[0]=v[4]=v[3]=0;}
    else if(!generalTrace&&(p.x==lo.x||p.x==hi.x)){v[0]=v[1]=v[4]=v[5]=v[2]=0;}
    if(!generalTrace&&(p.y==lo.y||p.y==hi.y)){v[0]=v[1]=v[2]=v[3]=v[5]=0;}
}
}

HydroCompatibleMomentSpline::NativeLayout::NativeLayout(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentBoundary* boundary,bool prepareMomentLift):locator(mesh),identity(domain.meshIdentity()),rz(domain.axisymmetric()),
    generalTrace(boundary&&boundary->usesGeneralTraceSpace()),permeable(boundary&&boundary->open()){
    const auto started=std::chrono::steady_clock::now();
    const auto H=hermite();const auto B=edgeInterior();
    const auto timing=[&](const char* stage){if(std::getenv("QUANTUM_MOMENT_SETUP_PROFILE"))std::cout<<"MOMENT_C1_SETUP stage="<<stage<<" cells="<<mesh.cellIndexNum()<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<std::endl;};
    const double infinity=std::numeric_limits<double>::infinity();lower={infinity,infinity};upper={-infinity,-infinity};
    slotByCell.assign(mesh.cellIndexNum(),-1);
    struct Node {Vector2D p;double scale;int coarse=-1;bool horizontal=false;std::array<Row,9> jets;int state=0;};
    std::vector<Node> nodes;std::vector<std::array<int,4>> corners;
    std::map<int,int> vertex;
    for(int index:domain.cellIndexByRow()){
        const auto& source=mesh.cells(mesh.cellIndexToLabel[index]);
        const auto lo=mesh.verts(source.cellToVert[0]).coord(),hi=mesh.verts(source.cellToVert[2]).coord(),h=hi-lo;
        require(h.x>0&&h.y>0,"invalid native C1 cell");
        slotByCell[index]=int(cells.size());cells.push_back({index,lo,h});std::array<int,4> ids;
        lower.x=std::min(lower.x,lo.x);lower.y=std::min(lower.y,lo.y);upper.x=std::max(upper.x,hi.x);upper.y=std::max(upper.y,hi.y);
        for(int j=0;j<2;j++)for(int i=0;i<2;i++){
            // Topological identity belongs to the mesh vertex. Rebuilding
            // corners from centre +/- half-size creates distinct rounded
            // copies on non-binary spacings and false hanging-node cycles.
            const int label=source.cellToVert[j?(i?2:3):i];const auto p=mesh.verts(label).coord();
            const auto inserted=vertex.emplace(label,int(nodes.size()));
            if(inserted.second)nodes.push_back({p,std::min(h.x,h.y)});else nodes[inserted.first->second].scale=std::min(nodes[inserted.first->second].scale,std::min(h.x,h.y));
            ids[2*j+i]=inserted.first->second;
        }corners.push_back(ids);
    }
    require(cells.size()==std::size_t(mesh.cellIndexNum())&&(!rz||lower.x>=0),"native C1 requires complete nonnegative Cartesian/RZ leaves");
    require(!generalTrace||!rz||lower.x==0,"general boundary traces in RZ currently require an axis-containing domain");
    // Complete C2 Q8 patches: nodal jets, three edge-interior moments for
    // each normal jet, and nine cell-interior directions. Edge interiors
    // are part of the same space, not a post-hoc correction of a Q5 trace.
    struct Edge {int ends[2];double fixed,start,length;bool horizontal;int parent=-1,state=0;std::array<std::array<Row,width>,3> coefficients;};
    std::vector<Edge> edges;std::vector<std::array<int,4>> cellEdges(cells.size());std::map<std::tuple<bool,double,double,double>,int> edgeIndex;
    for(int slot=0;slot<int(cells.size());slot++)for(int direction=0;direction<2;direction++)for(int side=0;side<2;side++){
        const bool horizontal=direction==1;
        const int a=corners[slot][horizontal?2*side:side],b=corners[slot][horizontal?2*side+1:side+2];
        const double fixed=horizontal?nodes[a].p.y:nodes[a].p.x,start=horizontal?nodes[a].p.x:nodes[a].p.y,
            end=horizontal?nodes[b].p.x:nodes[b].p.y,length=end-start;
        const auto found=edgeIndex.emplace(std::make_tuple(horizontal,fixed,start,end),int(edges.size()));
        if(found.second)edges.push_back({{a,b},fixed,start,length,horizontal});cellEdges[slot][2*direction+side]=found.first->second;
    }
    for(auto& edge:edges){const auto a=nodes[edge.ends[0]].p,b=nodes[edge.ends[1]].p;const auto incident=locator.incident((a+b)*.5);double longest=edge.length;
        for(int k=0;k<incident.count;k++)for(int candidate:cellEdges[slotByCell.at(incident.indices[k])]){
            const auto& parent=edges[candidate];const auto end=nodes[parent.ends[1]].p;
            if(parent.horizontal==edge.horizontal&&parent.fixed==edge.fixed&&parent.start<=edge.start
                &&(edge.horizontal?end.x>=b.x:end.y>=b.y)&&parent.length>longest){edge.parent=candidate;longest=parent.length;}
        }
    }
    for(auto& node:nodes){const auto incident=locator.incident(node.p);double longest=0;
        for(int k=0;k<incident.count;k++)for(int candidate:cellEdges[slotByCell.at(incident.indices[k])]){
            const auto& edge=edges[candidate];const double normal=edge.horizontal?node.p.y:node.p.x,tangent=edge.horizontal?node.p.x:node.p.y;
            const auto end=nodes[edge.ends[1]].p;
            if(normal==edge.fixed&&tangent>edge.start&&tangent<(edge.horizontal?end.x:end.y)&&edge.length>longest){node.coarse=candidate;longest=edge.length;}
        }
    }
    int nodal=0;std::function<void(int)> resolveNode,resolveEdge;
    resolveNode=[&](int index){auto& node=nodes[index];if(node.state==2)return;require(node.state==0,"cyclic native C1 hanging node");node.state=1;
        if(node.coarse<0){for(int a=0;a<3;a++)for(int b=0;b<3;b++){
            const bool xwall=node.p.x==upper.x||((!rz||lower.x>0)&&node.p.x==lower.x),ywall=node.p.y==lower.y||node.p.y==upper.y;
            // Throughflow permits distinct streamfunction constants on
            // disconnected walls. Do not set them all to zero. Physical
            // no-slip is imposed below on the same velocity polynomial.
            // XY removes one scalar gauge. Open RZ uses physical Psi, zero
            // value/first radial jet at the axis, with higher even regularity
            // below. Closed RZ retains its existing chi coordinates.
            if((!generalTrace&&((xwall&&a<2)||(ywall&&b<2)))||(rz&&node.p.x==0&&(generalTrace?a<2:a==1))
                ||(generalTrace&&!rz&&node.p.x==lower.x&&node.p.y==lower.y&&a==0&&b==0))continue;
            node.jets[3*a+b][nodal++]=1;
        }}else {resolveEdge(node.coarse);const auto& edge=edges[node.coarse];const double t=((edge.horizontal?node.p.x:node.p.y)-edge.start)/edge.length;
            const auto jets=oneDimensionalJets(t);
            for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int k=0;k<width;k++)
                add(node.jets[3*a+b],edge.coefficients[edge.horizontal?b:a][k],jets(edge.horizontal?a:b,k)*std::pow(node.scale/edge.length,a+b));
        }node.state=2;
    };
    resolveEdge=[&](int index){auto& edge=edges[index];if(edge.state==2)return;require(edge.state==0,"cyclic native C1 hanging edge");edge.state=1;
        if(edge.parent>=0){resolveEdge(edge.parent);const auto& parent=edges[edge.parent];const double a=(edge.start-parent.start)/parent.length,b=edge.length/parent.length;
            const auto restricted=restriction(a,b);
            for(int n=0;n<3;n++)for(int j=0;j<width;j++)for(int k=0;k<width;k++)if(restricted(j,k)!=0)
                add(edge.coefficients[n][j],parent.coefficients[n][k],double(restricted(j,k)*std::pow(static_cast<long double>(b),n)));
        }else {
            for(int end=0;end<2;end++){resolveNode(edge.ends[end]);const auto& node=nodes[edge.ends[end]];
                for(int n=0;n<3;n++)for(int d=0;d<3;d++)
                    add(edge.coefficients[n][3*end+d],node.jets[edge.horizontal?3*d+n:3*n+d],std::pow(edge.length/node.scale,n+d));
            }
            const bool axis=rz&&!edge.horizontal&&edge.fixed==0;
            const bool wall=edge.horizontal?(edge.fixed==lower.y||edge.fixed==upper.y):(edge.fixed==upper.x||(!axis&&edge.fixed==lower.x));
            for(int n=0;n<3;n++)if(!(!generalTrace&&wall&&n<2)&&!(axis&&(generalTrace?n<2:n==1)))for(int b=0;b<3;b++){
                edge.coefficients[n][6+b][nodal++]=1;
            }
        }edge.state=2;
    };
    for(int k=0;k<int(nodes.size());k++)resolveNode(k);for(int k=0;k<int(edges.size());k++)resolveEdge(k);
    dofs=nodal+9*int(cells.size());std::vector<Eigen::Triplet<double>> poly,energy,moment,target;
    cellInteriorCoordinates.resize(cells.size());
    for(int slot=0;slot<int(cells.size());slot++)for(int k=0;k<9;k++)cellInteriorCoordinates[slot][k]=nodal+9*slot+k;
    const int nv=int(space.mass().rows());int equation=nv;const auto& gauss=Geometry::Parametric::gaussEight();const GaussTwelve radial;
    faces.resize(domain.nativeFaceCount());std::vector<std::vector<int>> cellFaces(cells.size());
    for(const auto f:domain.faces(mesh)){
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const int owner=f.negativeCell>=0?f.negativeCell:f.positiveCell,slot=slotByCell.at(owner);
        const auto kind=boundary?boundary->faces()[f.nativeRow].kind:HydroCompatibleMomentBoundary::Kind::Wall;
        faces[f.nativeRow]={{a,edge.midpointCoord(),b},slot,f.horizontal,f.outwardSign!=0,
            f.outwardSign&&kind==HydroCompatibleMomentBoundary::Kind::Wall,
            f.outwardSign&&HydroCompatibleMomentBoundary::impermeable(kind)};cellFaces[slot].push_back(f.nativeRow);
        if(f.outwardSign&&kind==HydroCompatibleMomentBoundary::Kind::PressureOutlet){
            const int side=(f.horizontal?2:0)+(f.outwardSign>0?1:0);
            outletSpans[side].push_back({f.horizontal?a.x:a.y,f.horizontal?b.x:b.y});}
        if((generalTrace||!f.outwardSign)&&space.faces()[f.nativeRow].first>=0)for(int k=0;k<2;k++)target.emplace_back(space.faces()[f.nativeRow].first+k,space.faces()[f.nativeRow].first+k,1);
    }
    for(auto& spans:outletSpans)std::sort(spans.begin(),spans.end());
    std::vector<int> independentAxisRows;std::set<std::pair<double,int>> axisJets;
    for(int slot=0;slot<int(cells.size());slot++){
        const auto& c=cells[slot];const auto h=c.size;const double rc=c.lower.x+.5*h.x;
        std::array<Row,terms> rows;
        for(int j=0;j<2;j++)for(int i=0;i<2;i++)for(int a=0;a<3;a++)for(int b=0;b<3;b++){
            const auto& node=nodes[corners[slot][2*j+i]];const double scale=std::pow(h.x/node.scale,a)*std::pow(h.y/node.scale,b);
            add(rows[width*(3*i+a)+3*j+b],node.jets[3*a+b],scale);
        }
        for(int direction=0;direction<2;direction++)for(int side=0;side<2;side++)for(int n=0;n<3;n++){
            const bool horizontal=direction==1;const auto& edge=edges[cellEdges[slot][2*direction+side]];
            for(int k=6;k<width;k++)add(rows[horizontal?width*k+3*side+n:width*(3*side+n)+k],
                edge.coefficients[n][k],std::pow((horizontal?h.y:h.x)/edge.length,n));
        }
        for(int a=0;a<3;a++)for(int b=0;b<3;b++)rows[width*(6+a)+6+b][nodal+9*slot+3*a+b]=1;
        std::set<int> columns;for(const auto& row:rows)for(const auto& v:row)if(v.second!=0)columns.insert(v.first);
        std::vector<int> ids(columns.begin(),columns.end());Eigen::MatrixXd lift=Eigen::MatrixXd::Zero(terms,ids.size());
        for(int k=0;k<terms;k++)for(const auto& v:rows[k])if(v.second!=0){const int col=int(std::lower_bound(ids.begin(),ids.end(),v.first)-ids.begin());lift(k,col)=v.second;}
        // Assemble the canonical tensor coordinates directly. Endpoint jets,
        // edge interiors and cell interiors have disjoint owned rows; their
        // identity must not depend on cancellation in a power-basis roundtrip.
        for(int k=0;k<terms;k++)for(int col=0;col<int(ids.size());col++)
            if(lift(k,col)!=0)poly.emplace_back(terms*slot+k,ids[col],lift(k,col));
        Eigen::MatrixXd local=Eigen::MatrixXd::Zero(ids.size(),ids.size());
        Eigen::Matrix<double,6,terms> cellMoments=Eigen::Matrix<double,6,terms>::Zero();
        const double volume=h.x*h.y*measure(rc,rz);
        for(int a=0;a<12;a++)for(int b=0;b<12;b++){
            const double x=.5+.5*radial.x[a],y=.5+.5*radial.x[b];const Vector2D p{c.lower.x+x*h.x,c.lower.y+y*h.y};
            const auto v=physicalJets(c.lower,h,p,rz,generalTrace);const double w=.25*h.x*h.y*radial.w[a]*radial.w[b]*measure(p.x,rz);
            const Eigen::MatrixXd velocities=v.topRows<2>()*lift;local+=w*velocities.transpose()*velocities;
            cellMoments.row(0)+=w/volume*v.row(0);cellMoments.row(1)+=w/volume*v.row(1);
            cellMoments.row(2)+=w/volume*(-(y-.5)*v.row(0)+(x-.5)*v.row(1));
        }
        if(rz)for(int panel=0;panel<4;panel++)for(int a=0;a<8;a++)for(int b=0;b<8;b++){
            const double x=(panel+.5+.5*gauss.x[a])/4,y=.5+.5*gauss.x[b],s=(c.lower.x+x*h.x)/rc,s0=c.lower.x/rc;
            const Vector2D p{c.lower.x+x*h.x,c.lower.y+y*h.y};const auto scalar=scalarJets(c.lower,h,p);
            Eigen::RowVectorXd psi;
            if(generalTrace){psi=scalar.row(0)/(rc*rc);psi[0]-=1/(rc*rc);}
            else {psi=s*s*scalar.row(0);psi[0]=s*s*(scalar(0,0)-1)+(h.x/rc)*x*(2*s0+(h.x/rc)*x);}
            const double w=.25*gauss.w[a]*gauss.w[b]*rc/(4*std::sqrt(h.x*h.y));
            cellMoments.row(3)+=w*(x-.5)/s*psi;cellMoments.row(4)+=w*(x-.5)*(x-.5)/(s*s)*psi;cellMoments.row(5)+=w*(y-.5)/(s*s)*psi;
        }
        for(int a=0;a<int(ids.size());a++)for(int b=0;b<int(ids.size());b++)if(local(a,b)!=0)energy.emplace_back(ids[a],ids[b],local(a,b));
        const auto append=[&](int row,const Eigen::RowVectorXd& functional){const Eigen::RowVectorXd values=functional*lift;
            for(int k=0;k<int(ids.size());k++)if(values[k]!=0)moment.emplace_back(row,ids[k],values[k]);};
        const int first=space.firstCellVelocity(c.index);
        for(int k=0;k<space.cellVelocityCount();k++){append(first+k,cellMoments.row(k));target.emplace_back(first+k,first+k,1);}
        for(int row:cellFaces[slot]){
            const auto& f=faces[row];const auto& info=space.faces()[row];if((!generalTrace&&f.boundary)||info.first<0)continue;
            const double length=f.horizontal?f.chart.end.x-f.chart.begin.x:f.chart.end.y-f.chart.begin.y;
            Eigen::Matrix<double,2,terms> moments=Eigen::Matrix<double,2,terms>::Zero();
            for(int k=0;k<8;k++){
                const double t=.5+.5*gauss.x[k];const auto p=f.chart.begin+(f.chart.end-f.chart.begin)*t;
                const auto v=physicalJets(c.lower,h,p,rz,generalTrace);
                const double w=.5*length*gauss.w[k]*measure(p.x,rz)/info.area;
                moments.row(0)+=w*v.row(f.horizontal?1:0);
                moments.row(1)+=w*((f.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale*v.row(f.horizontal?1:0);
            }append(info.first,moments.row(0));append(info.first+1,moments.row(1));
            if(generalTrace&&f.impermeable){
                // The trace (or its r-weighted numerator away from the RZ
                // axis) has degree <=8. Twelve distinct points constrain
                // the entire trace, not a sampled penalty. Physical Psi
                // permits a nonzero constant on horizontal/cylindrical
                // walls alike; wall constants are not pinned separately.
                for(int k=0;k<12;k++){
                    const auto p=f.chart.begin+(f.chart.end-f.chart.begin)*(.5+.5*radial.x[k]);
                    const auto v=physicalJets(c.lower,h,p,rz,generalTrace);
                    if(f.wall||!f.horizontal)append(equation++,v.row(0));
                    if(f.wall||f.horizontal)append(equation++,v.row(1));
                }
            }
        }
        // Even axis regularity: closed chi has zero first/third radial jets;
        // physical Psi has the exact r^2 factor plus zero third/fifth jets.
        if(rz&&c.lower.x==0)for(int d:{3,5})for(int b=0;b<width;b++){
            if(!generalTrace&&d==5)continue;
            Eigen::RowVectorXd value=Eigen::RowVectorXd::Zero(terms);
            for(int a=0;a<width;a++)value[a*width+b]=a<6?H(d,a):B(d,a-6);
            if(!generalTrace){
                // chi_rrr along the axis is one C2 piecewise Q8 trace:
                // three interior coefficients per interval and one shared
                // value/first/second jet at each end. No-slip already fixes
                // its value and first jet at the two axial domain walls.
                bool keep=b>=6;
                if(b<6){const double z=b<3?c.lower.y:c.lower.y+h.y;const int jet=b%3;
                    keep=!(jet<2&&(z==lower.y||z==upper.y))&&axisJets.emplace(z,jet).second;}
                if(keep)independentAxisRows.push_back(equation);
            }
            append(equation++,value);
        }
    }
    Matrix kinetic(dofs,dofs);kinetic.setFromTriplets(energy.begin(),energy.end());
    kinetic=Matrix(.5*(kinetic+Matrix(kinetic.transpose())));
    Eigen::VectorXd scaling(dofs);for(int k=0;k<dofs;k++){require(kinetic.coeff(k,k)>0,"empty native C1 energy direction");scaling[k]=1/std::sqrt(kinetic.coeff(k,k));}
    polynomial.resize(terms*cells.size(),dofs);polynomial.setFromTriplets(poly.begin(),poly.end());polynomial=polynomial*scaling.asDiagonal();
    // Every independent native Hermite coordinate has an owned endpoint,
    // edge-interior or cell-interior coefficient. Recover it from a single
    // exact row, rather than fitting or minimizing the retained field again.
    std::vector<int> rowCount(polynomial.rows(),0);
    for(int col=0;col<polynomial.outerSize();col++)for(Matrix::InnerIterator it(polynomial,col);it;++it)if(it.value()!=0)++rowCount[it.row()];
    coordinateRows.assign(dofs,-1);coordinateValues=Eigen::VectorXd::Zero(dofs);
    for(int col=0;col<polynomial.outerSize();col++)for(Matrix::InnerIterator it(polynomial,col);it;++it)
        if(rowCount[it.row()]==1&&std::abs(it.value())>std::abs(coordinateValues[col])){coordinateRows[col]=int(it.row());coordinateValues[col]=it.value();}
    allConstraints.resize(equation,dofs);allConstraints.setFromTriplets(moment.begin(),moment.end());allConstraints=allConstraints*scaling.asDiagonal();
    momentMap=allConstraints.topRows(nv);
    allData.resize(equation,nv);allData.setFromTriplets(target.begin(),target.end());
    Eigen::VectorXd norms=Eigen::VectorXd::Zero(equation);for(int col=0;col<allConstraints.outerSize();col++)for(Matrix::InnerIterator it(allConstraints,col);it;++it)norms[it.row()]+=it.value()*it.value();
    for(int row=0;row<equation;row++)norms[row]=norms[row]>0?1/std::sqrt(norms[row]):1;
    allConstraints=norms.asDiagonal()*allConstraints;allData=norms.asDiagonal()*allData;
    kinematicConstraints=allConstraints.bottomRows(equation-nv);
    timing("assembled");
    // Exact field restriction only needs the geometry, coordinate map and
    // kinematics. Do not factor an unused moment-to-field minimization.
    if(!prepareMomentLift)return;
    // Conservation makes some retained moments redundant on solenoidal
    // fields. Every original moment is checked against the completed field.
    int rank=0;std::vector<int> selected(equation,-1);
    if(generalTrace){
        // Physical wall traces add exact identities between multiple edges,
        // face moments and D. Rank-reveal the ORIGINAL normalized rows with
        // global norm pivoting. Fill-reducing sparse QR can overcount these
        // dependencies; compressing its already inaccurate R cannot repair
        // the rank. This geometric setup is reused for every velocity/work
        // query; the kinetic saddle and all physical checks stay sparse.
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> independent{Eigen::MatrixXd(allConstraints.transpose())};independent.setThreshold(1e-10);
        rank=int(independent.rank());for(int k=0;k<rank;k++)selected[independent.colsPermutation().indices()[k]]=k;
    }else{
        // Use the native exact divergence identities instead of discovering
        // them again by a global rank factorization of the lifted space.
        // P1 divergence fixes both cell means. Constant divergence fixes
        // the mean fluxes on a spanning tree of the connected leaf graph;
        // all other face moments and retained cell moments are independent.
        std::vector<int> parent(cells.size());std::iota(parent.begin(),parent.end(),0);
        const auto root=[&](int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
        int treeEdges=0;
        for(const auto face:domain.faces(mesh))if(!face.outwardSign){
            const int first=space.faces()[face.nativeRow].first;require(first>=0,"interior moment face has no velocity rows");
            const int a=root(slotByCell.at(face.negativeCell)),b=root(slotByCell.at(face.positiveCell));
            if(a!=b){parent[a]=b;++treeEdges;}else selected[first]=0;
            selected[first+1]=0;
        }
        require(treeEdges+1==int(cells.size()),"native C1 closed leaf graph is disconnected");
        for(const auto& cell:cells){const int first=space.firstCellVelocity(cell.index);
            require(space.divergence().coeff(3*cell.index+1,first)!=0&&space.divergence().coeff(3*cell.index+2,first+1)!=0,
                "native C1 cell means are not owned by P1 divergence");
            for(int k=2;k<space.cellVelocityCount();k++)selected[first+k]=0;
        }
        for(int row:independentAxisRows)selected[row]=0;
        for(int row=0;row<equation;row++)if(selected[row]==0)selected[row]=rank++;
    }
    timing("independent_constraints");
    std::vector<Eigen::Triplet<double>> entries,targets;kinetic=scaling.asDiagonal()*kinetic*scaling.asDiagonal();
    for(int col=0;col<kinetic.outerSize();col++)for(Matrix::InnerIterator it(kinetic,col);it;++it)entries.emplace_back(it.row(),it.col(),it.value());
    for(int col=0;col<allConstraints.outerSize();col++)for(Matrix::InnerIterator it(allConstraints,col);it;++it)if(selected[it.row()]>=0){
        entries.emplace_back(dofs+selected[it.row()],it.col(),it.value());entries.emplace_back(it.col(),dofs+selected[it.row()],it.value());}
    for(int col=0;col<allData.outerSize();col++)for(Matrix::InnerIterator it(allData,col);it;++it)if(selected[it.row()]>=0)targets.emplace_back(selected[it.row()],it.col(),it.value());
    data.resize(rank,nv);data.setFromTriplets(targets.begin(),targets.end());system.resize(dofs+rank,dofs+rank);system.setFromTriplets(entries.begin(),entries.end());
    solver.compute(system);
    timing("factorized");
    if(solver.info()!=Eigen::Success){
        if(dofs>4096){std::ostringstream message;message<<"native C1 kinetic factorization failed: dofs="<<dofs<<" selected="<<rank;throw std::invalid_argument(message.str());}
        Eigen::JacobiSVD<Eigen::MatrixXd> constraintRank{Eigen::MatrixXd(allConstraints)},polynomialRank{Eigen::MatrixXd(polynomial)};
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> energyRank{Eigen::MatrixXd(kinetic),Eigen::EigenvaluesOnly};
        const auto rankAt=[](const Eigen::VectorXd& s){int count=0;for(int k=0;k<s.size();k++)if(s[k]>1e-10*s[0])count++;return count;};
        const Eigen::MatrixXd selectedConstraints=Eigen::MatrixXd(system).bottomLeftCorner(rank,dofs);
        Eigen::JacobiSVD<Eigen::MatrixXd> selectedRank(selectedConstraints);
        const Eigen::MatrixXd physical=norms.cwiseInverse().asDiagonal()*allConstraints;
        const double commuting=(space.divergence()*physical.topRows(nv)).norm();
        std::ostringstream message;message<<"native C1 constrained kinetic space lost rank: dofs="<<dofs<<" selected="<<rank
            <<" constraint_rank="<<rankAt(constraintRank.singularValues())<<" polynomial_rank="<<rankAt(polynomialRank.singularValues())
            <<" kinetic_min="<<energyRank.eigenvalues()[0]<<" kinetic_max="<<energyRank.eigenvalues().tail(1)[0]
            <<" selected_rank="<<rankAt(selectedRank.singularValues())<<" selected_min="<<selectedRank.singularValues().tail(1)[0]
            <<" constraint_min="<<constraintRank.singularValues().tail(1)[0]<<" commuting="<<commuting;throw std::invalid_argument(message.str());
    }
}

Eigen::VectorXd HydroCompatibleMomentSpline::NativeLayout::solve(const Eigen::VectorXd& rhs) const{
    std::lock_guard<std::mutex> lock(mutex);Eigen::VectorXd solution=solver.solve(rhs);
    require(solver.info()==Eigen::Success&&solution.allFinite(),"native C1 kinetic solve failed");
    require((system*solution-rhs).norm()<=2e-9*rhs.norm(),"native C1 kinetic solve residual unresolved");return solution;
}
Eigen::MatrixXd HydroCompatibleMomentSpline::NativeLayout::coefficients(const Eigen::VectorXd& u) const{
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(system.rows());rhs.tail(data.rows())=data*u;const auto solution=solve(rhs);
    const Eigen::VectorXd requested=allData*u,defect=allConstraints*solution.head(dofs)-requested;
    if(defect.norm()>2e-9*requested.norm()){
        std::ostringstream text;text<<"native C1 cannot retain every requested moment: relative="<<defect.norm()/requested.norm()<<" dofs="<<dofs<<" rank="<<data.rows();throw std::invalid_argument(text.str());}
    const Eigen::VectorXd flat=polynomial*solution.head(dofs);return Eigen::Map<const Eigen::MatrixXd>(flat.data(),terms,cells.size());
}
Eigen::VectorXd HydroCompatibleMomentSpline::NativeLayout::fieldCoordinates(const Eigen::MatrixXd& coefficients) const{
    require(coefficients.rows()==terms&&coefficients.cols()==int(cells.size())&&coefficients.allFinite(),"invalid complete native polynomial field");
    require(std::all_of(coordinateRows.begin(),coordinateRows.end(),[](int row){return row>=0;}),"native Hermite coordinate has no owned coefficient");
    const Eigen::Map<const Eigen::VectorXd> flat(coefficients.data(),coefficients.size());Eigen::VectorXd result(dofs);
    for(int col=0;col<dofs;col++)result[col]=flat[coordinateRows[col]]/coordinateValues[col];
    const Eigen::VectorXd reconstructed=polynomial*result;
    require((reconstructed-flat).norm()<=2e-11*flat.norm(),"complete field is not in this native Hermite space");return result;
}
Eigen::MatrixXd HydroCompatibleMomentSpline::NativeLayout::restrictCoefficients(const NativeLayout& source,const Eigen::MatrixXd& value) const{
    require(rz==source.rz&&generalTrace==source.generalTrace&&permeable==source.permeable&&lower.x==source.lower.x&&lower.y==source.lower.y&&upper.x==source.upper.x&&upper.y==source.upper.y,
        "native refinement changed physical domain or scalar representation");
    require(value.rows()==terms&&value.cols()==int(source.cells.size())&&value.allFinite(),"invalid native refinement coefficients");
    Eigen::MatrixXd result(terms,cells.size());
    for(int slot=0;slot<int(cells.size());slot++){
        const auto& to=cells[slot];const int fromSlot=source.locate(to.lower+to.size*.5);const auto& from=source.cells[fromSlot];
        const double eps=64*std::numeric_limits<double>::epsilon()*std::max({from.size.x,from.size.y,std::abs(from.lower.x),std::abs(from.lower.y)});
        require(to.lower.x>=from.lower.x-eps&&to.lower.y>=from.lower.y-eps&&to.lower.x+to.size.x<=from.lower.x+from.size.x+eps&&to.lower.y+to.size.y<=from.lower.y+from.size.y+eps,
            "native polynomial restriction requires refinement, not coarsening across source patches");
        const auto rx=restriction((to.lower.x-from.lower.x)/from.size.x,to.size.x/from.size.x),ry=restriction((to.lower.y-from.lower.y)/from.size.y,to.size.y/from.size.y);
        Eigen::Matrix<long double,width,width> old;
        for(int a=0;a<width;a++)for(int b=0;b<width;b++)old(a,b)=value(width*a+b,fromSlot);
        const Eigen::Matrix<long double,width,width> restricted=rx*old*ry.transpose();
        for(int a=0;a<width;a++)for(int b=0;b<width;b++)result(width*a+b,slot)=double(restricted(a,b));
    }
    return result;
}
int HydroCompatibleMomentSpline::NativeLayout::locate(Vector2D p) const{
    require(std::isfinite(p.x)&&std::isfinite(p.y),"nonfinite native C1 query");const auto incident=locator.incident(p);
    if(incident.count==0){std::ostringstream message;message.precision(17);message<<"native C1 query outside physical domain: ("<<p.x<<", "<<p.y<<")";throw std::invalid_argument(message.str());}
    return slotByCell.at(incident.indices[0]);
}
void HydroCompatibleMomentSpline::NativeLayout::essentialJet(Eigen::Matrix<double,6,1>& v,Vector2D p) const{
    wallJet(v,p,lower,upper,rz,generalTrace);
    if(!generalTrace||(p.x!=lower.x&&p.x!=upper.x&&p.y!=lower.y&&p.y!=upper.y))return;
    // Exact evaluation of essential wall traces, including their port-rim
    // endpoints. This removes only solve roundoff from identities already
    // imposed on the entire polynomial. Apply the SAME map to virtual work.
    for(const auto& f:faces)if(f.impermeable){
        if(f.horizontal&&p.y==f.chart.begin.y&&p.x>=f.chart.begin.x&&p.x<=f.chart.end.x){
            v[1]=v[3]=0;if(f.wall)v[0]=v[2]=v[5]=0;}
        if(!f.horizontal&&(!rz||f.chart.begin.x!=0)&&p.x==f.chart.begin.x&&p.y>=f.chart.begin.y&&p.y<=f.chart.end.y){
            v[0]=v[4]=0;if(f.wall)v[1]=v[2]=v[5]=0;}
    }
}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentSpline::NativeLayout::firstJet(const Eigen::MatrixXd& coefficients,Vector2D p) const{
    const int slot=locate(p);const auto& c=cells[slot];auto v=physicalValue(c.lower,c.size,p,rz,generalTrace,coefficients.col(slot).data());
    essentialJet(v,p);return {{v[0],v[1]},{v[2],v[3]},{v[4],v[5]},c.index};
}
std::vector<HydroCompatibleMomentSpline::NativeVelocity::CellSpeedSamples>
HydroCompatibleMomentSpline::NativeVelocity::cellSpeedSamples() const{
    require(bool(layout_),"moved complete native velocity");
    std::vector<CellSpeedSamples> result;result.reserve(layout_->cells.size());
    for(std::size_t slot=0;slot<layout_->cells.size();slot++){
        const auto& c=layout_->cells[slot];CellSpeedSamples row;row.cell=c.index;
        for(int i=0;i<=8;i++)for(int j=0;j<=8;j++){
            const Vector2D point{c.lower.x+c.size.x*(i/8.),c.lower.y+c.size.y*(j/8.)};
            auto value=physicalValue(c.lower,c.size,point,layout_->rz,layout_->generalTrace,coefficients_.col(slot).data());
            layout_->essentialJet(value,point);const double speed=std::hypot(value[0],value[1]);
            require(std::isfinite(speed),"nonfinite complete-field speed sample");
            row.fineMaximum=std::max(row.fineMaximum,speed);
            if(i%2==0&&j%2==0)row.coarseMaximum=std::max(row.coarseMaximum,speed);
        }
        result.push_back(row);
    }
    return result;
}
HydroCompatibleMomentSpline::FirstJet HydroCompatibleMomentSpline::NativeLayout::transportJet(const Eigen::MatrixXd& coefficients,Vector2D p) const{
    require(std::isfinite(p.x)&&std::isfinite(p.y),"nonfinite native transport query");
    const std::array<bool,2> exterior{p.x<lower.x||p.x>upper.x,p.y<lower.y||p.y>upper.y};
    if(!exterior[0]&&!exterior[1])return firstJet(coefficients,p);
    const Vector2D anchor{std::clamp(p.x,lower.x,upper.x),std::clamp(p.y,lower.y,upper.y)};
    for(int axis=0;axis<2;axis++)if(exterior[axis]){
        const int side=2*axis+((axis?p.y:p.x)>(axis?upper.y:upper.x)?1:0);const auto& spans=outletSpans[side];
        const double transverse=axis?anchor.x:anchor.y;
        const auto found=std::lower_bound(spans.begin(),spans.end(),transverse,[](const auto& ab,double x){return ab[1]<x;});
        require(found!=spans.end()&&transverse>=(*found)[0],"material transport crossed an undeclared outlet");
    }
    require(!rz||p.x>=0,"material transport crossed the RZ axis");
    const int slot=locate(anchor);const auto& c=cells[slot];
    auto v=physicalValue(c.lower,c.size,p,rz,generalTrace,coefficients.col(slot).data(),exterior);essentialJet(v,p);
    require(v.allFinite(),"nonfinite outlet streamfunction continuation");
    return {{v[0],v[1]},{v[2],v[3]},{v[4],v[5]},c.index};
}
Eigen::Matrix<double,6,81> HydroCompatibleMomentSpline::NativeLayout::cellJets(int slot,Vector2D p) const{
    const auto& c=cells.at(slot);return physicalJets(c.lower,c.size,p,rz,generalTrace);
}
Eigen::Matrix<double,3,9> HydroCompatibleMomentSpline::NativeLayout::cartesianAxis(int slot,int axis,double coordinate) const{
    require(!rz,"Cartesian material contraction requires Cartesian velocity jets");
    require(axis==0||axis==1,"invalid Cartesian material axis");
    const auto& c=cells.at(slot);const double h=axis?c.size.y:c.size.x,lower=axis?c.lower.y:c.lower.x;
    auto jets=oneDimensionalJets((coordinate-lower)/h);jets.row(1)/=h;jets.row(2)/=h*h;return jets;
}
InterfaceNormalFluxView HydroCompatibleMomentSpline::NativeLayout::normalFlux(const Eigen::MatrixXd& coefficients,bool velocityProducts) const{
    auto polys=std::make_shared<std::vector<std::array<std::vector<Geometry::UnitPolynomialTerm>,2>>>(faces.size());
    auto charts=std::make_shared<std::vector<InterfaceNormalFluxView::FaceChart>>();charts->reserve(faces.size());
    const auto H=hermite();const auto B=edgeInterior();Eigen::Matrix<double,width,width> powers=Eigen::Matrix<double,width,width>::Zero();powers.topLeftCorner<6,6>()=H;powers.rightCols<3>()=B;
    for(int row=0;row<int(faces.size());row++){
        const auto& f=faces[row];charts->push_back(f.chart);if(f.impermeable)continue;const auto& c=cells[f.cell];
        const double normal=f.horizontal?(f.chart.begin.y-c.lower.y)/c.size.y:(f.chart.begin.x-c.lower.x)/c.size.x;
        require(normal==0||normal==1,"native C1 face is not an actual cell edge");const int normalBasis=normal==0?0:3;
        const double a=f.horizontal?(f.chart.begin.x-c.lower.x)/c.size.x:(f.chart.begin.y-c.lower.y)/c.size.y;
        const double length=f.horizontal?(f.chart.end.x-f.chart.begin.x)/c.size.x:(f.chart.end.y-f.chart.begin.y)/c.size.y;
        for(int mode=0;mode<width;mode++){
            const double coefficient=coefficients(f.horizontal?width*mode+normalBasis:width*normalBasis+mode,f.cell);if(coefficient==0)continue;
            double primitive[width+2]{};for(int k=0;k<width;k++)primitive[k]=powers(k,mode);int degree=width-1;
            if(rz&&!generalTrace){if(f.horizontal){std::fill(std::begin(primitive),std::end(primitive),0.);
                for(int k=0;k<width;k++){primitive[k]+=c.lower.x*c.lower.x*powers(k,mode);primitive[k+1]+=2*c.lower.x*c.size.x*powers(k,mode);primitive[k+2]+=c.size.x*c.size.x*powers(k,mode);}degree=width+1;
            }else for(double& v:primitive)v*=f.chart.begin.x*f.chart.begin.x;}
            for(int half=0;half<2;half++){Geometry::UnitPolynomialTerm term;term.factor=(f.horizontal?-1:1)*.5*length*coefficient*(rz?CartesianMeasure::twoPi:1.);
            for(int k=1;k<=degree;k++)for(int j=0;j<k;j++){
                double choose=1;for(int q=1;q<=j;q++)choose*=double(k-q)/q;
                term.coefficients[j]+=k*primitive[k]*choose*std::pow(a+.5*half*length,k-1-j)*std::pow(.5*length,j);
            }(*polys)[row][half].push_back(term);}
        }
    }
    InterfaceNormalFluxView::Queries q;q.retainedMomentDegree=2;q.retainedVelocityProducts=velocityProducts;q.chart=[charts](int row){return charts->at(row);};
    q.polynomial=[polys](int row,int half){require(half==0||half==1,"invalid native C1 half face");return polys->at(row)[half];};
    q.integral=[poly=q.polynomial](int row,int half,double a,double b,double c,double d){return Geometry::integrateUnitPolynomial(poly(row,half),a,b,c,d);};
    q.density=[poly=q.polynomial](int row,int half,bool end,long double offset){return Geometry::evaluateUnitPolynomial(poly(row,half),double(end?1-offset:offset));};
    q.phaseIntegral=[](int,int,const LevelSetMappedSegment&)->long double{throw std::invalid_argument("native C1 exports parametric phase flux, not legacy geometry");};return InterfaceNormalFluxView(std::move(q));
}
HydroCompatibleMomentSpline::NativeFaceTrace HydroCompatibleMomentSpline::NativeLayout::faceTrace(int row) const{
    constexpr int degree=10;using Poly=std::array<long double,degree+1>;
    const auto& f=faces.at(row);const auto& c=cells[f.cell];NativeFaceTrace out;out.x.resize(degree+1,terms);out.y.resize(degree+1,terms);
    Eigen::Matrix<long double,width,width> powers=Eigen::Matrix<long double,width,width>::Zero();
    powers.topLeftCorner<6,6>()=hermite().cast<long double>();powers.rightCols<3>()=edgeInterior().cast<long double>();
    const bool quotient=rz&&generalTrace&&c.lower.x==0,physicalPsi=rz&&generalTrace&&c.lower.x!=0;
    const double normal=f.horizontal?(f.chart.begin.y-c.lower.y)/c.size.y:(f.chart.begin.x-c.lower.x)/c.size.x;
    auto normalJet=!f.horizontal&&quotient?axisQuotientJets(normal):oneDimensionalJets(normal);
    if(!f.horizontal&&quotient)normalJet/=c.size.x*c.size.x;
    normalJet.row(1)/=f.horizontal?c.size.y:c.size.x;
    if(f.horizontal&&quotient){const auto original=powers;powers.setZero();
        for(int i=0;i<width-2;i++)for(int j=2;j<width;j++)powers(i,j)=original(i+2,j)/(c.size.x*c.size.x);}
    const double origin=f.horizontal?(f.chart.begin.x-c.lower.x)/c.size.x:(f.chart.begin.y-c.lower.y)/c.size.y;
    const double length=f.horizontal?f.chart.end.x-f.chart.begin.x:f.chart.end.y-f.chart.begin.y;
    const double scale=length/(f.horizontal?c.size.x:c.size.y),r=f.chart.begin.x,dr=f.horizontal?length:0;
    const auto restrict=[&](const Poly& p){Poly q{};for(int k=0;k<=degree;k++){long double choose=1;
        for(int j=0;j<=k;j++){if(j)choose=choose*(k-j+1)/j;q[j]+=p[k]*choose*std::pow(static_cast<long double>(origin),k-j)*std::pow(static_cast<long double>(scale),j);}}return q;};
    const auto radial=[&](const Poly& p){Poly q{};require(dr==0||p.back()==0,"native face numerator exceeds degree ten");
        for(int k=0;k<=degree;k++){q[k]+=r*p[k];if(k<degree)q[k+1]+=dr*p[k];}return q;};
    for(int column=0;column<terms;column++){
        const int a=f.horizontal?column/width:column%width,b=f.horizontal?column%width:column/width;
        Poly s{},n{};for(int k=0;k<width;k++){s[k]=powers(k,a)*normalJet(0,b);n[k]=powers(k,a)*normalJet(1,b);}
        s=restrict(s);n=restrict(n);Poly tangent{};for(int k=1;k<=degree;k++)tangent[k-1]=k*s[k]/length;
        const Poly dx=f.horizontal?tangent:n,dy=f.horizontal?n:tangent;Poly ux{},uy{};
        if(physicalPsi){require(r>0,"physical streamfunction trace touches the axis");for(int k=0;k<=degree;k++){ux[k]=dy[k]/r;uy[k]=-dx[k]/r;}}
        else if(rz){ux=radial(dy);uy=radial(dx);for(int k=0;k<=degree;k++)uy[k]=-uy[k]-2*s[k];
            if(f.horizontal&&r>0){ux=radial(ux);uy=radial(uy);for(int k=0;k<=degree;k++){ux[k]/=r;uy[k]/=r;}}}
        else{ux=dy;for(int k=0;k<=degree;k++)uy[k]=-dx[k];}
        for(int k=0;k<=degree;k++){long double vx=0,vy=0,ratio=1;
            for(int j=0;j<=k;j++){if(j)ratio=ratio*(k-j+1)/(degree-j+1);vx+=ux[j]*ratio;vy+=uy[j]*ratio;}
            out.x(k,column)=double(vx);out.y(k,column)=double(vy);
        }
    }return out;
}
HydroCompatibleMomentSpline::NativeFaceProducts HydroCompatibleMomentSpline::NativeLayout::faceProducts(const Eigen::MatrixXd& coefficients,int row) const{
    constexpr int degree=InterfaceParametricTransport::FaceVolumes::productDegree;
    using Poly=std::array<long double,degree+1>;const auto& f=faces.at(row);const auto& c=cells[f.cell];
    NativeFaceProducts result;if(f.impermeable)return result;
    Eigen::Matrix<long double,width,width> powers=Eigen::Matrix<long double,width,width>::Zero();
    powers.topLeftCorner<6,6>()=hermite().cast<long double>();powers.rightCols<3>()=edgeInterior().cast<long double>();
    const bool quotient=rz&&generalTrace&&c.lower.x==0,physicalPsi=rz&&generalTrace&&c.lower.x!=0;
    const double normal=f.horizontal?(f.chart.begin.y-c.lower.y)/c.size.y:(f.chart.begin.x-c.lower.x)/c.size.x;
    auto normalJet=!f.horizontal&&quotient?axisQuotientJets(normal):oneDimensionalJets(normal);
    if(!f.horizontal&&quotient)normalJet/=c.size.x*c.size.x;
    normalJet.row(1)/=f.horizontal?c.size.y:c.size.x;
    if(f.horizontal&&quotient){const auto original=powers;powers.setZero();
        // Same analytic r^2 quotient as physicalJets, including the omitted
        // constrained axis value/first-derivative basis columns.
        for(int i=0;i<width-2;i++)for(int j=2;j<width;j++)powers(i,j)=original(i+2,j)/(c.size.x*c.size.x);
    }
    Poly s{},n{};for(int k=0;k<width;k++)for(int a=0;a<width;a++)for(int b=0;b<width;b++){
        const long double coefficient=coefficients(f.horizontal?width*a+b:width*b+a,f.cell)*powers(k,a);
        s[k]+=coefficient*normalJet(0,b);n[k]+=coefficient*normalJet(1,b);
    }
    const double origin=f.horizontal?(f.chart.begin.x-c.lower.x)/c.size.x:(f.chart.begin.y-c.lower.y)/c.size.y;
    const double length=f.horizontal?f.chart.end.x-f.chart.begin.x:f.chart.end.y-f.chart.begin.y;
    const double scale=length/(f.horizontal?c.size.x:c.size.y);
    const auto restrict=[&](const Poly& p){Poly out{};for(int k=0;k<=degree;k++){long double binomial=1;
        for(int j=0;j<=k;j++){if(j)binomial*=static_cast<long double>(k-j+1)/j;out[j]+=p[k]*binomial*std::pow(static_cast<long double>(origin),k-j)*std::pow(static_cast<long double>(scale),j);}}return out;};
    s=restrict(s);n=restrict(n);Poly tangent{};for(int k=1;k<=degree;k++)tangent[k-1]=k*s[k]/length;
    const Poly dx=f.horizontal?tangent:n,dy=f.horizontal?n:tangent;
    const double r=f.chart.begin.x,dr=f.horizontal?length:0;
    const auto radial=[&](const Poly& p){Poly out{};for(int k=0;k<=degree;k++){out[k]+=r*p[k];if(k<degree)out[k+1]+=dr*p[k];}return out;};
    Poly ux{},uy{};
    if(physicalPsi){const double divisor=f.horizontal?r:f.chart.begin.x;require(divisor>0,"physical streamfunction face quotient touches the axis");
        for(int k=0;k<=degree;k++){ux[k]=dy[k]/divisor;uy[k]=-dx[k]/divisor;}}
    else if(rz){ux=radial(dy);uy=radial(dx);for(int k=0;k<=degree;k++)uy[k]=-uy[k]-2*s[k];
        if(f.horizontal&&r>0){ux=radial(ux);uy=radial(uy);for(int k=0;k<=degree;k++){ux[k]/=r;uy[k]/=r;}}}
    else {ux=dy;for(int k=0;k<=degree;k++)uy[k]=-dx[k];}
    // ux,uy are the numerator of velocity with denominator d=r/r_min
    // only on off-axis horizontal RZ faces. The retained weight is B/d^2.
    Poly px=ux,py=uy;
    if(rz&&f.horizontal&&r>0){px=radial(ux);py=radial(uy);for(int k=0;k<=degree;k++){px[k]/=r;py[k]/=r;}}
    const auto bernstein=[&](const Poly& p,int order){InterfaceParametricTransport::FaceVolumes::Products out{};
        for(int k=0;k<=order;k++){long double value=0,ratio=1;for(int j=0;j<=k;j++){if(j)ratio*=static_cast<long double>(k-j+1)/(order-j+1);value+=p[j]*ratio;}out[k]=double(value);}
        return out;};
    result.x=bernstein(px,degree);result.y=bernstein(py,degree);
    // Multiply bounded Bernstein traces directly. Squaring raw Hermite
    // power coefficients first creates large cancelling energy coefficients.
    constexpr int numeratorDegree=degree/2;
    for(int k=numeratorDegree+1;k<=degree;k++)require(ux[k]==0&&uy[k]==0,"native velocity product exceeds its declared Bernstein space");
    const auto bx=bernstein(ux,numeratorDegree),by=bernstein(uy,numeratorDegree);
    const auto choose=[](int n,int k){long double v=1;for(int j=1;j<=k;j++)v=v*(n-j+1)/j;return v;};
    for(int k=0;k<=degree;k++){long double value=0;
        for(int i=std::max(0,k-numeratorDegree);i<=std::min(k,numeratorDegree);i++){
            const int j=k-i;value+=choose(numeratorDegree,i)*choose(numeratorDegree,j)/choose(degree,k)*(static_cast<long double>(bx[i])*bx[j]+static_cast<long double>(by[i])*by[j]);}
        result.speedSquared[k]=double(value);}
    return result;
}
HydroCompatibleMomentSpline::NativeLayout::Matrix HydroCompatibleMomentSpline::NativeLayout::normalBoundaryConstraints(const HydroCompatibleMomentBoundary& boundary) const{
    require(boundary.usesGeneralTraceSpace()==generalTrace&&boundary.open()==permeable&&boundary.faces().size()==faces.size(),"normal data belong to another native layout");
    if(boundary.nestedNormalConstraints())return boundary.normalConstraints()*momentMap;
    const auto& references=boundary.normalReferences();std::vector<Eigen::Triplet<double>> entries;const MomentSplineBasis::GaussTwelve rule;
    for(int k=0;k<int(references.size());k++){const auto& r=references[k];
        for(const auto& f:faces)if(f.boundary&&f.horizontal==r.horizontal){
            if((f.horizontal?f.chart.begin.y:f.chart.begin.x)!=r.plane)continue;
            const double a=std::max(r.lo,f.horizontal?f.chart.begin.x:f.chart.begin.y),b=std::min(r.hi,f.horizontal?f.chart.end.x:f.chart.end.y);if(!(b>a))continue;
            Eigen::Matrix<double,2,81> local=Eigen::Matrix<double,2,81>::Zero();
            for(int i=0;i<12;i++){const double t=a+.5*(b-a)*(1+rule.x[i]);const Vector2D p=f.horizontal?Vector2D{t,r.plane}:Vector2D{r.plane,t};
                const double weight=.5*(b-a)*rule.w[i]*(rz?CartesianMeasure::twoPi*p.x:1.)/r.area;
                const auto jets=cellJets(f.cell,p);local.row(0)+=weight*jets.row(f.horizontal?1:0);local.row(1)+=weight*(t-r.center)/r.scale*jets.row(f.horizontal?1:0);}
            for(int d=0;d<2;d++)for(int j=0;j<81;j++)if(local(d,j)!=0)entries.emplace_back(2*k+d,81*f.cell+j,local(d,j));
        }
    }
    Matrix local(2*references.size(),polynomial.rows());local.setFromTriplets(entries.begin(),entries.end());return local*polynomial;
}
Eigen::VectorXd HydroCompatibleMomentSpline::NativeLayout::nativeTranspose(const std::vector<PointForce>& forces) const{
    Eigen::VectorXd flat=Eigen::VectorXd::Zero(polynomial.rows());
    for(const auto& f:forces){Eigen::Matrix<double,6,1> force;force<<f.measure.x,f.measure.y,f.dx.x,f.dx.y,f.dy.x,f.dy.y;
        require(force.allFinite(),"nonfinite native C1 work");const int slot=locate(f.point);const auto& c=cells[slot];essentialJet(force,f.point);
        flat.segment<terms>(terms*slot)+=physicalJets(c.lower,c.size,f.point,rz,generalTrace).transpose()*force;
    }
    return polynomial.transpose()*flat;
}
HydroCompatibleMomentSpline::NativeLayout::Matrix HydroCompatibleMomentSpline::NativeLayout::pointWorkGram(const std::vector<PointForce>& forces) const{
    std::vector<Eigen::Triplet<double>> entries;entries.reserve(forces.size()*terms);
    for(int row=0;row<int(forces.size());row++){
        const auto& f=forces[row];Eigen::Matrix<double,6,1> force;
        force<<f.measure.x,f.measure.y,f.dx.x,f.dx.y,f.dy.x,f.dy.y;
        require(force.allFinite(),"nonfinite native C1 metric");const int slot=locate(f.point);const auto& cell=cells[slot];essentialJet(force,f.point);
        const Eigen::Matrix<double,terms,1> local=physicalJets(cell.lower,cell.size,f.point,rz,generalTrace).transpose()*force;
        for(int j=0;j<terms;j++)if(local[j]!=0)entries.emplace_back(row,terms*slot+j,local[j]);
    }
    Matrix flat(forces.size(),polynomial.rows());flat.setFromTriplets(entries.begin(),entries.end());
    const Matrix rows=flat*polynomial;return rows.transpose()*rows;
}
Eigen::VectorXd HydroCompatibleMomentSpline::NativeLayout::transpose(const Eigen::MatrixXd&,const std::vector<PointForce>& forces) const{
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(system.rows());rhs.head(dofs)=nativeTranspose(forces);const auto solution=solve(rhs);
    return data.transpose()*solution.tail(data.rows());
}
}
