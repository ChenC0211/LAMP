#include "Core/Geometry/ClosedCurveVolumeFit.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SVD>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace LAMP::Geometry::ClosedCurveVolumeFit {
namespace {
using Vec=Eigen::VectorXd;using Sparse=Eigen::SparseMatrix<double>;
using Entry=Eigen::Triplet<double>;
struct Cell {double x0,y0,x1,y1,cap,target;};
struct Derivative {int index;double x,y;};
struct Vertex{double x,y;std::vector<Derivative> gradient;};
// Preserve polygon order and winding. A contiguous arc wholly outside one
// rectangle half-plane can be replaced by its endpoint chord: both the arc
// and that chord clip to nothing. This also works when the rectangle contains
// no polygon vertex, or its intersection has several disconnected pieces.
class PolygonIndex {
    struct Node {int begin,end,left=-1,right=-1;double x0,y0,x1,y1;};
    const std::vector<Vertex>& polygon;std::vector<Node> nodes;
    int build(int begin,int end) {
        const auto& p=polygon[begin];Node node{begin,end,-1,-1,p.x,p.y,p.x,p.y};
        const int index=static_cast<int>(nodes.size());nodes.push_back(node);
        if(end-begin>8) {
            node.left=build(begin,(begin+end)/2);node.right=build((begin+end)/2,end);
            const auto& a=nodes[node.left];const auto& b=nodes[node.right];
            node.x0=std::min(a.x0,b.x0);node.x1=std::max(a.x1,b.x1);
            node.y0=std::min(a.y0,b.y0);node.y1=std::max(a.y1,b.y1);
        } else for(int i=begin+1;i<end;++i) {
            node.x0=std::min(node.x0,polygon[i].x);node.x1=std::max(node.x1,polygon[i].x);
            node.y0=std::min(node.y0,polygon[i].y);node.y1=std::max(node.y1,polygon[i].y);
        }
        nodes[index]=node;return index;
    }
    void select(int index,const Cell& cell,std::vector<Vertex>& out)const {
        const auto& node=nodes[index];
        if(node.x1<cell.x0 || node.x0>cell.x1 || node.y1<cell.y0 || node.y0>cell.y1) {
            out.push_back(polygon[node.begin]);if(node.end-node.begin>1)out.push_back(polygon[node.end-1]);
        } else if(node.left>=0) {select(node.left,cell,out);select(node.right,cell,out);}
        else for(int i=node.begin;i<node.end;++i)out.push_back(polygon[i]);
    }
public:
    explicit PolygonIndex(const std::vector<Vertex>& p):polygon(p){nodes.reserve(p.size()/2+1);build(0,static_cast<int>(p.size()));}
    std::vector<Vertex> select(const Cell& cell)const {std::vector<Vertex> out;select(0,cell,out);return out;}
};
struct Spline{
    struct Support {std::array<int,4> index;std::array<double,4> weight;};
    int count,samples;bool reflected,open;std::vector<Support> basis;
    std::vector<double> knots;std::vector<Vertex> closure;Vec referenceX,referenceY,nx,ny;
    int variable(int i)const{
        if(!open)return reflected?std::min(i,count-i):i;
        if(i==0 || i==count-1)return -1;
        return (reflected?std::min(i,count-1-i):i)-1;
    }
    Support openSupport(double t)const {
        const int span=t==1?count-1:std::clamp(static_cast<int>(std::upper_bound(knots.begin(),knots.end(),t)-knots.begin())-1,3,count-1);
        std::array<double,4> value{1,0,0,0},left{},right{};
        for(int j=1;j<=3;++j){left[j]=t-knots[span+1-j];right[j]=knots[span+j]-t;double saved=0;
            for(int r=0;r<j;++r){const double q=value[r]/(right[r+1]+left[j-r]);value[r]=saved+right[r+1]*q;saved=left[j-r]*q;}
            value[j]=saved;}
        Support result;for(int j=0;j<4;++j){result.index[j]=span-3+j;result.weight[j]=value[j];}return result;
    }
    double parameter(int i)const {return open?(knots[i+1]+knots[i+2]+knots[i+3])/3:double(i)/count;}
    std::array<double,4> roughness(int i)const {
        if(!open)return {-1,3,-3,1};
        // The repeated endpoint knots make plain control-point differences
        // incorrect there. Differentiate the actual clamped spline three
        // times and integrate its piecewise-constant third derivative.
        std::array<std::array<double,4>,4> d{};for(int j=0;j<4;++j)d[j][j]=1;
        for(int order=1;order<=3;++order)for(int j=0;j<4-order;++j){
            const double factor=double(4-order)/(knots[i+j+4]-knots[i+j+order]);
            for(int k=0;k<4;++k)d[j][k]=factor*(d[j+1][k]-d[j][k]);}
        const double scale=std::sqrt(knots[i+4]-knots[i+3])/std::pow(count-3,2.5);
        for(auto& x:d[0])x*=scale;return d[0];
    }
    Spline(int n,int m,bool reflection,bool finite):count(n),samples(m),reflected(reflection),open(finite),basis(m){
        if(open){knots.resize(n+4);for(int i=0;i<n+4;++i)knots[i]=i<=3?0:i>=n?1:double(i-3)/(n-3);}
        for(int k=0;k<m;++k){if(open){basis[k]=openSupport(double(k)/(m-1));continue;}
            double t=double(k)*n/m;int i=int(t);double u=t-i;
            double w[4]={std::pow(1-u,3)/6,(3*u*u*u-6*u*u+4)/6,(-3*u*u*u+3*u*u+3*u+1)/6,u*u*u/6};
            for(int j=0;j<4;++j){basis[k].index[j]=(i+j-1+n)%n;basis[k].weight[j]=w[j];}}}
    std::vector<Vertex> polygon(const Vec& c,bool derivatives)const{std::vector<Vertex> p;p.reserve(samples);
        Vec expanded(count);for(int i=0;i<count;++i)expanded[i]=variable(i)<0?0:c[variable(i)];
        const Vec x=referenceX+nx.cwiseProduct(expanded),y=referenceY+ny.cwiseProduct(expanded);
        for(int k=0;k<samples;++k){Vertex q{0,0,{}};
            for(int j=0;j<4;++j){const int i=basis[k].index[j];const double w=basis[k].weight[j];
                q.x+=w*x[i];q.y+=w*y[i];if(derivatives && w!=0 && variable(i)>=0)q.gradient.push_back({variable(i),w*nx[i],w*ny[i]});}
            std::sort(q.gradient.begin(),q.gradient.end(),[](const auto& a,const auto& b){return a.index<b.index;});
            p.push_back(std::move(q));}
        if(open && reflected){
            // Use one evaluated half for the exact discrete RZ reflection,
            // including derivatives used by the volume equations.
            for(int k=0;k<samples/2;++k){auto q=p[k];q.x=-q.x;for(auto& d:q.gradient)d.x=-d.x;p[samples-1-k]=std::move(q);}
            p[samples/2].x=0;for(auto& d:p[samples/2].gradient)d.x=0;
        }
        p.insert(p.end(),closure.begin(),closure.end());return p;}
};
std::pair<double,std::vector<Derivative>> area(const PolygonIndex& poly,const Cell& cell,bool deriv,bool rz){
    std::vector<Vertex> current=poly.select(cell);
    for(int edge=0;edge<4;++edge){if(current.empty())break;std::vector<Vertex> clipped;clipped.reserve(current.size()+4);
        const bool xAxis=edge<2,above=(edge%2)==0;const double boundary=edge==0?cell.x0:edge==1?cell.x1:edge==2?cell.y0:cell.y1;
        for(std::size_t i=0;i<current.size();++i){const auto& a=current[i];const auto& b=current[(i+1)%current.size()];
            const double av=xAxis?a.x:a.y,bv=xAxis?b.x:b.y;bool ai=above?av>=boundary:av<=boundary,bi=above?bv>=boundary:bv<=boundary;
            if(ai)clipped.push_back(a);if(ai==bi)continue;
            const double t=(boundary-av)/(bv-av);Vertex q{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y),{}};
            if(deriv){std::size_t ia=0,ib=0;
                while(ia<a.gradient.size() || ib<b.gradient.size()) {
                    const int ai=ia<a.gradient.size()?a.gradient[ia].index:std::numeric_limits<int>::max();
                    const int bi=ib<b.gradient.size()?b.gradient[ib].index:std::numeric_limits<int>::max();
                    const int index=std::min(ai,bi);Derivative da{index,0,0},db{index,0,0};
                    if(ai==index)da=a.gradient[ia++];if(bi==index)db=b.gradient[ib++];
                    const double dx=(1-t)*da.x+t*db.x,dy=(1-t)*da.y+t*db.y;
                    const double dt=-(xAxis?dx:dy)/(bv-av);
                    const double gx=xAxis?0:dx+(b.x-a.x)*dt,gy=xAxis?dy+(b.y-a.y)*dt:0;
                    if(gx!=0 || gy!=0)q.gradient.push_back({index,gx,gy});
                }}
            if(xAxis)q.x=boundary;else q.y=boundary;clipped.push_back(std::move(q));
        }current=std::move(clipped);
    }
    double value=0;std::vector<Derivative> gradient;
    for(std::size_t i=0;i<current.size();++i){const auto& a=current[i];const auto& b=current[(i+1)%current.size()];
        const double ax=a.x-cell.x0,ay=a.y-cell.y0,bx=b.x-cell.x0,by=b.y-cell.y0;const double cross=ax*by-ay*bx; const double factor=rz?6.2831853071795864769*(cell.x0/2+(ax+bx)/6):.5; value+=factor*cross;
        if(deriv){const double radial=rz?(6.2831853071795864769/6)*cross:0;
            for(const auto& d:a.gradient)gradient.push_back({d.index,factor*(d.x*by-d.y*bx)+radial*d.x,0});
            for(const auto& d:b.gradient)gradient.push_back({d.index,factor*(ax*d.y-ay*d.x)+radial*d.x,0});}}
    return {value,gradient};
}
Cell cell(const Rectangle& box,bool rz,long double target=0) {
    return {box.lower.x,box.lower.y,box.upper.x,box.upper.y,
        static_cast<double>(PlanarCellCut::capacity(box,rz)),static_cast<double>(target)};
}
bool simplePolygon(const std::vector<Vertex>& polygon)
{
    struct Segment {int index;double lo,hi,ylo,yhi;};std::vector<Segment> segments;
    const int count=static_cast<int>(polygon.size());
    for(int i=0;i<count;++i) {const auto& a=polygon[i];const auto& b=polygon[(i+1)%count];
        if(a.x==b.x && a.y==b.y)return false;
        segments.push_back({i,std::min(a.x,b.x),std::max(a.x,b.x),std::min(a.y,b.y),std::max(a.y,b.y)});}
    std::sort(segments.begin(),segments.end(),[](const auto& a,const auto& b){return a.lo<b.lo;});
    const auto side=[](const Vertex& a,const Vertex& b,const Vertex& p) {
        const long double dx=b.x-a.x,dy=b.y-a.y,x=p.x-a.x,y=p.y-a.y;
        const auto value=dx*y-dy*x,bound=64*std::numeric_limits<double>::epsilon()*(std::abs(dx*y)+std::abs(dy*x));
        return value>bound?1:value< -bound?-1:0;
    };
    std::vector<Segment> active;
    for(const auto& s:segments) {
        active.erase(std::remove_if(active.begin(),active.end(),[&](const auto& a){return a.hi<s.lo;}),active.end());
        for(const auto& other:active) {
            const int distance=std::abs(s.index-other.index);if(distance==1 || distance==count-1)continue;
            if(s.ylo>other.yhi || other.ylo>s.yhi)continue;
            const auto& a=polygon[s.index];const auto& b=polygon[(s.index+1)%count];
            const auto& c=polygon[other.index];const auto& d=polygon[(other.index+1)%count];
            const auto ac=side(a,b,c),ad=side(a,b,d),ca=side(c,d,a),cb=side(c,d,b);
            if(!((ac!=0 && ac==ad) || (ca!=0 && ca==cb)))return false;
        }
        active.push_back(s);
    }
    return true;
}
}
long double measure(const std::vector<Point>& polygon,const Rectangle& box,bool rz)
{
    if(polygon.size()<3)throw std::invalid_argument("[ClosedCurveVolumeFit] empty polygon");
    std::vector<Vertex> vertices;vertices.reserve(polygon.size());
    for(auto p:polygon) vertices.push_back({static_cast<double>(p.x),static_cast<double>(p.y),{}});
    return area(PolygonIndex(vertices),cell(box,rz),false,rz).first;
}
std::vector<long double> measureMany(const std::vector<Point>& polygon,const std::vector<Rectangle>& boxes,bool rz)
{
    if(polygon.size()<3)throw std::invalid_argument("[ClosedCurveVolumeFit] empty polygon");
    std::vector<Vertex> vertices;vertices.reserve(polygon.size());
    for(auto p:polygon)vertices.push_back({static_cast<double>(p.x),static_cast<double>(p.y),{}});
    const PolygonIndex index(vertices);std::vector<long double> result;result.reserve(boxes.size());
    for(const auto& box:boxes)result.push_back(area(index,cell(box,rz),false,rz).first);
    return result;
}
namespace {
Result fitImpl(const std::vector<Point>& trace,const std::vector<Constraint>& constraints,double spacing,bool rz,int samplesPerControl,bool reflected,
    bool open=false,const std::vector<Point>& fixedClosure={})
{
    if(trace.size()<8 || constraints.empty() || !std::isfinite(spacing) || spacing<=0 || samplesPerControl<8)
        throw std::invalid_argument("[ClosedCurveVolumeFit] invalid support");
    std::vector<Vertex> chain;for(auto p:trace) {
        if(!std::isfinite(p.x) || !std::isfinite(p.y))throw std::invalid_argument("[ClosedCurveVolumeFit] nonfinite trace");
        chain.push_back({static_cast<double>(p.x),static_cast<double>(p.y),{}});
    }
    std::vector<Vertex> closure;for(auto p:fixedClosure){
        if(!std::isfinite(p.x) || !std::isfinite(p.y))throw std::invalid_argument("[ClosedCurveVolumeFit] nonfinite boundary closure");
        closure.push_back({static_cast<double>(p.x),static_cast<double>(p.y),{}});}
    // Boundary vertices only close the area integral; they are not curve DOFs.
    auto orientation=chain;orientation.insert(orientation.end(),closure.begin(),closure.end());
    double signedArea=0;const auto origin=orientation.front();
    for(std::size_t i=0;i<orientation.size();++i){const auto& a=orientation[i];const auto& b=orientation[(i+1)%orientation.size()];
        signedArea+=(a.x-origin.x)*(b.y-origin.y)-(a.y-origin.y)*(b.x-origin.x);}
    if(signedArea==0)throw std::invalid_argument("[ClosedCurveVolumeFit] zero oriented area");
    if(open && !simplePolygon(orientation))throw std::invalid_argument("[ClosedCurveVolumeFit] invalid physical boundary closure");
    if(signedArea<0){std::reverse(chain.begin()+(reflected && !open?1:0),chain.end());std::reverse(closure.begin(),closure.end());}
    std::vector<double> arc{0};for(std::size_t i=0;i<chain.size()-(open?1:0);++i){const auto& a=chain[i];const auto& b=chain[(i+1)%chain.size()];
        const double length=std::hypot(b.x-a.x,b.y-a.y);
        if(!(length>0))throw std::invalid_argument("[ClosedCurveVolumeFit] duplicate trace node");
        arc.push_back(arc.back()+length);}
    int n=std::max({16,static_cast<int>((reflected?4:2)*constraints.size()),static_cast<int>(std::ceil(3*arc.back()/spacing))});
    if(reflected && n%2!=(open?1:0))++n;
    const int variables=open?(reflected?(n-1)/2:n-2):reflected?n/2+1:n;
    Spline spline(n,samplesPerControl*n+(open?1:0),reflected,open);spline.closure=closure;
    Sparse interpolation(n,n);std::vector<Entry> entries;Vec x(n),y(n);
    for(int i=0;i<n;++i){const double s=open?arc.back()*spline.parameter(i):arc.back()*i/n;
        const auto k=std::min(static_cast<std::size_t>(std::upper_bound(arc.begin(),arc.end(),s)-arc.begin()-1),chain.size()-(open?2:1));
        const auto& a=chain[k];const auto& b=chain[(k+1)%chain.size()];const double t=(s-arc[k])/(arc[k+1]-arc[k]);
        x[i]=a.x+t*(b.x-a.x);y[i]=a.y+t*(b.y-a.y);
        if(open){const auto basis=spline.openSupport(spline.parameter(i));for(int j=0;j<4;++j)entries.emplace_back(i,basis.index[j],basis.weight[j]);}
        else {entries.emplace_back(i,i,4./6);entries.emplace_back(i,(i+1)%n,1./6);entries.emplace_back(i,(i+n-1)%n,1./6);}}
    // Exact reflection in coefficient space supplies odd radial/even axial
    // parity. Axis endpoints may move along z but never away from r=0.
    const auto reflect=[&](Vec& px,Vec& py) {
        if(!reflected)return;
        if(open){px[n/2]=0;for(int i=0;i<n/2;++i){const int j=n-1-i;const double a=(px[i]-px[j])/2,b=(py[i]+py[j])/2;
            px[i]=a;px[j]=-a;py[i]=py[j]=b;}return;}
        px[0]=px[n/2]=0;
        for(int i=1;i<n/2;++i){const double a=(px[i]-px[n-i])/2,b=(py[i]+py[n-i])/2;
            px[i]=a;px[n-i]=-a;py[i]=py[n-i]=b;}
    };
    reflect(x,y);interpolation.setFromTriplets(entries.begin(),entries.end());
    if(open){Eigen::SparseLU<Sparse> solver;solver.compute(interpolation);if(solver.info()!=Eigen::Success)return {};
        spline.referenceX=solver.solve(x);spline.referenceY=solver.solve(y);}
    else {Eigen::SimplicialLDLT<Sparse> solver(interpolation);
        if(solver.info()!=Eigen::Success)return {};
        spline.referenceX=solver.solve(x);spline.referenceY=solver.solve(y);}
    reflect(spline.referenceX,spline.referenceY);
    if(open){spline.referenceX[0]=chain.front().x;spline.referenceY[0]=chain.front().y;
        spline.referenceX[n-1]=chain.back().x;spline.referenceY[n-1]=chain.back().y;}
    spline.nx.resize(n);spline.ny.resize(n);
    for(int i=0;i<n;++i){const int after=open?std::min(i+1,n-1):(i+1)%n,before=open?std::max(i-1,0):(i+n-1)%n;
        const double dx=spline.referenceX[after]-spline.referenceX[before],
        dy=spline.referenceY[after]-spline.referenceY[before],length=std::hypot(dx,dy);
        if(!(length>0))throw std::invalid_argument("[ClosedCurveVolumeFit] singular reference tangent");
        spline.nx[i]=dy/length;spline.ny[i]=-dx/length;}
    std::vector<Cell> cells;for(const auto& c:constraints) {
        const auto next=cell(c.box,rz,c.insideVolume);
        if(!std::isfinite(next.target) || next.target<0 || next.target>next.cap)
            throw std::invalid_argument("[ClosedCurveVolumeFit] target outside capacity");
        cells.push_back(next);
    }
    // Recover an initial feasible support without changing the reference,
    // normals, objective, inventory or equation set. A missing mixed cell has
    // a genuine trace witness. Place a nearby polygon-edge point at it using
    // only the original local normal displacement DOFs (at most eight).
    // Rank-one straight segments are allowed here if they reproduce the
    // requested displacement; ALL volume rows still face the full-rank gate.
    Vec coefficients=Vec::Zero(variables);Result result;result.controls=n;bool supported=false;
    for(int pass=0;pass<32;++pass) {
        const auto polygon=spline.polygon(coefficients,true);const PolygonIndex index(polygon);int missing=-1;
        for(int row=0;row<static_cast<int>(cells.size());++row){const auto& c=cells[row];
            if(c.target<=0 || c.target>=c.cap)continue;const auto value=area(index,c,false,rz).first;
            if(value<=0 || value>=c.cap){missing=row;break;}}
        if(missing<0){supported=true;break;}
        const auto& c=cells[missing];int witness=-1;double margin=0;
        for(int i=0;i<static_cast<int>(chain.size());++i){const auto& p=chain[i];const double candidate=std::min({p.x-c.x0,c.x1-p.x,p.y-c.y0,c.y1-p.y});
            if(candidate>margin){margin=candidate;witness=i;}}
        if(witness<0)return result;const auto& point=chain[witness];
        const int center=open?static_cast<int>(std::llround(arc[witness]/arc.back()*(spline.samples-1))):
            static_cast<int>(std::llround(arc[witness]/arc.back()*spline.samples))%spline.samples;
        int edge=-1;double fraction=0,best=std::numeric_limits<double>::infinity();
        for(int j=-2*samplesPerControl;j<=2*samplesPerControl;++j){const int k=open?center+j:(center+j+spline.samples)%spline.samples;
            if(open && (k<0 || k>=spline.samples-1))continue;
            const auto& a=polygon[k];const auto& b=polygon[(k+1)%spline.samples];const double dx=b.x-a.x,dy=b.y-a.y;
            const double length=dx*dx+dy*dy;if(!(length>0))continue;
            const double t=std::clamp(((point.x-a.x)*dx+(point.y-a.y)*dy)/length,0.,1.);
            const double distance=std::hypot(point.x-a.x-t*dx,point.y-a.y-t*dy);
            if(distance<best){best=distance;edge=k;fraction=t;}}
        if(edge<0)return result;
        const auto& a=polygon[edge];const auto& b=polygon[(edge+1)%spline.samples];std::vector<Derivative> gradient;
        for(const auto& g:a.gradient)gradient.push_back({g.index,(1-fraction)*g.x,(1-fraction)*g.y});
        for(const auto& g:b.gradient)gradient.push_back({g.index,fraction*g.x,fraction*g.y});
        std::sort(gradient.begin(),gradient.end(),[](const auto& x,const auto& y){return x.index<y.index;});
        std::vector<Derivative> compact;for(const auto& g:gradient){if(!compact.empty() && compact.back().index==g.index){compact.back().x+=g.x;compact.back().y+=g.y;}else compact.push_back(g);}
        if(compact.empty())return result;
        Eigen::MatrixXd matrix(2,compact.size());for(int j=0;j<matrix.cols();++j){matrix(0,j)=compact[j].x;matrix(1,j)=compact[j].y;}
        Eigen::Vector2d offset(point.x-a.x-fraction*(b.x-a.x),point.y-a.y-fraction*(b.y-a.y));
        Eigen::JacobiSVD<Eigen::MatrixXd> correction(matrix,Eigen::ComputeThinU|Eigen::ComputeThinV);correction.setThreshold(1e-11);
        const Vec delta=correction.solve(offset);if((matrix*delta-offset).norm()>1e-11*spacing)return result;
        for(int j=0;j<delta.size();++j)coefficients[compact[j].index]+=delta[j];
        result.maxInitialDisplacement=coefficients.lpNorm<Eigen::Infinity>();++result.initialSupportRepairs;
        if(result.maxInitialDisplacement>.25*spacing)return result;
    }
    if(!supported)return result;
    const int differences=open?n-3:n,fairRows=2*differences;
    Sparse fair(fairRows,variables);entries.clear();Vec fairBase=Vec::Zero(fairRows);
    for(int i=0;i<differences;++i){const auto weights=spline.roughness(i);for(int k=0;k<4;++k){const int j=(i+k)%n;
        if(spline.variable(j)>=0){entries.emplace_back(i,spline.variable(j),weights[k]*spline.nx[j]);entries.emplace_back(differences+i,spline.variable(j),weights[k]*spline.ny[j]);}
        fairBase[i]+=weights[k]*spline.referenceX[j];fairBase[differences+i]+=weights[k]*spline.referenceY[j];}}
    fair.setFromTriplets(entries.begin(),entries.end());
    for(int iteration=0;iteration<16;++iteration) {
        const auto polygon=spline.polygon(coefficients,true);const PolygonIndex index(polygon);
        Sparse jac(cells.size(),variables);Vec rhs(cells.size());entries.clear();
        for(int i=0;i<rhs.size();++i){const auto value=area(index,cells[i],true,rz);
            rhs[i]=(cells[i].target-value.first)/cells[i].cap;
            for(const auto& d:value.second)entries.emplace_back(i,d.index,d.x*spacing/cells[i].cap);}
        jac.setFromTriplets(entries.begin(),entries.end());jac.prune(0.0);
        result.maxFractionResidual=rhs.lpNorm<Eigen::Infinity>();result.iterations=iteration+1;
        // Audit independent volume equations without forming a dense nullspace.
        // QR pivot and SVD singular thresholds are different numerical rank
        // estimates; both require all original rows, with no row dropping.
        Sparse transpose=jac.transpose();double largestRow=0;
        for(int j=0;j<transpose.outerSize();++j)largestRow=std::max(largestRow,transpose.col(j).norm());
        Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> rankSolver;
        rankSolver.setPivotThreshold(1e-11*largestRow);rankSolver.compute(transpose);
        if(rankSolver.info()!=Eigen::Success)return result;
        result.rank=static_cast<int>(rankSolver.rank());
        if(result.rank!=rhs.size()) return result;
        // z=step/spacing, e=A*z+b. The augmented system minimizes ||e||^2
        // subject to J*z=rhs. No normal equations, volume penalty or shift.
        const int m=static_cast<int>(cells.size()),constraintOffset=fairRows+variables,size=constraintOffset+m;entries.clear();
        for(int i=0;i<fairRows;++i)entries.emplace_back(i,i,1);
        for(int j=0;j<fair.outerSize();++j)for(Sparse::InnerIterator a(fair,j);a;++a){
            const int row=static_cast<int>(a.row()),col=static_cast<int>(a.col());
            entries.emplace_back(row,fairRows+col,-a.value());entries.emplace_back(fairRows+col,row,-a.value());}
        for(int j=0;j<jac.outerSize();++j)for(Sparse::InnerIterator a(jac,j);a;++a){
            const int row=static_cast<int>(a.row()),col=static_cast<int>(a.col());
            entries.emplace_back(constraintOffset+row,fairRows+col,a.value());entries.emplace_back(fairRows+col,constraintOffset+row,a.value());}
        Sparse system(size,size);system.setFromTriplets(entries.begin(),entries.end());
        Vec load=Vec::Zero(size);load.head(fairRows)=(fairBase+fair*coefficients)/spacing;load.tail(m)=rhs;
        Eigen::SparseLU<Sparse> constrained;constrained.compute(system);
        if(constrained.info()!=Eigen::Success)return result;
        const Vec solution=constrained.solve(load);
        if(constrained.info()!=Eigen::Success || !solution.allFinite())return result;
        if((system*solution-load).lpNorm<Eigen::Infinity>()>1e-9*std::max(1.,load.lpNorm<Eigen::Infinity>()))return result;
        const Vec scaledStep=solution.segment(fairRows,variables);
        if((jac*scaledStep-rhs).lpNorm<Eigen::Infinity>()>1e-8)return result;
        const Vec step=spacing*scaledStep;
        if(result.maxFractionResidual<1e-13 && step.lpNorm<Eigen::Infinity>()<1e-7*spacing) {result.converged=true;break;}
        bool accepted=false;
        for(double alpha=1;alpha>=1./1024;alpha*=.5) {
            const Vec trial=coefficients+alpha*step;
            // This is an auxiliary local displacement, not permission to
            // change branch identity or search for another liquid component.
            if(trial.lpNorm<Eigen::Infinity>()>.25*spacing)continue;
            const auto trialPolygon=spline.polygon(trial,false);const PolygonIndex trialIndex(trialPolygon);double next=0;
            for(const auto& c:cells)next=std::max(next,std::abs(area(trialIndex,c,false,rz).first-c.target)/c.cap);
            if(next<std::max(1e-13L,result.maxFractionResidual)){coefficients=trial;accepted=true;break;}
        }
        if(!accepted)return result;
    }
    result.maxNormalDisplacement=coefficients.lpNorm<Eigen::Infinity>();
    if(result.converged) {
        const auto polygon=spline.polygon(coefficients,false);
        if(!simplePolygon(polygon)){result.converged=false;return result;}
        for(const auto& p:polygon)result.polygon.push_back({p.x,p.y});
        if(open){result.fixedPolygonVertices={0,static_cast<std::size_t>(spline.samples-1)};
            for(std::size_t i=spline.samples;i<polygon.size();++i)result.fixedPolygonVertices.push_back(i);}
    }
    return result;
}
}
Result fit(const std::vector<Point>& trace,const std::vector<Constraint>& constraints,double spacing,bool rz,int samplesPerControl)
{return fitImpl(trace,constraints,spacing,rz,samplesPerControl,false);}
Result fitAxisToAxis(const std::vector<Point>& trace,const std::vector<Constraint>& constraints,double spacing,int samplesPerControl)
{
    if(trace.size()<8 || trace.front().x!=0 || trace.back().x!=0 || trace.front().y==trace.back().y)
        throw std::invalid_argument("[ClosedCurveVolumeFit] axis trace requires two distinct r=0 endpoints");
    std::vector<Point> reflected=trace;
    for(std::size_t i=1;i+1<trace.size();++i)if(!(trace[i].x>0))
        throw std::invalid_argument("[ClosedCurveVolumeFit] axis trace interior must have positive radius");
    for(std::size_t i=trace.size()-2;i>0;--i)reflected.push_back({-trace[i].x,trace[i].y});
    return fitImpl(reflected,constraints,spacing,true,samplesPerControl,true);
}
Result fitOpen(const std::vector<Point>& trace,const std::vector<Point>& closure,const std::vector<Constraint>& constraints,
    double spacing,bool rz,int samplesPerControl)
{return fitImpl(trace,constraints,spacing,rz,samplesPerControl,false,true,closure);}
Result fitAxisToContact(const std::vector<Point>& input,const std::vector<Point>& inputClosure,const std::vector<Constraint>& constraints,
    double spacing,int samplesPerControl)
{
    if(input.size()<8 || inputClosure.empty() || (input.front().x==0)==(input.back().x==0))
        throw std::invalid_argument("[ClosedCurveVolumeFit] axis-contact trace/closure required");
    auto trace=input,closure=inputClosure;
    if(trace.front().x!=0){std::reverse(trace.begin(),trace.end());std::reverse(closure.begin(),closure.end());}
    for(std::size_t i=1;i<trace.size();++i)if(!(trace[i].x>0))
        throw std::invalid_argument("[ClosedCurveVolumeFit] axis-contact interior must have positive radius");
    std::size_t firstAxis=closure.size();
    for(std::size_t i=0;i<closure.size();++i){
        if(closure[i].x==0 && firstAxis==closure.size())firstAxis=i;
        if(closure[i].x<0 || (firstAxis<closure.size() && closure[i].x!=0))
            throw std::invalid_argument("[ClosedCurveVolumeFit] closure must end along the physical axis");}
    if(firstAxis==closure.size())throw std::invalid_argument("[ClosedCurveVolumeFit] closure has no axis boundary vertex");
    std::vector<Point> reflected;for(std::size_t i=trace.size()-1;i>0;--i)reflected.push_back({-trace[i].x,trace[i].y});
    reflected.insert(reflected.end(),trace.begin(),trace.end());
    closure.resize(firstAxis+1);for(std::size_t i=firstAxis;i>0;--i)closure.push_back({-closure[i-1].x,closure[i-1].y});
    return fitImpl(reflected,constraints,spacing,true,samplesPerControl,true,true,closure);
}
}
