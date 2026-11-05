#include "Modules/Hydrodynamics/HydroCharacteristic.h"
#include "Modules/Interface/HermiteGalsPublication.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <functional>
#include <iomanip>
#include <sstream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace LAMP::Hydrodynamics
{
namespace
{
using Matrix=LevelSetBoundaryTrace::Jacobian;
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[HydroCharacteristic] ")+why);}
void fp()
{
    require(std::fegetround()==FE_TONEAREST,"round-to-nearest required; mode not changed");
#if defined(_M_X64) || defined(__SSE2__)
    require((_mm_getcsr() & ((1U<<15)|(1U<<6)))==0,"gradual underflow required; mode not changed");
#endif
}
bool sameTime(Interface::HermiteGalsTimeLayer a,Interface::HermiteGalsTimeLayer b)
{a.time();b.time();return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
double coord(Vector2D v,int d){return d==0?v.x:v.y;}
// Failure-only context; no diagnostic value participates in a numerical branch.
void diagnosticNumber(std::ostringstream& out,const std::string& key,double value)
{out<<" "<<key<<"="<<std::setprecision(17)<<std::defaultfloat<<value
    <<" "<<key<<"_hex="<<std::hexfloat<<value<<std::defaultfloat;}
void diagnosticVector(std::ostringstream& out,const std::string& key,Vector2D value)
{diagnosticNumber(out,key+"x",value.x);diagnosticNumber(out,key+"y",value.y);}

bool finite(Vector2D v){return std::isfinite(v.x) && std::isfinite(v.y);}
bool finite(Matrix a){return std::isfinite(a.xx)&&std::isfinite(a.xy)&&std::isfinite(a.yx)&&std::isfinite(a.yy);}
Matrix add(Matrix a,Matrix b){return {a.xx+b.xx,a.xy+b.xy,a.yx+b.yx,a.yy+b.yy};}
Matrix scale(Matrix a,double c){return {c*a.xx,c*a.xy,c*a.yx,c*a.yy};}
Matrix outer(Vector2D a,Vector2D b){return {a.x*b.x,a.x*b.y,a.y*b.x,a.y*b.y};}
Matrix multiply(Matrix a,Matrix b){return LevelSetBoundaryTrace::Detail::multiply(a,b);}
Vector2D row(Matrix a,int d){return d==0?Vector2D{a.xx,a.xy}:Vector2D{a.yx,a.yy};}
bool contains(const HydroAdvectorPiece& b,Vector2D p)
{return p.x>=b.lower.x && p.x<=b.upper.x && p.y>=b.lower.y && p.y<=b.upper.y;}
bool onBoundary(const HydroAdvectorPiece& b,Vector2D p)
{return p.x==b.lower.x || p.x==b.upper.x || p.y==b.lower.y || p.y==b.upper.y;}
bool inward(const HydroAdvectorPiece& b,Vector2D p,Vector2D w)
{
    for(int d=0;d<2;++d)if((coord(p,d)==coord(b.lower,d)&&coord(w,d)<0)
        ||(coord(p,d)==coord(b.upper,d)&&coord(w,d)>0))return false;
    return true;
}
HydroCharacteristicJet mixEndpointJets(const HydroAdvectorCellJetPair& endpoints,
    double time,double begin,double end)
{
    const auto& p=endpoints.previous.jet;
    const auto& a=endpoints.arrival.jet;
    const double theta=(time-begin)/(end-begin);
    const auto mix=[&](double left,double right){return time==begin?left:time==end?right:(1-theta)*left+theta*right;};
    HydroCharacteristicJet result;
    result.spatial.value={mix(p.value.x,a.value.x),mix(p.value.y,a.value.y)};
    result.spatial.derivative={mix(p.derivative.xx,a.derivative.xx),mix(p.derivative.xy,a.derivative.xy),mix(p.derivative.yx,a.derivative.yx),mix(p.derivative.yy,a.derivative.yy)};
    result.timeDerivative=(a.value-p.value)*(1/(end-begin));
    require(finite(result.spatial.value)&&finite(result.spatial.derivative)&&finite(result.timeDerivative),"nonfinite actual slab jet");return result;
}
struct Node {Vector2D lo{0,0},hi{0,0};int left=-1,right=-1,begin=0,end=0;};
struct Sample {Vector2D w{0,0},elapsedDerivative{0,0};Matrix jacobian;};
struct Connection {int negative=-1,positive=-1,axis=0;double plane=0,lo=0,hi=0;};
// Roots of a*z^2+b*z+c, with a common scaling and stable small-root formula.
std::vector<double> quadraticRoots(double a,double b,double c)
{
    const double magnitude=std::max({std::abs(a),std::abs(b),std::abs(c)});
    if(magnitude==0)return {};
    a/=magnitude;b/=magnitude;c/=magnitude;
    if(a==0)return b==0?std::vector<double>{}:std::vector<double>{-c/b};
    const double discriminant=std::fma(b,b,-4*a*c);
    if(discriminant<0)return {};
    const double square=std::sqrt(discriminant),q=-.5*(b+std::copysign(square,b));
    if(q==0)return {-b/(2*a)};
    return {q/a,c/q};
}
// Cubic dense curve, not a fitted velocity polynomial. Derivative roots split
// [0,1] into monotone intervals; isolate ALL roots before choosing first exit.
std::vector<double> denseRoots(double a,double b,double c,double d,double minimum=0)
{
    const double scale=std::max({std::abs(a),std::abs(b),std::abs(c),std::abs(d)});
    require(std::isfinite(scale),"nonfinite dense polynomial");
    if(scale==0)return {};a/=scale;b/=scale;c/=scale;d/=scale;
    if(c==0&&d==0)return b==0?std::vector<double>{}:std::vector<double>{-a/b};
    const auto value=[&](double z){return std::fma(std::fma(std::fma(d,z,c),z,b),z,a);};
    const auto derivative=[&](double z){return std::fma(std::fma(3*d,z,2*c),z,b);};
    const auto work=[&](double z){return std::abs(a)+z*(std::abs(b)+z*(std::abs(c)+z*std::abs(d)));};
    const double eps=std::numeric_limits<double>::epsilon();
    std::vector<double> cuts{minimum,1},found;
    for(double z:quadraticRoots(3*d,2*c,b))if(z>minimum&&z<1)cuts.push_back(z);
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    for(double z:cuts)if(std::abs(value(z))<=64*eps*work(z))found.push_back(z);
    for(std::size_t i=1;i<cuts.size();++i) {
        double lo=cuts[i-1],hi=cuts[i],fl=value(lo),fh=value(hi);
        if(fl==0||fh==0||std::signbit(fl)==std::signbit(fh))continue;
        double z=b!=0?-a/b:lo+(hi-lo)*.5;
        if(!(z>lo&&z<hi))z=lo+(hi-lo)*.5;
        bool resolved=false;
        for(int iteration=0;iteration<128;++iteration) {
            const double f=value(z);
            if(f==0){resolved=true;break;}
            if(std::signbit(f)==std::signbit(fl)){lo=z;fl=f;}else{hi=z;fh=f;}
            const double df=derivative(z);double next=df!=0?z-f/df:lo+(hi-lo)*.5;
            if(!(next>lo&&next<hi))next=lo+(hi-lo)*.5;
            if(next==lo||next==hi){z=std::abs(fl)<std::abs(fh)?lo:hi;resolved=true;break;}
            if(std::abs(next-z)<=4*eps*std::max(std::abs(z),std::numeric_limits<double>::min())
                &&std::abs(value(next))<=64*eps*work(next)){z=next;resolved=true;break;}
            z=next;
        }
        require(resolved,"cubic event root did not resolve");found.push_back(z);
    }
    return found;
}

}
Detail::HydroTraceAgreement Detail::hydroQuarticTraceAgreement(const HydroQuarticTrace& left,const HydroQuarticTrace& right)
{
    // Inverse of the degree-four Bernstein evaluation matrix at physical
    // quarter nodes. Work includes every absolute conversion term, including
    // cancellation; there is no fixed absolute floor for tiny velocities.
    constexpr double convert[5][5]={{1,0,0,0,0},{-13./12,4,-3,4./3,-1./4},
        {13./18,-32./9,20./3,-32./9,13./18},{-1./4,4./3,-3,4,-13./12},{0,0,0,0,1}};
    double work=0,jump=0;
    for(int k=0;k<6;++k)for(int b=0;b<5;++b) {
        double a=0,c=0;
        for(int node=0;node<5;++node){const double x=convert[b][node]*left[k][node],y=convert[b][node]*right[k][node];
            a+=x;c+=y;work+=std::abs(x)+std::abs(y);}
        jump=std::max(jump,std::abs(a-c));
    }
    require(std::isfinite(work)&&std::isfinite(jump),"nonfinite connection certificate");
    return {jump<=128*std::numeric_limits<double>::epsilon()*work,work==0?0:jump/work};
}
struct HydroCharacteristicGeometry::Data
{
    explicit Data(HydroAdvectorReconstruction r):reconstruction(std::move(r)){}
    HydroAdvectorReconstruction reconstruction;
    std::vector<int> order;
    std::vector<Node> tree;
    std::vector<Connection> connections;
    std::vector<std::array<std::vector<int>,4>> incident;
};
HydroCharacteristicGeometry::HydroCharacteristicGeometry(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroCharacteristicGeometry HydroCharacteristicGeometry::prepare(const QuadtreeMesh& mesh,HydroAdvectorReconstruction r)
{
    r.momentumGeometry().validateMesh(mesh);auto d=std::make_shared<Data>(std::move(r));
    const auto& pieces=d->reconstruction.pieces();require(!pieces.empty(),"empty actual piece domain");
    d->order.resize(pieces.size());std::iota(d->order.begin(),d->order.end(),0);
    const std::function<int(int,int)> build=[&](int begin,int end) {
        Node node;node.begin=begin;node.end=end;node.lo=pieces[d->order[begin]].lower;node.hi=pieces[d->order[begin]].upper;
        for(int i=begin+1;i<end;++i){const auto& p=pieces[d->order[i]];node.lo={std::min(node.lo.x,p.lower.x),std::min(node.lo.y,p.lower.y)};node.hi={std::max(node.hi.x,p.upper.x),std::max(node.hi.y,p.upper.y)};}
        const int id=static_cast<int>(d->tree.size());d->tree.push_back(node);
        if(end-begin>8){const int dimension=node.hi.x-node.lo.x>=node.hi.y-node.lo.y?0:1,middle=(begin+end)/2;
            std::nth_element(d->order.begin()+begin,d->order.begin()+middle,d->order.begin()+end,[&](int a,int b){return coord(pieces[a].lower+pieces[a].upper,dimension)<coord(pieces[b].lower+pieces[b].upper,dimension);});
            const int left=build(begin,middle),right=build(middle,end);d->tree[id].left=left;d->tree[id].right=right;}
        return id;
    };
    build(0,static_cast<int>(pieces.size()));
    d->incident.resize(pieces.size());
    for(int axis=0;axis<2;++axis) {
        struct Face {int piece;double lo,hi;};
        std::map<double,std::array<std::vector<Face>,2>> planes;
        for(int i=0;i<static_cast<int>(pieces.size());++i)for(int side=0;side<2;++side)
            planes[coord(side?pieces[i].upper:pieces[i].lower,axis)][side].push_back({i,coord(pieces[i].lower,1-axis),coord(pieces[i].upper,1-axis)});
        for(auto& item:planes) {
            auto& lower=item.second[0];auto& upper=item.second[1];
            const auto less=[](const Face& a,const Face& b){return a.lo<b.lo;};
            std::sort(lower.begin(),lower.end(),less);std::sort(upper.begin(),upper.end(),less);
            std::size_t i=0,j=0;
            while(i<upper.size()&&j<lower.size()) {
                const auto a=upper[i],b=lower[j];const double lo=std::max(a.lo,b.lo),hi=std::min(a.hi,b.hi);
                if(lo<hi){const int id=static_cast<int>(d->connections.size());d->connections.push_back({a.piece,b.piece,axis,item.first,lo,hi});
                    d->incident[a.piece][2*axis+1].push_back(id);d->incident[b.piece][2*axis].push_back(id);}
                if(a.hi<=b.hi)++i;if(b.hi<=a.hi)++j;
            }
        }
    }
    return HydroCharacteristicGeometry(std::move(d));
}
const HydroAdvectorReconstruction& HydroCharacteristicGeometry::reconstruction() const noexcept{return data_->reconstruction;}
bool HydroCharacteristicGeometry::sameIdentity(const HydroCharacteristicGeometry& other) const noexcept{return data_==other.data_;}
void HydroCharacteristicGeometry::validateFor(const QuadtreeMesh& mesh) const{reconstruction().momentumGeometry().validateMesh(mesh);}
std::vector<int> HydroCharacteristicGeometry::locatePieces(Vector2D p) const
{
    require(finite(p),"nonfinite location");std::vector<int> stack{0},found;
    while(!stack.empty()) {const int i=stack.back();stack.pop_back();const auto& n=data_->tree[i];
        if(p.x<n.lo.x||p.x>n.hi.x||p.y<n.lo.y||p.y>n.hi.y)continue;
        if(n.left>=0){stack.push_back(n.left);stack.push_back(n.right);}
        else for(int j=n.begin;j<n.end;++j){const int index=data_->order[j];if(contains(reconstruction().pieces()[index],p))found.push_back(index);}}
    std::sort(found.begin(),found.end());return found;
}
struct HydroCharacteristicSource::Data
{
    Data(HydroCharacteristicGeometry g,HydroSpatialAdvector p,HydroSpatialAdvector a):geometry(std::move(g)),previous(std::move(p)),arrival(std::move(a)){}
    HydroCharacteristicGeometry geometry;
    HydroSpatialAdvector previous,arrival;
    std::vector<bool> transparent;
    HydroCharacteristicContinuityReport continuity;
};
struct HydroCharacteristicNativeStarts::Data
{
    struct Slot
    {
        mutable std::once_flag ready;
        mutable HydroAdvectorCellJetPair endpoints;
    };
    struct Entry
    {
        Entry(Vector2D q,std::vector<int> incident)
            :point(q),pieces(std::move(incident)),jets(pieces.size()) {}
        Vector2D point;
        std::vector<int> pieces;
        // Fixed-size construction: once_flags and their results never move.
        std::vector<Slot> jets;
    };
    Data(HydroCharacteristicSource actual,const Interface::HermiteGalsPublication& old)
        :source(std::move(actual)),original(old),sourceIdentity(old.state().identity()),
         sourceClock(old.state().timeLayer()),previousClock(source.previous().timeLayer()),
         arrivalClock(source.arrival().timeLayer()) {}
    HydroCharacteristicSource source;
    Interface::HermiteGalsPublication original;
    Interface::HermiteGalsIdentity sourceIdentity;
    Interface::HermiteGalsTimeLayer sourceClock,previousClock,arrivalClock;
    std::map<std::pair<double,double>,std::shared_ptr<const Entry>> entries;
};
HydroCharacteristicNativeStarts::HydroCharacteristicNativeStarts(std::shared_ptr<const Data> d):data_(std::move(d)){}
std::size_t HydroCharacteristicNativeStarts::pointCount() const
{
    fp();require(static_cast<bool>(data_),"unbound native starts");
    const auto count=data_->entries.size();fp();return count;
}
bool HydroCharacteristicNativeStarts::usesSource(const HydroCharacteristicSource& source) const
{
    fp();require(static_cast<bool>(data_),"unbound native starts");
    const bool same=data_->source.sameIdentity(source);fp();return same;
}
HydroCharacteristicSource::HydroCharacteristicSource(std::shared_ptr<const Data> d):data_(std::move(d)){}
HydroCharacteristicSource HydroCharacteristicGeometry::linear(const QuadtreeMesh& mesh,HydroSpatialAdvector p,HydroSpatialAdvector a) const
{
    validateFor(mesh);const auto& g=reconstruction().momentumGeometry();
    p.validateFor(mesh,g,p.geometryIdentity(),p.timeLayer());a.validateFor(mesh,g,a.geometryIdentity(),a.timeLayer());
    require(p.reconstruction().sameIdentity(reconstruction())&&a.reconstruction().sameIdentity(reconstruction()),"different actual reconstruction");
    require((!p.boundary()&&!a.boundary())||(p.boundary()&&a.boundary()&&p.boundary()->sameIdentity(*a.boundary())),"different endpoint boundary law");
    using Position=Interface::HermiteGalsTimeLayer::Position;const auto pc=p.timeLayer(),ac=a.timeLayer();
    require((pc.position==Position::Initial||pc.position==Position::End)&&ac.position==Position::End
        &&ac.step>pc.step&&ac.step-pc.step==1&&ac.begin==pc.time()&&ac.end>ac.begin,"invalid full endpoint clock interval");
    auto result=std::make_shared<HydroCharacteristicSource::Data>(*this,std::move(p),std::move(a));
    const auto& pieces=reconstruction().pieces();
    for(const auto& edge:data_->connections) {
        const auto& left=pieces[edge.negative];const auto& right=pieces[edge.positive];
        const double hx=std::max(left.upper.x-left.lower.x,right.upper.x-right.lower.x);
        const double hy=std::max(left.upper.y-left.lower.y,right.upper.y-right.lower.y);
        bool transparent=true;double relative=0;
        for(const auto* endpoint:{&result->previous,&result->arrival}) {
            std::array<Detail::HydroQuarticTrace,2> traces{};
            for(int side=0;side<2;++side)for(int node=0;node<5;++node) {
                const double theta=node*.25;
                const double tangent=(1-theta)*edge.lo+theta*edge.hi;
                const Vector2D x=edge.axis==0?Vector2D{edge.plane,tangent}:Vector2D{tangent,edge.plane};
                const auto jet=endpoint->pieceVelocityJet(side?edge.positive:edge.negative,x).jet;
                const auto v=jet.value;const auto j=jet.derivative;
                const std::array<double,6> trace{v.x,v.y,hx*j.xx,hy*j.xy,hx*j.yx,hy*j.yy};
                for(int k=0;k<6;++k)traces[side][k][node]=trace[k];
            }
            const auto certificate=Detail::hydroQuarticTraceAgreement(traces[0],traces[1]);
            transparent=transparent&&certificate.continuous;relative=std::max(relative,certificate.relativeJump);
        }
        result->transparent.push_back(transparent);++result->continuity.connections;
        if(transparent){++result->continuity.transparentConnections;result->continuity.maxTransparentRelativeJump=std::max(result->continuity.maxTransparentRelativeJump,relative);}
    }
    return HydroCharacteristicSource(std::move(result));
}
const HydroCharacteristicGeometry& HydroCharacteristicSource::geometry() const noexcept{return data_->geometry;}
const HydroSpatialAdvector& HydroCharacteristicSource::previous() const noexcept{return data_->previous;}
const HydroSpatialAdvector& HydroCharacteristicSource::arrival() const noexcept{return data_->arrival;}
const HydroCharacteristicContinuityReport& HydroCharacteristicSource::continuityReport() const noexcept{return data_->continuity;}
bool HydroCharacteristicSource::sameIdentity(const HydroCharacteristicSource& other) const noexcept{return data_==other.data_;}
void HydroCharacteristicSource::validateFor(const QuadtreeMesh& mesh,HydroAdvectorIdentity p,HydroAdvectorIdentity a) const
{geometry().validateFor(mesh);require(p==previous().identity()&&a==arrival().identity(),"different actual endpoint identity");}
HydroCharacteristicJet HydroCharacteristicSource::sample(int piece,Vector2D x,double time) const
{
    const double begin=previous().timeLayer().time(),end=arrival().timeLayer().time();
    require(std::isfinite(time)&&time>=begin&&time<=end,"sample time outside actual slab");
    const auto endpoints=previous().pieceVelocityJetPair(arrival(),piece,x);
    return mixEndpointJets(endpoints,time,begin,end);
}
HydroCharacteristicNativeStarts HydroCharacteristicSource::prepareNativeStarts(
    const QuadtreeMesh& mesh,const Interface::HermiteGalsPublication& original) const
{
    fp();original.validateMesh(mesh);geometry().validateFor(mesh);
    const auto& g=geometry().reconstruction().momentumGeometry();
    previous().validateFor(mesh,g,previous().geometryIdentity(),previous().timeLayer());
    arrival().validateFor(mesh,g,arrival().geometryIdentity(),arrival().timeLayer());
    require(previous().reconstruction().sameIdentity(geometry().reconstruction())
        &&arrival().reconstruction().sameIdentity(geometry().reconstruction()),"different native-start reconstruction");
    using Position=Interface::HermiteGalsTimeLayer::Position;
    const auto pc=previous().timeLayer(),ac=arrival().timeLayer();
    require((pc.position==Position::Initial||pc.position==Position::End)&&ac.position==Position::End
        &&ac.step>pc.step&&ac.step-pc.step==1&&ac.begin==pc.time()&&ac.end>ac.begin,
        "invalid native-start endpoint clock interval");
    require(previous().geometryIdentity()==original.state().identity()
        &&sameTime(pc,original.state().timeLayer()),"different original publication for native starts");
    auto prepared=std::make_shared<HydroCharacteristicNativeStarts::Data>(*this,original);
    const auto add=[&](Vector2D q) {
        require(finite(q),"nonfinite native-start point");
        const auto key=std::make_pair(q.x,q.y);
        if(prepared->entries.find(key)==prepared->entries.end())
            prepared->entries.emplace(key,std::make_shared<HydroCharacteristicNativeStarts::Data::Entry>(
                q,geometry().locatePieces(q)));
    };
    const auto vertexPoint=[&](int v){return mesh.verts(mesh.vertIndexToLabel[v]).coord();};
    const auto& primary=original.state().primary();
    for(const auto& v:primary.vertices)add(vertexPoint(v.vertex));
    for(const auto& edge:primary.edges) {
        const auto a=vertexPoint(edge.vertices[0]),b=vertexPoint(edge.vertices[1]);
        add(a+(b-a)*.5);
    }
    for(int cell:primary.geometryCells)add(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord());
    fp();return HydroCharacteristicNativeStarts(std::move(prepared));
}
HydroCharacteristicResult HydroCharacteristicSource::traceBackward(Vector2D point,double at,double until,HydroCharacteristicConfig config,int chosen) const
{return traceImpl(point,at,until,config,nullptr,chosen,false);}
HydroCharacteristicResult HydroCharacteristicSource::traceBackward(Vector2D point,double at,double until,
    HydroCharacteristicConfig config,const HydroCharacteristicNativeStarts& starts,int chosen) const
{return traceImpl(point,at,until,config,&starts,chosen,false);}
HydroCharacteristicResult HydroCharacteristicSource::traceForward(Vector2D point,double at,double until,
    HydroCharacteristicConfig config,int chosen) const
{return traceImpl(point,at,until,config,nullptr,chosen,true);}
HydroCharacteristicResult HydroCharacteristicSource::traceImpl(Vector2D point,double at,double until,
    HydroCharacteristicConfig config,const HydroCharacteristicNativeStarts* starts,int chosen,bool forward) const
{
    fp();
    if(starts)require(starts->data_&&starts->data_->source.sameIdentity(*this),"different native-start source");
    require(finite(point)&&std::isfinite(at)&&std::isfinite(until)
        &&(forward?(at>=previous().timeLayer().time()&&until<=arrival().timeLayer().time()&&at<=until)
                   :(until>=previous().timeLayer().time()&&at<=arrival().timeLayer().time()&&until<=at)),
        "invalid trace point or interval");
    require(chosen>=-1,"invalid initial piece selector");
    require(std::isfinite(config.maxStep)&&config.maxStep>0&&config.maxSubsteps>0&&config.maxEvents>0,"invalid trace work budget");
    const auto& pieces=geometry().reconstruction().pieces();const double total=forward?until-at:at-until;
    const auto diagnosticArrival=point;
    const HydroCharacteristicNativeStarts::Data::Entry* nativeEntry=nullptr;
    if(starts) {
        const auto found=starts->data_->entries.find(std::make_pair(point.x,point.y));
        if(found!=starts->data_->entries.end())nativeEntry=found->second.get();
    }
    HydroCharacteristicResult result;result.point_=point;result.time_=at;
    const auto& topology=*geometry().data_;
    // At a C1 skeleton the incident closed pieces are charts of one field.
    // Resolve only through actual face connections at this point; corner-only
    // contact or a connection without the full C1 certificate is not an alias.
    const auto star=[&](int start,Vector2D x) {
        std::vector<int> found{start};
        for(std::size_t k=0;k<found.size();++k)for(int side=0;side<4;++side)
            for(int id:topology.incident[found[k]][side]) {
                const auto& e=topology.connections[id];const double tangent=coord(x,1-e.axis);
                if(coord(x,e.axis)!=e.plane||tangent<e.lo||tangent>e.hi)continue;
                require(data_->transparent[id],"uncertified skeleton connection");
                const int other=side%2?e.positive:e.negative;
                if(std::find(found.begin(),found.end(),other)==found.end())found.push_back(other);
            }
        std::sort(found.begin(),found.end());return found;
    };
    const auto containingChart=[&](int piece,Vector2D x) {
        if(!contains(pieces[piece],x)) {
            // A tied dense root may put a quadrature point a few ulps across
            // the chart recorded for its span. Evaluate the unchanged point
            // on its actual certified C1 chart; do not clamp the trajectory,
            // extrapolate a polynomial or average velocity/J values.
            const auto& box=pieces[piece];
            const Vector2D anchor{std::clamp(x.x,box.lower.x,box.upper.x),std::clamp(x.y,box.lower.y,box.upper.y)};
            for(int axis=0;axis<2;++axis) {
                try {require(std::abs(coord(x-anchor,axis))<=128*std::numeric_limits<double>::epsilon()
                    *std::max({std::abs(coord(x,axis)),std::abs(coord(anchor,axis)),coord(box.upper-box.lower,axis)}),
                    "quadrature point outside resolved connection bracket");}
                catch(const std::invalid_argument& failure) {
                    std::ostringstream out;out<<failure.what()<<" [bracket piece="<<piece<<" axis="<<axis;
                    diagnosticVector(out,"arrival",diagnosticArrival);diagnosticNumber(out,"at",at);diagnosticNumber(out,"until",until);
                    diagnosticNumber(out,"source_begin",previous().timeLayer().time());diagnosticNumber(out,"source_end",arrival().timeLayer().time());
                    diagnosticVector(out,"x",x);diagnosticVector(out,"anchor",anchor);
                    diagnosticVector(out,"lower",box.lower);diagnosticVector(out,"upper",box.upper);
                    for(int d=0;d<2;++d) {
                        diagnosticNumber(out,"gap"+std::to_string(d),coord(x-anchor,d));
                        diagnosticNumber(out,"tol"+std::to_string(d),128*std::numeric_limits<double>::epsilon()
                            *std::max({std::abs(coord(x,d)),std::abs(coord(anchor,d)),coord(box.upper-box.lower,d)}));
                    }
                    const auto actual=geometry().locatePieces(x);out<<" actual_count="<<actual.size();
                    for(std::size_t i=0;i<actual.size();++i)out<<" actual"<<i<<"="<<actual[i];
                    out<<"]";throw std::invalid_argument(out.str());
                }
            }
            int actual=-1;for(int chart:star(piece,anchor))if(contains(pieces[chart],x)){actual=chart;break;}
            require(actual>=0,"quadrature point outside certified C1 charts");piece=actual;
        }
        return piece;
    };
    const auto sampled=[&](int piece,Vector2D x,double elapsed) {
        require(elapsed>=0&&elapsed<=total,"elapsed time outside trace interval");
        try {
            const int actual=containingChart(piece,x);
            const double time=elapsed==total?until:(forward?at+elapsed:at-elapsed);
            HydroCharacteristicJet j;
            if(nativeEntry&&elapsed==0&&x.x==diagnosticArrival.x&&x.y==diagnosticArrival.y) {
                const auto found=std::lower_bound(nativeEntry->pieces.begin(),nativeEntry->pieces.end(),actual);
                if(found!=nativeEntry->pieces.end()&&*found==actual) {
                    fp();const auto& slot=nativeEntry->jets[static_cast<std::size_t>(found-nativeEntry->pieces.begin())];
                    std::call_once(slot.ready,[&] {
                        slot.endpoints=previous().pieceVelocityJetPair(arrival(),actual,nativeEntry->point);
                    });
                    const double begin=previous().timeLayer().time(),end=arrival().timeLayer().time();
                    require(std::isfinite(time)&&time>=begin&&time<=end,"sample time outside actual slab");
                    j=mixEndpointJets(slot.endpoints,time,begin,end);fp();
                } else j=sample(actual,x,time);
            } else j=sample(actual,x,time);
            // d(direction*u(t))/d elapsed = u_t in either direction.
            return forward?Sample{j.spatial.value,j.timeDerivative,j.spatial.derivative}
                          :Sample{j.spatial.value*(-1),j.timeDerivative,scale(j.spatial.derivative,-1)};
        } catch(const std::invalid_argument& failure) {
            std::ostringstream out;out<<failure.what()<<" [sample piece="<<piece;
            diagnosticNumber(out,"elapsed",elapsed);out<<"]";throw std::invalid_argument(out.str());
        }
    };
    const auto entering=[&](int start,Vector2D x,Vector2D direction) {
        const auto aliases=star(start,x);
        for(int i:aliases)if(contains(pieces[i],x)&&inward(pieces[i],x,direction))return i;
        return -1;
    };
    auto candidates=nativeEntry?nativeEntry->pieces:geometry().locatePieces(point);require(!candidates.empty(),"arrival outside selected domain");
    if(chosen>=0){require(std::find(candidates.begin(),candidates.end(),chosen)!=candidates.end(),"initial piece does not contain arrival");
        require(inward(pieces[chosen],point,sampled(chosen,point,0).w),"initial branch points outside its closure");}
    else {
        std::vector<int> entered;for(int i:candidates)if(inward(pieces[i],point,sampled(i,point,0).w))entered.push_back(i);
        if(entered.empty()) {
            const auto aliases=star(candidates.front(),point);
            require(aliases==candidates,"ambiguous disconnected arrival");
            const int i=candidates.front();const auto w=sampled(i,point,0).w;int side=-1;
            for(int d=0;d<2;++d)for(int sign:{-1,1})if(coord(point,d)==coord(sign<0?pieces[i].lower:pieces[i].upper,d)&&sign*coord(w,d)>0) {
                const int candidateSide=2*d+(sign>0);bool chartSide=false;
                // A physical exit may coincide with a transparent chart edge.
                // Only a certified connection covering this exact point is an
                // internal chart side; it is not a second boundary event.
                for(int id:topology.incident[i][candidateSide]) {
                    const auto& edge=topology.connections[id];const double tangent=coord(point,1-edge.axis);
                    if(coord(point,edge.axis)!=edge.plane||tangent<edge.lo||tangent>edge.hi)continue;
                    require(data_->transparent[id],"uncertified skeleton connection");chartSide=true;
                }
                if(chartSide)continue;
                require(side<0,"simultaneous initial boundary events");side=candidateSide;
            }
            require(side>=0,"unresolved initial boundary direction");result.boundary_=HydroCharacteristicBoundary{side/2,side%2?1:-1,i};result.initialPiece_=result.finalPiece_=i;fp();return result;
        }
        if(entered.size()>1){const auto aliases=star(entered.front(),point);
            for(int i:entered)require(std::binary_search(aliases.begin(),aliases.end(),i),"ambiguous disconnected arrival");}
        chosen=entered.front();
    }
    result.initialPiece_=chosen;result.oneSided_=onBoundary(pieces[chosen],point);
    double elapsed=0;Matrix derivative=Matrix::identity();Vector2D elapsedGradient{0,0};
    int timeSlot=1;double stop=std::min(total,config.maxStep);
    const double eps=std::numeric_limits<double>::epsilon();
    struct Span {int piece=-1,axis=-1;double lo=0,hi=1;};
    struct Walk {int piece=-1,next=-1,side=-1,crossings=0;double fraction=1;std::vector<Span> spans;};
    // Inspect each traversed connection: connectivity through a different
    // route is not evidence that this face has no actual jump.
    const auto neighbor=[&](int piece,int side,Vector2D x,Vector2D direction) {
        int found=-1;
        for(int id:topology.incident[piece][side]) {
            const auto& edge=topology.connections[id];const double t=coord(x,1-edge.axis);
            if(t<edge.lo||t>edge.hi)continue;
            const int other=side%2?edge.positive:edge.negative;
            if(!inward(pieces[other],x,direction))continue;
            if(found>=0){const auto aliases=star(piece,x);
                require(data_->transparent[found]&&data_->transparent[id]
                    &&std::binary_search(aliases.begin(),aliases.end(),other),"ambiguous successor at connection endpoint");}
            else found=id;
        }
        return found;
    };
    const auto walk=[&](int start,Vector2D origin,Vector2D b,Vector2D c,Vector2D d,Vector2D endpoint,bool capture=false,const char* site="linear-stage") {
        Walk out;out.piece=start;double minimum=0;
        // The first nonzero curve coefficient determines departure from a
        // skeleton even when instantaneous velocity is tangent or zero.
        const Vector2D initialDirection{b.x!=0?b.x:c.x!=0?c.x:d.x,b.y!=0?b.y:c.y!=0?c.y:d.y};
        if(!inward(pieces[start],origin,initialDirection)) {
            const int next=entering(start,origin,initialDirection);
            require(next>=0,"unresolved initial dense-curve departure");out.piece=next;
        }
        for(int crossing=0;crossing<=config.maxEvents;++crossing) {
            const auto& box=pieces[out.piece];int side=-1;double first=1;std::vector<int> simultaneous;
            for(int axis=0;axis<2;++axis)for(int sign:{-1,1}) {
                const double plane=coord(sign<0?box.lower:box.upper,axis);
                for(double z:denseRoots(coord(origin,axis)-plane,coord(b,axis),coord(c,axis),coord(d,axis),minimum)) {
                    if(!(z>minimum&&z<=1+64*eps))continue;z=std::min(z,1.);
                    const auto velocity=b+c*(2*z)+d*(3*z*z);
                    // Scale the derivative guard by the terms evaluated at this root,
                    // so a near-start transverse crossing is not classified as grazing.
                    const double normal=coord(velocity,axis),work=std::abs(coord(b,axis))+2*std::abs(coord(c,axis))*std::abs(z)+3*std::abs(coord(d,axis))*z*z;
                    if(std::abs(normal)<=128*eps*work) {
                        auto x=origin+b*z+c*(z*z)+d*(z*z*z);(axis==0?x.x:x.y)=plane;
                        const int id=neighbor(out.piece,2*axis+(sign>0),x,velocity);
                        require(id>=0&&data_->transparent[id],"grazing or unresolved dense-trajectory event");continue;
                    }
                    if(sign*normal<=0)continue;
                    if(side>=0&&std::abs(z-first)<=64*eps) {
                        const int other=2*axis+(sign>0);
                        if(other!=side&&std::find(simultaneous.begin(),simultaneous.end(),other)==simultaneous.end())simultaneous.push_back(other);
                    } else if(side<0||z<first){first=z;side=2*axis+(sign>0);simultaneous.clear();}
                }
            }
            if(side<0){
                try {out.piece=containingChart(out.piece,endpoint);}
                catch(const std::invalid_argument& failure) {
                    std::ostringstream diagnostic;diagnostic<<failure.what()<<" [walk site="<<site<<" start="<<start
                        <<" current="<<out.piece<<" capture="<<capture<<" crossings="<<out.crossings;
                    diagnosticNumber(diagnostic,"minimum",minimum);diagnosticNumber(diagnostic,"elapsed",elapsed);
                    diagnosticVector(diagnostic,"origin",origin);diagnosticVector(diagnostic,"b",b);
                    diagnosticVector(diagnostic,"c",c);diagnosticVector(diagnostic,"d",d);diagnosticVector(diagnostic,"endpoint",endpoint);
                    for(int axis=0;axis<2;++axis)for(int sign:{-1,1}) {
                        const auto key="r"+std::to_string(axis)+(sign<0?"l":"u");
                        const double plane=coord(sign<0?box.lower:box.upper,axis);diagnosticNumber(diagnostic,key+"plane",plane);
                        const auto roots=denseRoots(coord(origin,axis)-plane,coord(b,axis),coord(c,axis),coord(d,axis),minimum);
                        diagnostic<<" "<<key<<"count="<<roots.size();
                        for(std::size_t i=0;i<roots.size();++i)diagnosticNumber(diagnostic,key+std::to_string(i),roots[i]);
                    }
                    diagnostic<<"]";throw std::invalid_argument(diagnostic.str());
                }
                if(capture)out.spans.push_back({out.piece,-1,minimum,1});return out;}
            require(++out.crossings<=config.maxEvents,"event work budget exhausted");
            const int axis=side/2;auto x=origin+b*first+c*(first*first)+d*(first*first*first);
            (axis==0?x.x:x.y)=coord(side%2?box.upper:box.lower,axis);
            const auto velocity=b+c*(2*first)+d*(3*first*first);
            if(!simultaneous.empty()) {
                for(int other:simultaneous)(other/2==0?x.x:x.y)=coord(other%2?box.upper:box.lower,other/2);
                const int next=entering(out.piece,x,velocity);
                require(next>=0,"simultaneous physical-boundary events");
                if(capture)out.spans.push_back({out.piece,first==1?-1:axis,minimum,first});
                if(first==1){
                    // Tied terminal roots can round to opposite sides of the
                    // corner. Keep the actual endpoint and its containing C1
                    // chart, just as for a single terminal connection.
                    for(int chart:star(out.piece,x))if(contains(pieces[chart],endpoint)){out.piece=chart;return out;}
                    require(false,"terminal skeleton has no containing chart");
                }
                out.piece=next;minimum=first;
                continue;
            }
            const int id=neighbor(out.piece,side,x,velocity);
            int cornerSuccessor=-1;
            if(id<0)for(int connection:topology.incident[out.piece][side]) {
                const auto& edge=topology.connections[connection];
                const double tangent=coord(x,1-edge.axis);
                if((tangent!=edge.lo&&tangent!=edge.hi)||!data_->transparent[connection])continue;
                // A tiny tangential displacement can round onto this exact
                // corner while its dense derivative points into a diagonal
                // chart. The single face neighbour then fails its inward
                // test; this remains an internal event in the certified C1
                // star, not an exit from the actual source domain.
                cornerSuccessor=entering(out.piece,x,velocity);
                if(cornerSuccessor>=0)break;
            }
            // A root indistinguishable from the endpoint may lie on either
            // side by rounding. Keep the actual stage/endpoint sample in its
            // containing branch across a certified C1 face. Physical hits
            // retain the snapped boundary point and no interior derivative.
            if(first==1&&id>=0&&data_->transparent[id]&&contains(box,endpoint)
                &&!contains(pieces[side%2?topology.connections[id].positive:topology.connections[id].negative],endpoint)) {
                --out.crossings;if(capture)out.spans.push_back({out.piece,-1,minimum,1});return out;
            }
            if(capture)out.spans.push_back({out.piece,first==1?-1:axis,minimum,first});
            if(cornerSuccessor<0&&(id<0||!data_->transparent[id])) {
                out.side=side;out.fraction=first;
                if(id>=0)out.next=side%2?topology.connections[id].positive:topology.connections[id].negative;
                return out;
            }
            out.piece=cornerSuccessor>=0?cornerSuccessor:
                side%2?topology.connections[id].positive:topology.connections[id].negative;
            minimum=first;
            if(first==1)return out;
        }
        throw std::invalid_argument("[HydroCharacteristic] event work budget exhausted");
    };
    const auto stagePiece=[&](Vector2D x) {
        const auto& current=pieces[chosen];
        // A straight stage segment from the closed rectangle to its strict
        // interior stays in this piece. Boundary endpoints retain event walk.
        if(contains(current,point)&&x.x>current.lower.x&&x.x<current.upper.x
            &&x.y>current.lower.y&&x.y<current.upper.y)return chosen;
        const auto path=walk(chosen,point,x-point,{0,0},{0,0},x);
        if(path.side>=0&&path.fraction<1)return -1;
        return contains(pieces[path.piece],x)?path.piece:-1;
    };
    for(int step=0;step<=config.maxSubsteps;++step) {
        result.finalPiece_=chosen;
        if(elapsed==total) {
            // If the final point is exactly an internal event, its event time
            // still depends on the start. Finish at FIXED end time on
            // the declared outgoing branch, including a zero-duration tail.
            const auto w=sampled(chosen,point,elapsed).w;
            derivative=add(derivative,scale(outer(w,elapsedGradient),-1));
            require(finite(derivative),"nonfinite final map derivative");
            result.point_=point;result.time_=until;result.derivative_=derivative;fp();return result;
        }
        require(step<config.maxSubsteps,"substep work budget exhausted");
        if(elapsed==stop){stop=std::min(total,(++timeSlot)*config.maxStep);require(stop>elapsed,"unresolved fixed time slot");}
        auto first=sampled(chosen,point,elapsed);
        if(!inward(pieces[chosen],point,first.w)) {
            const int next=entering(chosen,point,first.w);
            require(next>=0,"initial branch points outside its closure");
            chosen=next;first=sampled(chosen,point,elapsed);
        }
        const auto& box=pieces[chosen];
        require(inward(box,point,first.w),"initial branch points outside its closure");
        const auto firstDerivative=add(multiply(first.jacobian,derivative),outer(first.elapsedDerivative,elapsedGradient));
        // An event splits an existing time slot. Its remainder ends at the
        // original fixed clock; it does not shift every later RK stage.
        const double remaining=stop-elapsed;double h=remaining;
        Vector2D hGradient=elapsedGradient*(-1);Sample second,third,fourth;
        bool valid=false;
        for(int reduction=0;reduction<128;++reduction) {
            const auto p2=point+first.w*(h/3);
            const int piece2=stagePiece(p2);
            if(piece2>=0) {
                second=sampled(piece2,p2,elapsed+h/3);
                const auto p3=point+second.w*(2*h/3);
                const int piece3=stagePiece(p3);
                if(piece3>=0) {
                    third=sampled(piece3,p3,elapsed+2*h/3);
                    const auto p4=point+third.w*(h*.5);
                    const int piece4=stagePiece(p4);
                    if(piece4>=0){fourth=sampled(piece4,p4,elapsed+h*.5);valid=true;break;}
                }
            }
            h*=.5;hGradient=hGradient*.5;
            if(!(h>0&&elapsed+h>elapsed))throw std::invalid_argument("[HydroCharacteristic] unresolved RK stage at piece boundary: x="+std::to_string(point.x)+", y="+std::to_string(point.y)+", elapsed="+std::to_string(elapsed)+", piece="+std::to_string(chosen));
        }
        require(valid,"RK stages cannot remain across resolved C1 connections");
        const auto stageDerivative=[&](const Sample& atStage,const Sample& prior,Matrix priorDerivative,double fraction) {
            const auto position=add(add(derivative,scale(priorDerivative,h*fraction)),outer(prior.w,hGradient*fraction));
            return add(multiply(atStage.jacobian,position),outer(atStage.elapsedDerivative,elapsedGradient+hGradient*fraction));
        };
        const auto secondDerivative=stageDerivative(second,first,firstDerivative,1./3);
        const auto thirdDerivative=stageDerivative(third,second,secondDerivative,2./3);
        const auto fourthDerivative=stageDerivative(fourth,third,thirdDerivative,.5);
        const auto d2=second.w-first.w,d3=third.w-first.w,d4=fourth.w-first.w;
        const auto a2=d2*2.5-d3+d4*(2./3);auto a3=d2*(-2.5)+d3*1.75-d4*(2./3);
        const auto j2=add(secondDerivative,scale(firstDerivative,-1));
        const auto j3=add(thirdDerivative,scale(firstDerivative,-1));
        const auto j4=add(fourthDerivative,scale(firstDerivative,-1));
        const auto a2Derivative=add(add(scale(j2,2.5),scale(j3,-1)),scale(j4,2./3));
        auto a3Derivative=add(add(scale(j2,-2.5),scale(j3,1.75)),scale(j4,-2./3));
        auto endpoint=first.w+(third.w-first.w)*.75;
        auto endpointDerivative=add(firstDerivative,scale(j3,.75));
        // Complete interior slots use one Picard integral along the quadratic
        // order-two predictor. Integrate each actual Q4 piece's polynomial
        // moments; retain the original Gauss chart selection at rounded ties.
        // Physical/non-C1 hits retain the cubic event map below; no integral
        // samples an extension outside the actual source domain.
        const auto predictor2=d2*1.5;
        const auto predictor2Derivative=scale(j2,1.5);
        const auto predicted=walk(chosen,point,first.w*h,predictor2*h,{0,0},point+(first.w+predictor2)*h,true,"quadratic-predictor");
        if(predicted.side<0) {
            constexpr double nodes[9]={-.96816023950762608984,-.83603110732663579430,-.61337143270059039731,-.32425342340380892904,0,.32425342340380892904,.61337143270059039731,.83603110732663579430,.96816023950762608984};
            constexpr double weights[9]={.08127438836157441197,.18064816069485740406,.26061069640293546232,.31234707704000284007,.33023935500125976316,.31234707704000284007,.26061069640293546232,.18064816069485740406,.08127438836157441197};
            const auto fixedCurveDerivative=[&](double theta) {
                return add(add(derivative,scale(add(scale(firstDerivative,theta),scale(predictor2Derivative,theta*theta)),h)),outer(first.w*theta+predictor2*(theta*theta),hGradient));
            };
            // These internal boundaries join certified C1 pieces. Their
            // Leibniz terms cancel analytically before evaluation, so the
            // derivative is the integral at fixed theta. Differentiating
            // each moving Gauss node/width separately divides by the normal
            // crossing speed and loses this cancellation near tangency.
            // The fixed-theta polynomial moments retain this cancellation;
            // predictor, pieces, point map and physical event rule stay fixed.
            const auto P1=add(scale(firstDerivative,h),outer(first.w,hGradient));
            const auto P2=add(scale(predictor2Derivative,h),outer(predictor2,hGradient));
            const auto curve1=first.w*h,curve2=predictor2*h;
            const double clockWidth=arrival().timeLayer().time()-previous().timeLayer().time();
            const long double beta0=forward?(at+elapsed-previous().timeLayer().time())/clockWidth
                                           :(at-elapsed-previous().timeLayer().time())/clockWidth;
            const long double betaSlope=h/clockWidth,inverseClock=1/clockWidth;
            const std::array<double,4> anchorDerivative{{firstDerivative.xx,firstDerivative.xy,
                firstDerivative.yx,firstDerivative.yy}};
            const std::array<double,2> elapsedDirection{{elapsedGradient.x,elapsedGradient.y}},
                stepDirection{{hGradient.x,hGradient.y}};
            Vector2D mean=first.w;Matrix meanDerivative=firstDerivative;
            for(const auto& span:predicted.spans) {
                const double width=span.hi-span.lo;
                bool recordedChart=true;
                // On an exactly tangent chart boundary retain the original
                // endpoint-basis jet arithmetic: converting its zero normal
                // trace to powers can leave a cancellation-sized departure.
                for(int axis=0;axis<2;++axis)if(coord(curve1,axis)==0&&coord(curve2,axis)==0
                    &&(coord(point,axis)==coord(pieces[span.piece].lower,axis)
                        ||coord(point,axis)==coord(pieces[span.piece].upper,axis)))recordedChart=false;
                for(int q=0;recordedChart&&q<9;++q) {
                    const double z=(1+nodes[q])*.5,theta=(1-z)*span.lo+z*span.hi;
                    const auto x=point+(first.w*theta+predictor2*(theta*theta))*h;
                    if(!contains(pieces[span.piece],x)){recordedChart=false;break;}
                }
                // Polynomial integration belongs to this actual closed piece.
                // Test the same physical quadratic as the endpoint API, also
                // at extrema that the original GL nodes need not encounter.
                for(int axis=0;recordedChart&&axis<2;++axis) {
                    const long double a=coord(point,axis),b=coord(curve1,axis),c=coord(curve2,axis);
                    const long double lo=span.lo,hi=span.hi;
                    const long double lower=coord(pieces[span.piece].lower,axis),
                        upper=coord(pieces[span.piece].upper,axis);
                    const auto inside=[&](long double theta) {
                        const long double x=a+theta*(b+theta*c);
                        return lower<=x&&x<=upper;
                    };
                    recordedChart=inside(lo)&&inside(hi);
                    if(recordedChart&&c!=0) {
                        const long double stationary=-b/(2*c);
                        if(lo<stationary&&stationary<hi)recordedChart=inside(stationary);
                    }
                }
                if(recordedChart) {
                    const auto integrals=previous().pieceCurveIntegralPair(arrival(),
                        span.piece,point,curve1,curve2,derivative,P1,P2,span.lo,span.hi);
                    const auto& p=integrals.previous;const auto& a=integrals.arrival;
                    const auto mixed=[&](long double p0,long double p1,long double a0,long double a1) {
                        return forward?(1-beta0)*p0+beta0*a0+betaSlope*(a1-p1)
                                      :(1-beta0)*p0+beta0*a0+betaSlope*(p1-a1);
                    };
                    // Center in the integral's long-double theta measure before
                    // converting the span difference; constant anchors are not
                    // repeatedly summed or weighted by rounded quadrature nodes.
                    const long double measure=static_cast<long double>(span.hi)-span.lo;
                    std::array<double,2> difference{};
                    std::array<double,4> differenceDerivative{};
                    for(int r=0;r<2;++r)difference[r]=static_cast<double>(
                        forward?mixed(p.V0[r],p.V1[r],a.V0[r],a.V1[r])-measure*coord(first.w,r)
                               :-mixed(p.V0[r],p.V1[r],a.V0[r],a.V1[r])-measure*coord(first.w,r));
                    for(int k=0;k<4;++k) {
                        const int r=k/2,c=k%2;
                        const long double value=(forward?mixed(p.K0[k],p.K1[k],a.K0[k],a.K1[k])
                                                        :-mixed(p.K0[k],p.K1[k],a.K0[k],a.K1[k]))
                            +(a.V0[r]-p.V0[r])*inverseClock*elapsedDirection[c]
                            +(a.V1[r]-p.V1[r])*inverseClock*stepDirection[c];
                        differenceDerivative[k]=static_cast<double>(value-measure*anchorDerivative[k]);
                    }
                    const Vector2D delta{difference[0],difference[1]};
                    const Matrix deltaDerivative{differenceDerivative[0],differenceDerivative[1],
                        differenceDerivative[2],differenceDerivative[3]};
                    require(finite(delta)&&finite(deltaDerivative),"nonfinite piece curve integral");
                    mean=mean+delta;meanDerivative=add(meanDerivative,deltaDerivative);
                } else {
                // A tied root selected a neighbouring actual C1 chart at one
                // of the ORIGINAL points. Keep its original containingChart
                // correction and ordered GL9 arithmetic; no analytic failure
                // is caught or used to choose another integration method.
                for(int q=0;q<9;++q) {
                    const double z=(1+nodes[q])*.5,theta=(1-z)*span.lo+z*span.hi,weight=weights[q]*.5;
                    const auto x=point+(first.w*theta+predictor2*(theta*theta))*h;
                    const auto dx=fixedCurveDerivative(theta);
                    const auto value=[&]() {
                        try {return sampled(span.piece,x,elapsed+h*theta);}
                        catch(const std::invalid_argument& failure) {
                            std::ostringstream out;out<<failure.what()<<" [gl9 q="<<q<<" chosen="<<chosen
                                <<" final_piece="<<predicted.piece<<" span_count="<<predicted.spans.size();
                            diagnosticNumber(out,"elapsed",elapsed);diagnosticNumber(out,"h",h);diagnosticNumber(out,"theta",theta);
                            diagnosticVector(out,"origin",point);diagnosticVector(out,"w",first.w);diagnosticVector(out,"p2",predictor2);
                            diagnosticVector(out,"endpoint",point+(first.w+predictor2)*h);
                            for(std::size_t i=0;i<predicted.spans.size();++i) {
                                const auto& item=predicted.spans[i];const auto key="s"+std::to_string(i);
                                out<<" "<<key<<"piece="<<item.piece<<" "<<key<<"axis="<<item.axis;
                                diagnosticNumber(out,key+"lo",item.lo);diagnosticNumber(out,key+"hi",item.hi);
                                diagnosticVector(out,key+"lower",pieces[item.piece].lower);diagnosticVector(out,key+"upper",pieces[item.piece].upper);
                            }
                            out<<"]";throw std::invalid_argument(out.str());
                        }
                    }();
                    const auto valueDerivative=add(multiply(value.jacobian,dx),outer(value.elapsedDerivative,elapsedGradient+hGradient*theta));
                    const auto difference=value.w-first.w;
                    mean=mean+difference*(width*weight);
                    meanDerivative=add(meanDerivative,scale(add(valueDerivative,scale(firstDerivative,-1)),width*weight));
                }
                }
            }
            const auto correction=mean-endpoint;
            const auto correctionDerivative=add(meanDerivative,scale(endpointDerivative,-1));
            const auto vectorNorm=[](Vector2D v){return std::max(std::abs(v.x),std::abs(v.y));};
            const auto matrixNorm=[](Matrix a){return std::max({std::abs(a.xx),std::abs(a.xy),std::abs(a.yx),std::abs(a.yy)});};
            // A roundoff-sized algebraic correction must not move an exactly
            // represented terminal connection by one ulp. Retain the original
            // point AND its derivative in this branch, using their own scales.
            const double rounding=32*eps*(9*predicted.spans.size()+4);
            if(vectorNorm(correction)>rounding*(vectorNorm(mean)+vectorNorm(endpoint))
                ||h*matrixNorm(correctionDerivative)>rounding*(matrixNorm(derivative)+h*(matrixNorm(meanDerivative)+matrixNorm(endpointDerivative)))) {
                a3=a3+correction;
                a3Derivative=add(a3Derivative,correctionDerivative);
                endpoint=mean;endpointDerivative=meanDerivative;
            }
        }
        const auto path=walk(chosen,point,first.w*h,a2*h,a3*h,point+endpoint*h,false,"final-cubic");
        result.events_+=path.crossings;require(result.events_<=config.maxEvents,"event work budget exhausted");
        chosen=path.piece;
        const double event=h*path.fraction;const int eventSide=path.side;
        if(eventSide<0) {
            point=point+endpoint*h;
            require(contains(pieces[chosen],point),"dense root missed an actual piece exit");
            derivative=add(add(derivative,scale(endpointDerivative,h)),outer(endpoint,hGradient));
            elapsedGradient=elapsedGradient+hGradient;elapsed=h==remaining?stop:elapsed+h;
            if(elapsed==total&&onBoundary(pieces[chosen],point))result.oneSided_=true;
        } else {
            ++result.splitEvents_;const auto& eventBox=pieces[chosen];
            const int direction=eventSide/2,sign=eventSide%2?1:-1;
            const double ratio=event/h;
            const auto denseVelocity=first.w+a2*(2*ratio)+a3*(3*ratio*ratio);
            const auto q=add(add(derivative,scale(firstDerivative,event)),
                add(add(scale(a2Derivative,event*ratio),scale(a3Derivative,event*ratio*ratio)),
                    scale(outer(a2*(ratio*ratio)+a3*(2*ratio*ratio*ratio),hGradient),-1)));
            const auto eventGradient=row(q,direction)*(-1/coord(denseVelocity,direction));
            derivative=add(q,outer(denseVelocity,eventGradient));elapsedGradient=elapsedGradient+eventGradient;
            if(direction==0){derivative.xx=0;derivative.xy=0;}else{derivative.yx=0;derivative.yy=0;}
            point=point+first.w*event+a2*(event*ratio)+a3*(event*ratio*ratio);
            const double plane=coord(sign<0?eventBox.lower:eventBox.upper,direction);
            require(std::abs(coord(point,direction)-plane)<=128*eps*std::max({std::abs(plane),coord(eventBox.upper-eventBox.lower,direction),std::abs(coord(point,direction))}),"unresolved dense event root");
            (direction==0?point.x:point.y)=plane;
            const double next=event==remaining?stop:elapsed+event;
            require(next>elapsed&&next<=total,"unresolved event time advance");elapsed=next;
            if(elapsed==total)result.oneSided_=true;
            if(path.next<0) {
                result.point_=point;result.time_=forward?at+elapsed:at-elapsed;result.boundary_=HydroCharacteristicBoundary{direction,sign,chosen};result.finalPiece_=chosen;result.substeps_=step+1;fp();return result;
            }
            const auto w=sampled(path.next,point,elapsed).w;
            require(sign*coord(w,direction)>0&&inward(pieces[path.next],point,w),"unresolved or grazing successor piece");chosen=path.next;
        }
        require(finite(point)&&finite(derivative)&&finite(elapsedGradient),"nonfinite composite event map");result.substeps_=step+1;
    }
    throw std::invalid_argument("[HydroCharacteristic] substep work budget exhausted");
}
}
