#pragma once
#include "Modules/Interface/Compatibility/InterfaceParametricTransport.h"
#include <map>
#include <memory>
#include <tuple>
#include <sstream>
#include <exception>

// Shared transient first-jet world sheet and regular crossing charts. Phase
// volume and surface transport use the SAME trajectory, roots and enclosures.
// These are numerical interpolation/error estimates, not interval certificates.
namespace InterfaceParametricTransport::WorldSheet {
using namespace Quantum::Geometry::Parametric;
inline void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
inline Point tangent(const MapJet& m,Point d){return {m.dx.x*d.x+m.dy.x*d.y,m.dx.y*d.x+m.dy.y*d.y};}
inline double choose(int n,int k){double v=1;for(int j=1;j<=k;++j)v*=double(n-j+1)/j;return v;}
struct Bounds {double lo=std::numeric_limits<double>::infinity(),hi=-std::numeric_limits<double>::infinity();};
inline bool strict(Bounds a){return a.lo>0 || a.hi<0;}
struct Sheet {
    double a=0,b=1,c=0,d=1;
    double errorS=0,errorT=0,positionEnclosure=0;
    // Local material degree eight and local time degree four. Temporary
    // first-jet prediction only; accepted geometry remains the one Curve.
    std::array<std::array<Point,5>,9> p{};
    std::array<std::array<Point,5>,9> controls() const {
        std::array<std::array<Point,5>,9> out{};
        for(int i=0;i<9;++i)for(int j=0;j<5;++j)
            for(int k=0;k<=i;++k)for(int l=0;l<=j;++l)
                out[i][j]+=p[k][l]*(choose(i,k)/choose(8,k)*choose(j,l)/choose(4,l));
        return out;
    }
    Point value(double s,double t,int ds=0,int dt=0) const {
        Point out;for(int i=8;i>=ds;--i)out=out*s+localValue(p[i],t,dt)*(ds?i:1);return out;
    }
    Sheet restrict(double loS,double hiS,double loT,double hiT) const {
        Sheet out;out.a=a+(b-a)*loS;out.b=a+(b-a)*hiS;out.c=c+(d-c)*loT;out.d=c+(d-c)*hiT;
        out.errorS=errorS*(hiS-loS);out.errorT=errorT*(hiT-loT);
        out.positionEnclosure=positionEnclosure;
        std::array<std::array<Point,5>,9> spatial{};
        for(int j=0;j<5;++j){PolynomialArc<9> arc;arc.a=0;arc.b=1;for(int i=0;i<9;++i)arc.p[i]=p[i][j];
            const auto sub=restrictedArc(arc,loS,hiS);for(int i=0;i<9;++i)spatial[i][j]=sub.p[i];}
        for(int i=0;i<9;++i){Arc arc{0,1,spatial[i]};out.p[i]=restrictedArc(arc,loT,hiT).p;}return out;
    }
};
struct SheetBounds {Bounds x,y,xs,ys,xt,yt;};
inline SheetBounds enclose(const Sheet& sheet) {
    const auto b=sheet.controls();SheetBounds out;
    const auto add=[](Bounds& bound,double x){bound.lo=std::min(bound.lo,x);bound.hi=std::max(bound.hi,x);};
    for(int i=0;i<9;++i)for(int j=0;j<5;++j) {
        add(out.x,b[i][j].x);add(out.y,b[i][j].y);
        // Differentiate retained local coefficients before converting to
        // Bernstein. Subtracting global coordinates would erase short arcs.
        if(i<8){Point v;for(int k=0;k<=i;++k)for(int l=0;l<=j;++l)
            v+=sheet.p[k+1][l]*((k+1)*choose(i,k)/choose(7,k)*choose(j,l)/choose(4,l));
            add(out.xs,v.x);add(out.ys,v.y);}
        if(j<4){Point v;for(int k=0;k<=i;++k)for(int l=0;l<=j;++l)
            v+=sheet.p[k][l+1]*((l+1)*choose(i,k)/choose(8,k)*choose(j,l)/choose(3,l));
            add(out.xt,v.x);add(out.yt,v.y);}
    }
    const auto pad=[](Bounds& bound,double error){const double roundoff=128*std::numeric_limits<double>::epsilon()*std::max(std::abs(bound.lo),std::abs(bound.hi));bound.lo-=error+roundoff;bound.hi+=error+roundoff;};
    pad(out.xs,sheet.errorS);pad(out.ys,sheet.errorS);pad(out.xt,sheet.errorT);pad(out.yt,sheet.errorT);
    return out;
}
inline const std::array<std::array<double,5>,5>& timeBasis() {
    static const auto basis=[](){std::array<std::array<double,5>,5> result{};
        for(int i=0;i<5;++i){auto& p=result[i];p[0]=1;int degree=0;
            for(int j=0;j<5;++j)if(j!=i){std::array<double,5> next{};for(int k=0;k<=degree;++k){next[k]-=p[k]*double(j)/(i-j);next[k+1]+=p[k]*4./(i-j);}p=next;++degree;}}
        return result;}();return basis;
}
using Cache=std::map<std::tuple<double,double,double>,MapJet>;
struct Sampler {
    const Curve& old;const Flow& flow;double dt;FluxReport& report;Cache cache;
    const Cache* shared=nullptr;
    MapJet atDuration(Point p,double duration) {
        // Cache the actual argument passed to the flow, not a normalized
        // time reconstructed by division. Source slices and face charts
        // then share identical queries without changing any trajectory.
        const auto key=std::make_tuple(p.x,p.y,duration);const auto found=cache.find(key);
        if(found!=cache.end()){++report.mapCacheHits;return found->second;}
        if(shared){const auto ready=shared->find(key);if(ready!=shared->end()){++report.mapCacheHits;return ready->second;}}
        MapJet value;try{value=flow({p.x,p.y},duration);}catch(const std::exception& e){std::ostringstream message;message<<std::setprecision(17)<<e.what()<<"; phase-flux duration="<<duration;throw std::runtime_error(message.str());}
        cache.emplace(key,value);++report.mapQueries;return value;
    }
    MapJet at(Point p,double t) {return atDuration(p,dt*t);}
    struct Sample {Point p,ds,dt;double f;};
    Sample at(const Sheet& sheet,double s,double t,int axis=0,double plane=0) {
        const double time=t==0?sheet.c:t==1?sheet.d:sheet.c+(sheet.d-sheet.c)*t;
        const auto& arc=old.arcs[old.index(sheet.a)];const double width=arc.b-arc.a;
        const double fromStart=(sheet.a-arc.a)+(sheet.b-sheet.a)*s;
        const double fromEnd=(arc.b-sheet.b)+(sheet.b-sheet.a)*(1-s);
        // Evaluate the local coordinate from its nearer endpoint. Tiny
        // material intervals next to an owned axis pole must not overshoot
        // t=1 through addition; a represented endpoint uses its shared node.
        const double q=fromStart<=fromEnd?fromStart/width:1-fromEnd/width;
        Point p=arc.sampleFraction(q),jet=arc.sampleFraction(q,1);
        auto equation=arc.coordinate(axis);equation[0]-=plane;double relative=polynomial(equation,q);
        if(q==0){p=old.sample(arc.a);jet=old.sample(arc.a,1);relative=(axis?p.y:p.x)-plane;}
        if(q==1) {
            const double endpoint=arc.b==1&&!old.boundary.open()?0:arc.b;
            p=old.sample(endpoint);jet=old.sample(endpoint,1);relative=(axis?p.y:p.x)-plane;
        }
        const auto m=at(p,time);const auto v=tangent(m,jet)*(sheet.b-sheet.a);
        return {{m.point.x,m.point.y},v,{m.timeDerivative.x*dt*(sheet.d-sheet.c),m.timeDerivative.y*dt*(sheet.d-sheet.c)},
            relative+(axis?m.displacement.y:m.displacement.x)};
    }
};

using SliceKey=std::tuple<double,double,double>;
using Slices=std::map<SliceKey,Source>;
inline double estimateSheet(Sheet& sheet,Sampler& sample){
    double maximumError=0;sheet.errorS=sheet.errorT=sheet.positionEnclosure=0;
    for(double s:{0.,.2113248654051871,.7886751345948129,1.})for(double t:{0.,.2113248654051871,.7886751345948129,1.}){
        const auto exact=sample.at(sheet,s,t);const double error=norm(sheet.value(s,t)-exact.p);
        sheet.positionEnclosure=std::max(sheet.positionEnclosure,error);
        if(s>0&&s<1)maximumError=std::max(maximumError,error);
        sheet.errorS=std::max(sheet.errorS,norm(sheet.value(s,t,1)-exact.ds));
        sheet.errorT=std::max(sheet.errorT,norm(sheet.value(s,t,0,1)-exact.dt));
    }return maximumError;
}
struct SheetBuilder {
    const Curve& old;const InterfaceParametricProjection::Grid& grid;Sampler& sample;FluxReport& report;
    Slices slices;const Slices* common=nullptr;std::vector<Sheet> out;
    const Source& slice(double a,double b,double t) {
        using Key=SliceKey;
        const Key key{a,b,t};const auto found=slices.find(key);if(found!=slices.end())return found->second;
        Source curve;curve.boundary=old.boundary;const auto& global=common?*common:slices;const auto full=global.find(Key{0,1,t});
        if(full!=global.end()) {
            const auto knots=full->second.partition(a,b);for(std::size_t k=1;k<knots.size();++k)
                curve.arcs.push_back(restrictedArc(full->second.arcs[full->second.index(knots[k-1])],knots[k-1],knots[k]));
        } else {
            Curve patch;patch.boundary=old.boundary;const auto knots=old.partition(a,b);for(std::size_t k=1;k<knots.size();++k)
                patch.arcs.push_back(restrictedArc(old.arcs[old.index(knots[k-1])],knots[k-1],knots[k]));
            const Flow cachedFlow=[&](Vector2D p,double duration){return sample.atDuration({p.x,p.y},duration);};
            SourceReport diagnostics;curve=transportSource(patch,grid,sample.dt*t,cachedFlow,diagnostics);
        }
        return slices.emplace(key,std::move(curve)).first->second;
    }
    std::vector<Sheet> pieces(double a,double b,double c,double d) {
        std::array<const Source*,5> stage{};std::vector<double> knots;
        for(int i=0;i<5;++i){stage[i]=&slice(a,b,c+(d-c)*(.25*i));for(const auto& arc:stage[i]->arcs){knots.push_back(arc.a);knots.push_back(arc.b);}}
        std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
        std::vector<Sheet> result;
        for(std::size_t k=1;k<knots.size();++k) {
            Sheet sheet;sheet.a=knots[k-1];sheet.b=knots[k];sheet.c=c;sheet.d=d;
            std::array<PolynomialArc<9>,5> arcs;
            for(int i=0;i<5;++i)arcs[i]=restrictedArc(stage[i]->arcs[stage[i]->index(sheet.a)],sheet.a,sheet.b);
            for(int n=0;n<9;++n){sheet.p[n][0]=arcs[0].p[n];
                for(int i=1;i<5;++i)for(int j=0;j<5;++j)sheet.p[n][j]+=(arcs[i].p[n]-arcs[0].p[n])*timeBasis()[i][j];}
            result.push_back(std::move(sheet));
        }
        return result;
    }
    void append(Sheet sheet,int depth) {
        // Shared material nodes retain their own endpoint incidence;
        // refining time cannot remove an old spatial join.
        const double maximumError=estimateSheet(sheet,sample);
        // The sheet only locates regular charts. Refine the unresolved
        // material-time rectangle, not every unrelated part of the loop.
        if(maximumError<=8e-13){report.sheetPositionError=std::max(report.sheetPositionError,maximumError);++report.spaceTimePatches;out.push_back(std::move(sheet));}
        else {
            require(depth<12,"space-time first-jet interpolation unresolved");const double middle=.5*(sheet.c+sheet.d);
            for(auto& child:pieces(sheet.a,sheet.b,sheet.c,middle))append(std::move(child),depth+1);
            for(auto& child:pieces(sheet.a,sheet.b,middle,sheet.d))append(std::move(child),depth+1);
        }
    }
};
inline std::vector<Sheet> worldSheet(const Curve& old,const Source& end,const InterfaceParametricProjection::Grid& grid,Sampler& sample,FluxReport& report) {
    SheetBuilder root{old,grid,sample,report};root.slices.emplace(SliceKey{0,1,1},end);
    Source initial;initial.boundary=old.boundary;for(const auto& arc:old.arcs){PolynomialArc<9> p;p.a=arc.a;p.b=arc.b;std::copy(arc.p.begin(),arc.p.end(),p.p.begin());initial.arcs.push_back(p);}
    root.slices.emplace(SliceKey{0,1,0},std::move(initial));const auto seeds=root.pieces(0,1,0,1);
    struct Part {std::vector<Sheet> sheets;Cache cache;FluxReport report;std::exception_ptr failure;};
    std::vector<Part> parts(seeds.size());
    // The same top-level material intervals own all their temporal children.
    // Each task has private mutable slices/cache and reads immutable common
    // stages. Execution order cannot change polynomial or integration order.
#pragma omp parallel for schedule(dynamic,1) if(seeds.size()>16)
    for(int i=0;i<int(seeds.size());++i){auto& part=parts[i];try{
        Sampler local{old,sample.flow,sample.dt,part.report,{},&sample.cache};
        SheetBuilder builder{old,grid,local,part.report,{},&root.slices};builder.append(seeds[i],0);
        part.sheets=std::move(builder.out);part.cache=std::move(local.cache);
    }catch(...){part.failure=std::current_exception();}}
    for(const auto& part:parts)if(part.failure)std::rethrow_exception(part.failure);
    std::vector<Sheet> out;
    for(auto& part:parts){report.spaceTimePatches+=part.report.spaceTimePatches;report.sheetPositionError=std::max(report.sheetPositionError,part.report.sheetPositionError);
        report.mapQueries+=part.report.mapQueries;report.mapCacheHits+=part.report.mapCacheHits;
        sample.cache.merge(part.cache);for(auto& sheet:part.sheets)out.push_back(std::move(sheet));}
    return out;
}

struct Crossing {double time=0,z=0,dz=0,orientation=0,dtime=0,dependent=0;Point point,materialJet,timeJet;};
struct GraphEnvelope {Bounds z;double timeVariation=0;std::array<double,7> chart{};};
struct Chart {
    std::function<Crossing(double)> curve;
    std::function<GraphEnvelope(double,double)> envelope;
    double a=0,b=1;bool independentTime=false;
};
struct GraphRoot {double coordinate=0;Sampler::Sample value;};
template<class Sample,class Derivative>
GraphRoot bracketedRoot(double lo,double hi,const Sample& sample,const Derivative& derivative) {
    const auto left=sample(lo),right=sample(hi);const double fl=left.f,fr=right.f;
    require((fl<=0&&fr>=0)||(fl>=0&&fr<=0),"world-sheet root is not bracketed");
    GraphRoot best=std::abs(fl)<=std::abs(fr)?GraphRoot{lo,left}:GraphRoot{hi,right};
    if(best.value.f==0)return best;
    double q=std::clamp(lo-fl*(hi-lo)/(fr-fl),lo,hi);
    // Both edge events and interior crossings use this same bracket solve.
    // Every other iteration bisects unless Newton needs an adjacent-float
    // bracket test. Return the best evaluated coordinate, including endpoints,
    // never an unchecked update.
    for(int iteration=0;iteration<2*std::numeric_limits<double>::digits;iteration++){
        const auto value=sample(q);const double residual=value.f;
        if(std::abs(residual)<std::abs(best.value.f))best={q,value};
        if(residual==0)break;
        if((residual<0)==(fl<0))lo=q;else hi=q;
        const double middle=.5*(lo+hi);if(!(middle>lo&&middle<hi))break;
        const double slope=derivative(value),newton=q-residual/slope;
        double next=iteration%2?middle:newton;
        // If Newton rounds back to q, test its adjacent representable
        // coordinate toward the root. Opposite signs there close the SAME
        // bracket to floating-point resolution, without dozens of global
        // bisections. No residual-based early exit or looser root budget.
        if(newton==q && std::isfinite(slope) && slope!=0){
            const double toward=(residual>0)==(slope>0)?lo:hi;
            const double adjacent=std::nextafter(q,toward);
            if(adjacent>lo && adjacent<hi)next=adjacent;
        }
        if(!(next>lo&&next<hi))next=middle;q=next;
    }
    return best;
}
inline std::vector<Chart> regularCharts(const Sheet& sheet,const SheetBounds& bounds,Sampler& sampler,int axis,double plane,FluxReport& report){
    const bool horizontal=axis==1;const auto fs=axis?bounds.ys:bounds.xs,ft=axis?bounds.yt:bounds.xt;
    const bool sGraph=strict(fs),tGraph=strict(ft);require(sGraph||tGraph,"nonregular world-sheet chart");
    const auto normal=[axis](Point p){return axis?p.y:p.x;};const auto transverse=[axis](Point p){return axis?p.x:p.y;};
        // Solve in the direction with the larger guaranteed derivative.
        // At a temporal return this selects material position as the unknown,
        // so ds/dt goes through zero rather than dt/ds becoming infinite.
        // The same criterion also treats spatial tangencies and short arcs.
        const auto strength=[](Bounds b){return strict(b)?std::min(std::abs(b.lo),std::abs(b.hi)):0.;};
        const bool independentTime=sGraph && (!tGraph || strength(fs)>=strength(ft));
        using GraphSample=Sampler::Sample;
        const auto at=[&sampler,sheet,independentTime,axis,plane](double q,double v) {
            const double s=independentTime?v:q,t=independentTime?q:v;
            return sampler.at(sheet,s,t,axis,plane);
        };
        const auto edgePolynomial=[&](double dependent) {
            std::array<double,9> p{};
            if(independentTime)for(int j=0;j<5;++j){double value=0;for(int i=8;i>=0;--i)
                value=value*dependent+(normal(sheet.p[i][j])-(i==0 && j==0?plane:0.));p[j]=value;}
            else for(int i=0;i<9;++i){double value=0;for(int j=4;j>=0;--j)
                value=value*dependent+(normal(sheet.p[i][j])-(i==0 && j==0?plane:0.));p[i]=value;}
            return p;
        };
        std::vector<double> cuts{0,1};std::vector<std::array<double,4>> edgeEvents;
        for(double edge:{0.,1.}) {
            const auto equation=edgePolynomial(edge);std::array<double,9> derivative{};
            for(int k=1;k<9;++k)derivative[k-1]=k*equation[k];
            auto seeds=roots(derivative,0,1);seeds.insert(seeds.begin(),0);seeds.push_back(1);
            std::vector<double> brackets=seeds;
            for(std::size_t k=1;k<seeds.size();++k)brackets.push_back(.5*(seeds[k-1]+seeds[k]));
            std::sort(brackets.begin(),brackets.end());
            const auto slope=[&](double q){const auto value=at(q,edge);return normal(independentTime?value.dt:value.ds);};
            std::vector<double> monotone{0,1};
            for(std::size_t k=1;k<brackets.size();++k) {
                double lo=brackets[k-1],hi=brackets[k],fl=slope(lo),fr=slope(hi);
                if(fl==0)monotone.push_back(lo);if(fr==0)monotone.push_back(hi);
                if(fl*fr>=0)continue;
                for(int iteration=0;iteration<52;++iteration){const double middle=.5*(lo+hi);if(middle==lo || middle==hi)break;const double f=slope(middle);
                    if(f==0){lo=hi=middle;break;}if((f<0)==(fl<0))lo=middle;else hi=middle;}
                monotone.push_back(.5*(lo+hi));
            }
            std::sort(monotone.begin(),monotone.end());monotone.erase(std::unique(monotone.begin(),monotone.end()),monotone.end());
            // Keep the old face-relative polynomial and the integrated
            // displacement separately. An absolute mapped coordinate has a
            // rounding plateau near a tangency and cannot define this root.
            for(std::size_t k=1;k<monotone.size();++k) {
                double lo=monotone[k-1],hi=monotone[k];
                const double fl=at(lo,edge).f,fr=at(hi,edge).f;
                if(fl==0)cuts.push_back(lo);if(fr==0)cuts.push_back(hi);
                if(fl*fr>=0)continue;
                const auto solved=bracketedRoot(lo,hi,[&](double q){return at(q,edge);},
                    [&](const GraphSample& value){return normal(independentTime?value.dt:value.ds);});
                require(std::abs(solved.value.f)<=64*std::numeric_limits<double>::epsilon(),"space-time edge root exceeds the coordinate enclosure");
                edgeEvents.push_back({edge,solved.coordinate,solved.value.f,normal(independentTime?solved.value.dt:solved.value.ds)});
                cuts.push_back(solved.coordinate);
            }
        }
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        const auto solveCrossing=[&sampler,&report,sheet,independentTime,at,edgeEvents,normal,transverse,axis,plane,horizontal](double q) {
            auto left=at(q,0),right=at(q,1);const double fl=left.f,fr=right.f;
            constexpr double coordinate=64*std::numeric_limits<double>::epsilon();
            if(!(fl*fr<=0 || std::abs(fl)<=coordinate || std::abs(fr)<=coordinate)) {
                std::ostringstream message;message<<std::setprecision(17)<<"space-time crossing left its graph chart: axis="<<axis<<" plane="<<plane
                    <<" material="<<sheet.a<<':'<<sheet.b<<" time="<<sheet.c<<':'<<sheet.d<<" time_chart="<<independentTime
                    <<" q="<<q<<" edge_f="<<fl<<':'<<fr;
                for(const auto& event:edgeEvents)message<<" edge_event="<<event[0]<<','<<event[1]<<','<<event[2]<<','<<event[3];
                for(double t:{0.,.25,.5,.75,1.})for(double edge:{0.,1.}){const auto value=at(t,edge);
                    message<<" edge_sample="<<edge<<','<<t<<','<<value.f<<','<<normal(independentTime?value.dt:value.ds);}
                throw std::runtime_error(message.str());
            }
            const bool endpointIncidence=fl*fr>0;
            const GraphRoot solved=endpointIncidence?(std::abs(fl)<=std::abs(fr)?GraphRoot{0,left}:GraphRoot{1,right}):
                bracketedRoot(0.,1.,[&](double v){return at(q,v);},
                    [&](const GraphSample& value){return normal(independentTime?value.ds:value.dt);});
            const double v=solved.coordinate;const auto& value=solved.value;
            const double residual=std::abs(value.f);report.crossingResidual=std::max(report.crossingResidual,residual);
            if(residual>coordinate){std::ostringstream message;message<<std::setprecision(17)<<"space-time intersection exceeds the coordinate enclosure: axis="<<axis<<" plane="<<plane
                <<" material="<<sheet.a<<':'<<sheet.b<<" time="<<sheet.c<<':'<<sheet.d<<" time_chart="<<independentTime
                <<" q="<<q<<" v="<<v<<" residual="<<residual<<" edge_f="<<fl<<':'<<fr;throw std::runtime_error(message.str());}
            const Point dependent=independentTime?value.ds:value.dt,independent=independentTime?value.dt:value.ds;
            const double den=normal(dependent);require(den!=0,"space-time implicit derivative is singular");
            const double dv=-normal(independent)/den,theta=independentTime?q:v;
            return Crossing{sheet.c+(sheet.d-sheet.c)*theta,transverse(value.p),transverse(independent+dependent*dv),
                (horizontal?-1.:1.)*(independentTime?-1.:1.)*std::copysign(1.,den),
                (sheet.d-sheet.c)*(independentTime?1:dv),v,value.p,value.ds,value.dt};
        };
        // Every strip and retained face functional queries the same graph.
        // Keep only completed roots; envelope and quadrature queries share
        // their exact argument, with no interpolation or key quantization.
        auto samples=std::make_shared<std::map<double,Crossing>>();
        const auto crossing=[solveCrossing,samples](double q){
            const auto found=samples->find(q);if(found!=samples->end())return found->second;
            const auto value=solveCrossing(q);samples->emplace(q,value);return value;
        };
        const auto dep=independentTime?fs:ft,ind=independentTime?ft:fs;
        const double slopeBound=std::max(std::abs(ind.lo),std::abs(ind.hi))/strength(dep);
        const auto envelope=[sheet,independentTime,axis,crossing,slopeBound,denominator=strength(dep)](double a,double b) {
            const auto center=crossing(.5*(a+b));
            const double radius=.5*(b-a)*slopeBound+64*std::numeric_limits<double>::epsilon()/denominator;
            const double lo=std::clamp(center.dependent-radius,0.,1.),hi=std::clamp(center.dependent+radius,0.,1.);
            const auto sub=independentTime?sheet.restrict(lo,hi,a,b):sheet.restrict(a,b,lo,hi);
            const auto enclosed=enclose(sub);auto z=axis?enclosed.x:enclosed.y;
            const double pad=sheet.positionEnclosure+64*std::numeric_limits<double>::epsilon();z.lo-=pad;z.hi+=pad;
            return GraphEnvelope{z,(sheet.d-sheet.c)*(b-a)*(independentTime?1:slopeBound),
                {sheet.a,sheet.b,sheet.c,sheet.d,denominator,slopeBound,pad}};
        };
    std::vector<Chart> result;
    for(std::size_t k=1;k<cuts.size();++k){
        const double a=cuts[k-1],b=cuts[k],middle=.5*(a+b);
        const double fl=at(middle,0).f,fr=at(middle,1).f;
        // Half-open dependent charts own a persistent shared-node crossing
        // exactly once. Such a crossing can carry conormal work even when
        // its normal material speed is zero.
        if(!((fl<0&&fr>0)||(fl>0&&fr<0)||(fl==0&&fr!=0)))continue;
        ++report.crossingIntervals;if(independentTime)++report.timeCharts;else ++report.materialCharts;
        result.push_back({crossing,envelope,a,b,independentTime});
    }return result;
}

SurfaceFlux measureSurfaceFlux(const std::vector<Sheet>&,const std::vector<SheetBounds>&,
    Sampler&,const InterfaceParametricProjection::Grid&);
}
