#include "Modules/Interface/Compatibility/InterfaceTemporalPhaseFlux.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/QuarticBernsteinRoots.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <sstream>
#include <iomanip>
#include <set>
#include <type_traits>

namespace InterfaceTemporalPhaseFlux
{
namespace
{
using Patch=Quantum::Geometry::ActiveQuarticHermitePatch;
using Column=Quantum::Geometry::ActiveQuarticHermiteSensitivity;
using CubicColumn=Quantum::Geometry::ActiveHermiteNodalSensitivity;
template<class P> constexpr bool quartic=std::is_same_v<P,Patch>;
using Key=std::tuple<int,int,int>;
using Columns=std::map<Key,long double>;
struct Half { int row, half; Vector2D a,b; };
void require(bool ok,const char* message)
{ if (!ok) throw std::runtime_error(std::string("[InterfaceTemporalPhaseFlux] ")+message); }
double mix(double a,double b,double t)
{
    if (a==b || (std::isnan(a) && std::isnan(b))) return a;
    require(std::isfinite(a) && std::isfinite(b),"time interpolation changed an unused/finite data identity");
    return std::fma(t,b-a,a);
}
template<class P,class Operation> P transform(const std::array<const P*,6>& inputs,Operation operation)
{
    P p=*inputs[0];
    for (const auto* input:inputs) {
        require(input->activeCells==p.activeCells && input->nodes.size()==p.nodes.size()
            && input->cornerSectors.size()==p.cornerSectors.size(),"trajectory changed scalar data layout");
        if constexpr (quartic<P>) require(input->edgeMidpoints.size()==p.edgeMidpoints.size()
            && input->cellCenters.size()==p.cellCenters.size(),"trajectory changed Q4 data layout");
        for (std::size_t k=0;k<p.cornerSectors.size();++k) {
            const auto& a=p.cornerSectors[k];const auto& b=input->cornerSectors[k];
            require(a.cellIndex==b.cellIndex && a.vertexIndex==b.vertexIndex && a.phaseSign==b.phaseSign,
                "time integration requires fixed corner-sector identities");
            require(input->nodes.at(a.vertexIndex).value==0,"time trajectory moved its pinned zero");
        }
    }
    const auto apply=[&](const auto& datum) {
        std::array<double,6> values{};
        for (int k=0;k<6;++k) values[k]=datum(*inputs[k]);
        bool identical=true;
        for (double value:values) identical&=value==values[0] || (std::isnan(value) && std::isnan(values[0]));
        if (identical) return values[0];
        for (double value:values) require(std::isfinite(value),"trajectory mixed finite and unused data");
        const double result=operation(values);
        require(std::isfinite(result),"nonfinite trajectory data");
        return result;
    };
    for (std::size_t i=0;i<p.nodes.size();++i) {
        auto& n=p.nodes[i];
        n.value=apply([&](const P& q){return q.nodes[i].value;});
        n.dx=apply([&](const P& q){return q.nodes[i].dx;});
        n.dy=apply([&](const P& q){return q.nodes[i].dy;});
        n.dxy=apply([&](const P& q){return q.nodes[i].dxy;});
    }
    // Sector jets evolve along the prescribed physical predictor. They are
    // NOT correction unknowns; their endpoint is checked against that predictor.
    for (std::size_t i=0;i<p.cornerSectors.size();++i) {
        auto& n=p.cornerSectors[i].dof;
        n.value=apply([&](const P& q){return q.cornerSectors[i].dof.value;});
        n.dx=apply([&](const P& q){return q.cornerSectors[i].dof.dx;});
        n.dy=apply([&](const P& q){return q.cornerSectors[i].dof.dy;});
        n.dxy=apply([&](const P& q){return q.cornerSectors[i].dof.dxy;});
    }
    if constexpr (quartic<P>) {
    for (std::size_t i=0;i<p.edgeMidpoints.size();++i) {
        for (const auto* input:inputs) require(input->edgeMidpoints[i].vertices==p.edgeMidpoints[i].vertices,
            "trajectory changed edge identity");
        p.edgeMidpoints[i].value=apply([&](const P& q){return q.edgeMidpoints[i].value;});
        p.edgeMidpoints[i].normalDerivative=apply([&](const P& q){return q.edgeMidpoints[i].normalDerivative;});
    }
    for (std::size_t i=0;i<p.cellCenters.size();++i)
        p.cellCenters[i]=apply([&](const P& q){return q.cellCenters[i];});
    }
    return p;
}
template<std::size_t N> long double bernstein(std::array<long double,N> b,long double theta)
{
    for (int size=static_cast<int>(N)-1;size>0;--size)
        for (int i=0;i<size;++i) b[i]=std::fma(theta,b[i+1]-b[i],b[i]);
    return b[0];
}
template<class P> std::array<P,5> controls(const P& old,const P& current,const PhaseTrajectory<P>* trajectory)
{
    std::array<P,5> result;
    result[0]=old;result[4]=current;
    if (!trajectory) {
        const std::array<const P*,6> inputs{{&old,&current,&old,&old,&old,&old}};
        for (int k=1;k<4;++k) result[k]=transform(inputs,[&](const auto& x){return mix(x[0],x[1],.25*k);});
        return result;
    }
    const auto& p=trajectory->samples;
    require(p[4].cornerSectors.size()==current.cornerSectors.size(),"correction changed corner-sector inventory");
    for (std::size_t k=0;k<current.cornerSectors.size();++k) {
        const auto& a=p[4].cornerSectors[k];const auto& b=current.cornerSectors[k];
        require(a.cellIndex==b.cellIndex && a.vertexIndex==b.vertexIndex && a.phaseSign==b.phaseSign
            && a.dof.value==b.dof.value && a.dof.dx==b.dof.dx && a.dof.dy==b.dof.dy && a.dof.dxy==b.dof.dxy,
            "correction changed the prescribed endpoint sector jet");
    }
    const std::array<const P*,6> inputs{{&p[0],&p[1],&p[2],&p[3],&p[4],&current}};
    constexpr long double weights[3][5]={{-13.L/12,4,-3,4.L/3,-.25L},
        {13.L/18,-32.L/9,20.L/3,-32.L/9,13.L/18},{-.25L,4.L/3,-3,4,-13.L/12}};
    for (int k=1;k<4;++k) result[k]=transform(inputs,[&](const auto& x) {
        long double increment=.25L*k*(static_cast<long double>(x[5])-x[4]);
        for (int j=1;j<5;++j) increment+=weights[k-1][j]*(static_cast<long double>(x[j])-x[0]);
        return static_cast<double>(x[0]+increment);
    });
    return result;
}
template<class P> P interpolate(const std::array<P,5>& controls,double t)
{
    const std::array<const P*,6> inputs{{&controls[0],&controls[1],&controls[2],&controls[3],&controls[4],&controls[0]}};
    return transform(inputs,[&](const auto& x){return static_cast<double>(bernstein<5>({{x[0],x[1],x[2],x[3],x[4]}},t));});
}
void finish(const Columns& input,std::vector<Column>& output)
{
    for (const auto& item:input) {
        require(std::isfinite(item.second),"nonfinite integrated phase-flux derivative");
        if (item.second!=0) output.push_back({static_cast<Column::Kind>(std::get<0>(item.first)),std::get<1>(item.first),
                std::get<2>(item.first),item.second});
    }
}
void finish(const Columns& input,std::vector<CubicColumn>& output)
{
    std::map<int,std::array<long double,4>> nodes;
    for (const auto& entry:input) {
        require(std::get<0>(entry.first)==-1 && std::isfinite(entry.second),"invalid cubic flux column");
        nodes[std::get<1>(entry.first)][std::get<2>(entry.first)]+=entry.second;
    }
    for (const auto& entry:nodes) output.push_back({entry.first,entry.second});
}
void accumulate(Columns& target,const std::vector<Column>& columns,long double factor)
{
    for (const auto& c:columns) target[{static_cast<int>(c.kind),c.index,c.component}]+=factor*c.coefficient;
}
void accumulate(Columns& target,const std::vector<CubicColumn>& columns,long double factor)
{
    for (const auto& c:columns) for (int k=0;k<4;++k) target[{-1,c.vertexIndex,k}]+=factor*c.coefficients[k];
}
template<class P> auto valueVariation(const LevelSetMappedSegment& trace,long double position)
{
    constexpr auto fixed=Quantum::Geometry::CornerSectorVariation::HoldFixed;
    if constexpr (quartic<P>) return trace.nativeQuarticValueVariation(false,position,fixed);
    else return trace.nativeCubicValueVariation(false,position,fixed);
}
template<class P> bool sameZeros(const LevelSetMappedSegment& a,const LevelSetMappedSegment& b)
{
    if constexpr (quartic<P>) return a.sameNativeQuarticTraceZeros(b);
    else return a.sameNativeCubicTraceZeros(b);
}
template<class P> auto phaseVariation(const LevelSetMappedSegment& trace)
{
    constexpr auto fixed=Quantum::Geometry::CornerSectorVariation::HoldFixed;
    if constexpr (quartic<P>) return trace.nativeQuarticNegativePhaseVariation(fixed);
    else return trace.nativeNegativePhaseVariation(fixed);
}
struct Rule
{
    std::vector<long double> flux,absolute;
    std::vector<Columns> derivative;
    explicit Rule(int faces):flux(faces,0),absolute(faces,0),derivative(faces){}
};
template<class P>
auto integrateImpl(const QuadtreeMesh& mesh,const InterfaceTransportDomain& transport,
    const P& oldPatch,const P& currentPatch,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& oldInput,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& currentInput,
    const InterfaceNormalFluxView& oldVelocity,const InterfaceNormalFluxView& currentVelocity,
    double dt,const Options& options,const PhaseTrajectory<P>* trajectory)
{
    const auto& domain=transport.geometry(); domain.validateMeshIdentity(mesh);
    require(oldInput && currentInput,"two scalar sources are required");
    InterfacePhaseGeometry::validateLayout(mesh,domain,*oldInput);
    InterfacePhaseGeometry::validateLayout(mesh,domain,*currentInput);
    const auto nativeSource=[&](const P& patch,const auto& input) {
        return input->usesNativeLeaves()?input:std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
            LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,domain.axisymmetric()));
    };
    const auto oldSource=nativeSource(oldPatch,oldInput),currentSource=nativeSource(currentPatch,currentInput);
    InterfacePhaseGeometry::validateLayout(mesh,domain,*oldSource);
    InterfacePhaseGeometry::validateLayout(mesh,domain,*currentSource);
    require(std::isfinite(dt) && dt>0 && std::isfinite(options.relativeTolerance) && options.relativeTolerance>0
        && std::isfinite(options.capacityTolerance) && options.capacityTolerance>0
        && options.maxDepth>=0 && options.maxEvaluations>=12,"invalid time integration controls");
    const auto temporal=controls(oldPatch,currentPatch,trajectory);
    std::array<std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>,5> controlSource;
    controlSource[0]=oldSource;controlSource[4]=currentSource;
    for (int k=1;k<4;++k) controlSource[k]=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,temporal[k],domain.axisymmetric()));
    const int faceCount=domain.nativeFaceCount(),cellCount=mesh.cellIndexNum();
    const long double nan=std::numeric_limits<long double>::quiet_NaN();
    std::conditional_t<quartic<P>,Result,CubicResult> result;
    auto& integral=result.integral; auto& budget=integral.volume;
    budget.positiveFluxByNativeFace.assign(faceCount,nan);
    budget.outwardFluxByCell.assign(cellCount,nan);
    integral.absoluteTransportByCell.assign(cellCount,nan);
    integral.integrationUncertaintyByCell.assign(cellCount,nan);
    integral.errorByNativeFace.assign(faceCount,nan);
    result.derivative.positiveFluxByNativeFace.resize(faceCount);
    result.derivative.outwardFluxByCell.resize(cellCount);
    std::vector<long double> capacities(cellCount,nan),faceTolerance(faceCount,0),absolute(faceCount,0);
    std::vector<Columns> derivative(faceCount),cellDerivative(cellCount);
    for (std::size_t row=0;row<domain.cellIndexByRow().size();++row) {
        const int c=domain.cellIndexByRow()[row]; capacities[c]=domain.cellVolumes()[row];
        budget.outwardFluxByCell[c]=integral.absoluteTransportByCell[c]=0;
        integral.integrationUncertaintyByCell[c]=0;
    }
    std::vector<Half> dynamic;
    std::vector<double> events{0,1};
    const std::array<double,8> x8{{-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,
        .1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363}};
    const std::array<double,8> w8{{.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,
        .3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763}};
    // A strict whole-cell certificate covers each incident half-face. Cache
    // it once per temporal control cell; pure regions need no root inventory.
    std::vector<int> commonCellSign(cellCount,2);
    const auto cellSign=[&](int cell) {
        if (cell<0) return 0;
        auto& sign=commonCellSign[cell];
        if (sign!=2) return sign;
        const auto center=mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord();
        sign=controlSource[0]->closedCellSign(center);
        for (int k=1;k<5 && sign!=0;++k) if (controlSource[k]->closedCellSign(center)!=sign) sign=0;
        return sign;
    };
    for (const auto& face:domain.faces(mesh)) {
        if (face.negativeCell<0 && face.positiveCell<0) continue;
        const auto oldChart=oldVelocity.faceChart(face.nativeRow),chart=currentVelocity.faceChart(face.nativeRow);
        const auto same=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
        require(same(oldChart.begin,chart.begin) && same(oldChart.midpoint,chart.midpoint)
            && same(oldChart.end,chart.end),"velocity time layers have different face charts");
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal?a.x>b.x:a.y>b.y) std::swap(a,b);
        require(same(chart.begin,a) && same(chart.end,b) && same(chart.midpoint,face.center),
            "time flux does not use the actual native face");
        auto& flux=budget.positiveFluxByNativeFace[face.nativeRow]; flux=0;
        integral.errorByNativeFace[face.nativeRow]=0;
        const auto capacity=std::min(face.negativeCell>=0?capacities[face.negativeCell]:std::numeric_limits<long double>::infinity(),
            face.positiveCell>=0?capacities[face.positiveCell]:std::numeric_limits<long double>::infinity());
        faceTolerance[face.nativeRow]=options.capacityTolerance*capacity;
        const bool closed=transport.boundary(face.nativeRow)==InterfaceTransportDomain::Boundary::StationaryWall
            || (domain.axisymmetric() && !face.horizontal && face.center.x==0);
        int uniformCellSign=closed?0:cellSign(face.negativeCell);
        if (!closed && uniformCellSign==0) uniformCellSign=cellSign(face.positiveCell);
        for (int half=0;half<2;++half) {
            const auto begin=half?chart.midpoint:chart.begin,end=half?chart.end:chart.midpoint;
            if (closed) {
                require(oldVelocity.integratedNormalFlux(face.nativeRow,half,0,1)==0
                    && currentVelocity.integratedNormalFlux(face.nativeRow,half,0,1)==0,
                    "wall/axis time flux is not zero");
                continue;
            }
            if (uniformCellSign!=0) {
                if (uniformCellSign<0) {
                    const long double old=oldVelocity.integratedNormalFlux(face.nativeRow,half,0,1);
                    const long double current=currentVelocity.integratedNormalFlux(face.nativeRow,half,0,1);
                    flux+=.5L*dt*(old+current);absolute[face.nativeRow]+=.5L*dt*(std::abs(old)+std::abs(current));
                }
                continue;
            }
            const auto sign=[](const LevelSetMappedSegment& phase) {
                if (!phase.roots().empty() || phase.ranges().empty()) return 0;
                const bool negative=phase.ranges().front().negative;
                for (const auto& range:phase.ranges()) if (range.negative!=negative) return 0;
                return negative?-1:1;
            };
            std::vector<LevelSetMappedSegment> traces;
            std::array<int,5> signs{};
            for (int k=0;k<5;++k) {
                traces.emplace_back(!trajectory && k>0 && k<4?oldSource:controlSource[k],begin,end);
                signs[k]=sign(traces.back());
            }
            const auto& oldTrace=traces[0];const auto& newTrace=traces[4];
            if (!trajectory) for (int k=1;k<4;++k) signs[k]=signs[0]==signs[4]?signs[0]:0;
            const auto temporalValues=[&](long double position) {
                std::array<double,5> values{};
                values[0]=valueVariation<P>(oldTrace,position).value;
                values[4]=valueVariation<P>(newTrace,position).value;
                for (int k=1;k<4;++k) values[k]=trajectory?valueVariation<P>(traces[k],position).value
                    :mix(values[0],values[4],.25*k);
                return values;
            };
            bool uniform=signs[0]!=0;
            for (int k=1;k<5;++k) uniform&=signs[k]==signs[0];
            if (uniform) {
                if (signs[0]<0) {
                    const long double old=oldVelocity.integratedNormalFlux(face.nativeRow,half,0,1);
                    const long double current=currentVelocity.integratedNormalFlux(face.nativeRow,half,0,1);
                    flux+=.5L*dt*(old+current); absolute[face.nativeRow]+=.5L*dt*(std::abs(old)+std::abs(current));
                }
                continue;
            }
            bool proportional=signs[0]!=0;
            if (!trajectory) proportional=proportional && signs[4]!=0 && sameZeros<P>(oldTrace,newTrace);
            else for (int k=1;k<5 && proportional;++k)
                proportional=signs[k]!=0 && sameZeros<P>(oldTrace,traces[k]);
            if (proportional) {
                // Whole-face time crossings need the distributional moving-
                // event derivative in addition to ordinary spatial roots.
                const auto values=temporalValues(.5L);
                const auto roots=Quantum::Geometry::isolateQuarticBernstein(values);
                require(!roots.identicallyZero,"unresolved zero temporal trace");
                std::vector<double> cuts{0};
                for (double root:roots.roots) if (root>0 && root<1) cuts.push_back(root);
                cuts.push_back(1);
                const long double old=oldVelocity.integratedNormalFlux(face.nativeRow,half,0,1);
                const long double current=currentVelocity.integratedNormalFlux(face.nativeRow,half,0,1);
                for (std::size_t k=1;k<cuts.size();++k) {
                    const long double lo=cuts[k-1],hi=cuts[k];
                    const auto middle=bernstein<5>({{values[0],values[1],values[2],values[3],values[4]}},.5L*(lo+hi));
                    require(middle!=0,"unresolved temporal phase interval");
                    if (middle>0) continue;
                    const long double b=.5L*(hi-lo)*(hi+lo),a=hi-lo-b;
                    flux+=dt*(a*old+b*current); absolute[face.nativeRow]+=dt*(a*std::abs(old)+b*std::abs(current));
                }
                for (double theta:roots.roots) for (std::size_t q=0;q<x8.size();++q) {
                    const double position=.5+.5*x8[q];
                    const auto local=temporalValues(position);
                    std::array<long double,4> slope{};
                    for (int k=0;k<4;++k) slope[k]=4*(local[k+1]-local[k]);
                    const long double rate=std::abs(bernstein(slope,theta));
                    require(std::isfinite(rate) && rate>0,"unresolved temporal jump derivative");
                    const long double density=(1-theta)*oldVelocity.normalFluxDensity(face.nativeRow,half,false,position)
                        +theta*currentVelocity.normalFluxDensity(face.nativeRow,half,false,position);
                    const long double endpoint=(theta==0 || theta==1)?.5L:1.L;
                    const long double factor=-endpoint*dt*.5L*w8[q]*theta*density/rate;
                    accumulate(derivative[face.nativeRow],valueVariation<P>(newTrace,position).independent,factor);
                }
                continue;
            }
            dynamic.push_back({face.nativeRow,half,begin,end});
            for (auto point:{begin,end}) {
                std::array<double,5> values{};
                for (int k=0;k<5;++k) values[k]=controlSource[k]->value(point);
                const auto roots=Quantum::Geometry::isolateQuarticBernstein(values);
                if (roots.identicallyZero) {
                    bool pinned=false;
                    for (const auto& sector:oldPatch.cornerSectors)
                        pinned|=same(point,mesh.verts(mesh.vertIndexToLabel[sector.vertexIndex]).coord());
                    require(pinned,"stationary zero endpoint needs an explicit incidence policy");
                    continue; // A declared fixed contact is an endpoint for the entire step.
                }
                for (double event:roots.roots) if (event>0 && event<1) events.push_back(event);
            }
        }
    }
    std::sort(events.begin(),events.end()); events.erase(std::unique(events.begin(),events.end()),events.end());
    std::set<int> queryOwners;
    for (const auto& half:dynamic) {
        const auto face=domain.face(mesh,half.row);
        if (face.negativeCell>=0) queryOwners.insert(face.negativeCell);
        if (face.positiveCell>=0) queryOwners.insert(face.positiveCell);
    }
    const std::vector<int> queryCells(queryOwners.begin(),queryOwners.end());
    const std::array<double,4> x4{{-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526}};
    const std::array<double,4> w4{{.34785484513745385,.65214515486254615,.65214515486254615,.34785484513745385}};
    const auto rule=[&](double lo,double hi,const auto& nodes,const auto& weights,bool sensitivity) {
        Rule sum(faceCount);
        const double middle=.5*lo+.5*hi,half=.5*(hi-lo);
        for (std::size_t q=0;q<nodes.size();++q) {
            require(integral.evaluations<options.maxEvaluations,"time quadrature work limit exhausted");
            const double theta=std::fma(half,nodes[q],middle);
            require(theta>lo && theta<hi,"time quadrature interval has no representable interior sample");
            const auto patch=interpolate(temporal,theta);
            auto source=[&]() {
                if constexpr (quartic<P>) if (patch.cornerSectors.empty())
                    return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
                        LevelSetHermiteGeometry::Snapshot::fromNativeQuarticQueryCells(mesh,patch,queryCells,domain.axisymmetric()));
                return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
                    LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,domain.axisymmetric()));
            }();
            ++integral.evaluations;
            const long double weight=static_cast<long double>(dt)*half*weights[q];
            for (const auto& face:dynamic) {
                const LevelSetMappedSegment phase(source,face.a,face.b);
                const long double old=oldVelocity.integratedNormalFlux(face.row,face.half,phase);
                const long double current=currentVelocity.integratedNormalFlux(face.row,face.half,phase);
                const long double flux=(1-theta)*old+theta*current;
                require(std::isfinite(flux),"nonfinite time-sampled phase flux");
                sum.flux[face.row]+=weight*flux; sum.absolute[face.row]+=weight*std::abs(flux);
                if (!sensitivity) continue;
                for (const auto& root:phaseVariation<P>(phase)) {
                    const long double density=(1-theta)*oldVelocity.normalFluxDensity(face.row,face.half,root.fromEnd,root.offset)
                        +theta*currentVelocity.normalFluxDensity(face.row,face.half,root.fromEnd,root.offset);
                    accumulate(sum.derivative[face.row],root.independent,weight*theta*density);
                }
            }
        }
        return sum;
    };
    struct Interval {
        double lo,hi;int depth;
        Rule fine;
        std::vector<long double> error;
    };
    std::vector<Interval> intervals;
    const auto evaluate=[&](double lo,double hi,int depth) {
        const auto coarse=rule(lo,hi,x4,w4,false);
        auto fine=rule(lo,hi,x8,w8,true);
        std::vector<long double> error(faceCount,0);
        for (int row=0;row<faceCount;++row) error[row]=std::abs(fine.flux[row]-coarse.flux[row]);
        return Interval{lo,hi,depth,std::move(fine),std::move(error)};
    };
    if (!dynamic.empty()) {
        for (std::size_t k=1;k<events.size();++k) intervals.push_back(evaluate(events[k-1],events[k],0));
        // The original tolerance applies to the COMPLETE face time integral.
        // In particular, a tiny interval adjacent to a tangency need not meet
        // an artificial relative accuracy on its vanishing local flux.
        for (;;) {
            std::vector<long double> error(faceCount,0),magnitude=absolute;
            for (const auto& interval:intervals) for (int row=0;row<faceCount;++row) {
                error[row]+=interval.error[row];magnitude[row]+=interval.fine.absolute[row];
            }
            int worstFace=-1;long double worstRatio=1;
            for (const auto& face:dynamic) {
                const int row=face.row;
                const long double tolerance=faceTolerance[row]+options.relativeTolerance*magnitude[row];
                const long double ratio=error[row]/tolerance;
                if (ratio>worstRatio) {worstRatio=ratio;worstFace=row;}
            }
            if (worstFace<0) break;
            std::size_t selected=0;
            for (std::size_t k=1;k<intervals.size();++k)
                if (intervals[k].error[worstFace]>intervals[selected].error[worstFace]) selected=k;
            const auto lo=intervals[selected].lo,hi=intervals[selected].hi;
            const int depth=intervals[selected].depth;
            if (depth>=options.maxDepth) {
                std::ostringstream reason;
                reason<<std::setprecision(17)<<"[InterfaceTemporalPhaseFlux] complete face quadrature did not converge"
                    <<"; face="<<worstFace<<", theta=["<<lo<<','<<hi<<"], error="<<error[worstFace]
                    <<", tolerance="<<faceTolerance[worstFace]+options.relativeTolerance*magnitude[worstFace];
                throw std::runtime_error(reason.str());
            }
            const double middle=.5*lo+.5*hi;
            auto left=evaluate(lo,middle,depth+1),right=evaluate(middle,hi,depth+1);
            intervals[selected]=std::move(left);intervals.push_back(std::move(right));
        }
    }
    integral.intervals=static_cast<int>(intervals.size());
    for (const auto& interval:intervals) for (int row=0;row<faceCount;++row) {
        if (!std::isfinite(budget.positiveFluxByNativeFace[row])) continue;
        budget.positiveFluxByNativeFace[row]+=interval.fine.flux[row];absolute[row]+=interval.fine.absolute[row];
        integral.errorByNativeFace[row]+=interval.error[row];
        for (const auto& column:interval.fine.derivative[row]) derivative[row][column.first]+=column.second;
    }
    // Conservative operation-count estimate for the weighted sample sums and
    // shared-face accumulation; it does not certify root/geometry error.
    const long double neps=(8.L*integral.evaluations+32)*std::numeric_limits<double>::epsilon();
    require(neps<1,"time integration roundoff estimate is unbounded");
    const long double gamma=neps/(1-neps);
    for (const auto& face:domain.faces(mesh)) {
        if (face.negativeCell<0 && face.positiveCell<0) continue;
        const auto row=face.nativeRow; const auto flux=budget.positiveFluxByNativeFace[row];
        require(std::isfinite(flux),"nonfinite integrated face flux");
        finish(derivative[row],result.derivative.positiveFluxByNativeFace[row]);
        for (int side=0;side<2;++side) {
            const int cell=side?face.positiveCell:face.negativeCell;
            if (cell<0) continue;
            budget.outwardFluxByCell[cell]+=(side?-1:1)*flux;
            integral.absoluteTransportByCell[cell]+=absolute[row];
            integral.integrationUncertaintyByCell[cell]+=integral.errorByNativeFace[row]+gamma*absolute[row];
            for (const auto& column:derivative[row]) cellDerivative[cell][column.first]+=(side?-1:1)*column.second;
        }
        budget.boundaryOutward+=face.outwardSign*flux; budget.absoluteFaceFlux+=std::abs(flux); ++budget.activeFaces;
    }
    for (int c:domain.cellIndexByRow()) {
        require(std::isfinite(budget.outwardFluxByCell[c]) && std::isfinite(integral.absoluteTransportByCell[c])
            && std::isfinite(integral.integrationUncertaintyByCell[c]) && integral.integrationUncertaintyByCell[c]>=0,
            "nonfinite cell-integrated flux or numerical uncertainty");
        budget.summedCellOutward+=budget.outwardFluxByCell[c];
        finish(cellDerivative[c],result.derivative.outwardFluxByCell[c]);
    }
    budget.assemblyClosure=budget.summedCellOutward-budget.boundaryOutward;
    return result;
}
} // namespace
Result integrate(const QuadtreeMesh& mesh,const InterfaceTransportDomain& domain,
    const Patch& oldPatch,const Patch& currentPatch,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& oldSource,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& currentSource,
    const InterfaceNormalFluxView& oldVelocity,const InterfaceNormalFluxView& currentVelocity,
    double dt,const Options& options,const Trajectory* trajectory)
{ return integrateImpl(mesh,domain,oldPatch,currentPatch,oldSource,currentSource,oldVelocity,currentVelocity,dt,options,trajectory); }
CubicResult integrate(const QuadtreeMesh& mesh,const InterfaceTransportDomain& domain,
    const Quantum::Geometry::ActiveHermitePatch& oldPatch,const Quantum::Geometry::ActiveHermitePatch& currentPatch,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& oldSource,
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& currentSource,
    const InterfaceNormalFluxView& oldVelocity,const InterfaceNormalFluxView& currentVelocity,
    double dt,const Options& options,const CubicTrajectory* trajectory)
{ return integrateImpl(mesh,domain,oldPatch,currentPatch,oldSource,currentSource,oldVelocity,currentVelocity,dt,options,trajectory); }
}
