#include "Modules/Electric/ElectricVofMaxwellStress.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using Clock=LAMP::TimeLayer;
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[ElectricVofMaxwellStress] ")+why);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
double rounded(const FV::Dyadic& value)
{const auto result=value.roundedDouble();require(std::isfinite(result),"nonfinite derived Maxwell scalar");return result;}
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
double rounded(const FV::NativeSmallExact::Product& value)
{const auto result=value.roundedDouble();require(std::isfinite(result),"nonfinite derived Maxwell scalar");return result;}
double rounded(const FV::NativeSmallExact::Sum& value)
{const auto result=value.roundedDouble();require(std::isfinite(result),"nonfinite derived Maxwell scalar");return result;}
#endif
double interpolate(long double a,long double b,long double t)
{
    require(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(t)&&t>=0&&t<=1,"invalid RT0 interpolation");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if constexpr(FV::NativeSmallExact::supported){
        const auto first=FV::NativeSmallExact::product({1-t,a});
        const auto second=FV::NativeSmallExact::product({t,b});
        return rounded(FV::NativeSmallExact::combine(first,1,second,1));
    }
#endif
    auto result=FV::Dyadic::product({1-t,a});result.add(FV::Dyadic::product({t,b}));return rounded(result);
}
FV::Dyadic representedDisplacementDivergence(const ElectricVofMaxwellCell& c,bool rz)
{
    // The same represented whole-side densities define D=(a+bx,c+dz)
    // in XY and D=(a/r+br,c+dz) in RZ. Native reciprocals are retained
    // finite coefficients; all following products/sums stay exact.
    const auto side=[&](int i){const auto value=c.sideDisplacement[i]/c.sidePhysicalMeasure[i];
        require(std::isfinite(value),"unrepresentable whole side displacement density");return value;};
    const long double x0=c.lower.x,x1=c.upper.x,z0=c.lower.y,z1=c.upper.y;
    const long double inverseWidth=1/(x1-x0),inverseHeight=1/(z1-z0);
    require(std::isfinite(inverseWidth)&&std::isfinite(inverseHeight),"unrepresentable RT0 divergence width reciprocal");
    auto value=FV::Dyadic::product({side(3),inverseHeight});
    value.add(FV::Dyadic::product({side(2),inverseHeight}),-1);
    if(!rz){
        value.add(FV::Dyadic::product({side(1),inverseWidth}));
        value.add(FV::Dyadic::product({side(0),inverseWidth}),-1);
    }else if(x0==0){
        // Regular axis reconstruction is Dr=b*r, with div(D)=2*b+d.
        value.add(FV::Dyadic::product({2.L,side(1),inverseWidth}));
    }else{
        const long double inverseRadiusSum=1/(x0+x1);
        require(std::isfinite(inverseRadiusSum),"unrepresentable RT0 divergence radius reciprocal");
        value.add(FV::Dyadic::product({2.L,x1,side(1),inverseWidth,inverseRadiusSum}));
        value.add(FV::Dyadic::product({2.L,x0,side(0),inverseWidth,inverseRadiusSum}),-1);
    }
    return value;
}
}
struct ElectricVofMaxwellStress::Data
{
    ElectricVofField source;
    bool rz=false;
    std::vector<ElectricVofMaxwellCell> cells;
    std::vector<FV::Dyadic> representedDivergence;
    ElectricVofMaxwellReport report;
    explicit Data(const ElectricVofField& field):source(field),rz(field.geometry().domain().axisymmetric()){}
};
ElectricVofMaxwellStress::ElectricVofMaxwellStress(std::shared_ptr<const Data> data):data_(std::move(data)){}
ElectricVofMaxwellStress ElectricVofMaxwellStress::prepare(const QuadtreeMesh& mesh,const ElectricVofField& field)
{
    field.validateFor(mesh,field.geometry(),field.sourceCharge(),field.candidateTime());
    require(sameClock(field.geometryTime(),field.candidateTime()),"stress requires coeval full geometry/candidate clock");
    auto data=std::make_shared<Data>(field);const auto& domain=field.geometry().domain();
    const int count=mesh.cellIndexNum();require(field.epsilonByCell().size()==static_cast<std::size_t>(count),"foreign native epsilon size");
    data->cells.resize(count);data->report.cells=domain.cellIndexByRow().size();
    std::vector<std::array<FV::Dyadic,4>> displacement(count),measure(count);
    for(int ci=0;ci<count;++ci){auto& c=data->cells[ci];c.nativeCell=ci;
        if(domain.cellRowByIndex()[ci]<0)continue;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[ci]);c.lower=mesh.verts(cell.cellToVert[0]).coord();c.upper=mesh.verts(cell.cellToVert[2]).coord();
        const long double dx=static_cast<long double>(c.upper.x)-c.lower.x,dy=static_cast<long double>(c.upper.y)-c.lower.y;
        require(dx>0&&dy>0&&std::isfinite(dx)&&std::isfinite(dy)&&(!data->rz||c.lower.x>=0),"invalid actual RT0 rectangle");
        c.epsilon=field.epsilonByCell()[ci];require(std::isfinite(c.epsilon)&&c.epsilon>0,"nonpositive actual cell epsilon");
        if(data->rz){c.sidePhysicalMeasure[0]=CartesianMeasure::twoPi*static_cast<long double>(c.lower.x)*dy;
            c.sidePhysicalMeasure[1]=CartesianMeasure::twoPi*static_cast<long double>(c.upper.x)*dy;
            c.sidePhysicalMeasure[2]=c.sidePhysicalMeasure[3]=CartesianMeasure::twoPi*dx*(static_cast<long double>(c.lower.x)+c.upper.x)/2;
            if(c.lower.x==0)++data->report.axisCells;
        }else c.sidePhysicalMeasure={dy,dy,dx,dx};
        for(int side=0;side<4;++side)require(std::isfinite(c.sidePhysicalMeasure[side])&&c.sidePhysicalMeasure[side]>=0,"unrepresentable true whole side measure");
    }
    require(field.faces().size()==static_cast<std::size_t>(domain.nativeFaceCount()),"foreign actual native D faces");
    for(int row=0;row<domain.nativeFaceCount();++row){const auto f=domain.face(mesh,row);const auto& actual=field.faces()[row];
        if(f.negativeCell<0&&f.positiveCell<0){
            require(actual.nativeRow==row&&actual.negativeCell<0&&actual.positiveCell<0&&actual.physicalMeasure==0&&
                actual.epsilon==0&&actual.conductivity==0&&actual.implicitCoefficient==0&&actual.gradientCoefficients.empty()&&
                actual.gradientConstant==0&&actual.integratedField==0&&actual.displacement==0&&actual.current==0,
                "both-excluded native Maxwell face must retain the actual zero field sentinel");
            continue;
        }
        require(actual.nativeRow==row&&actual.negativeCell==f.negativeCell&&actual.positiveCell==f.positiveCell&&
            actual.physicalMeasure==f.area&&std::isfinite(actual.displacement),"different actual native D incidence/measure");
        for(int ci:{f.negativeCell,f.positiveCell})if(ci>=0){const auto& c=data->cells[ci];int side=-1;
            if(f.horizontal){if(f.center.y==c.lower.y)side=2;else if(f.center.y==c.upper.y)side=3;}
            else{if(f.center.x==c.lower.x)side=0;else if(f.center.x==c.upper.x)side=1;}
            require(side>=0,"actual native D face does not lie on its cell side");
            displacement[ci][side].add(actual.displacement);measure[ci][side].add(actual.physicalMeasure);
        }
    }
    for(int ci=0;ci<count;++ci)if(domain.cellRowByIndex()[ci]>=0)for(int side=0;side<4;++side){auto& c=data->cells[ci];
        const auto area=c.sidePhysicalMeasure[side],actual=measure[ci][side].rounded();
        c.sideDisplacement[side]=displacement[ci][side].rounded();
        require(std::isfinite(actual)&&std::isfinite(c.sideDisplacement[side]),"unrepresentable actual whole side D sum");
        if(area==0)require(data->rz&&side==0&&c.lower.x==0&&actual==0&&displacement[ci][side].compare(FV::Dyadic{})==0,"true axis requires zero actual D/area");
        else{const auto defect=std::abs(actual-area)/area;require(defect<=1e-12L,"split-face measures do not cover true whole side");
            data->report.maximumRelativeSideMeasureDefect=std::max(data->report.maximumRelativeSideMeasureDefect,static_cast<double>(defect));}
    }
    data->representedDivergence.resize(count);
    for(int ci=0;ci<count;++ci)if(domain.cellRowByIndex()[ci]>=0)
        data->representedDivergence[ci]=representedDisplacementDivergence(data->cells[ci],data->rz);
    return ElectricVofMaxwellStress(std::move(data));
}
const ElectricVofField& ElectricVofMaxwellStress::field() const
{require(bool(data_),"consumed stress owner");return data_->source;}
Clock ElectricVofMaxwellStress::timeLayer() const{return field().candidateTime();}
bool ElectricVofMaxwellStress::axisymmetric() const{require(bool(data_),"consumed stress owner");return data_->rz;}
const std::vector<ElectricVofMaxwellCell>& ElectricVofMaxwellStress::cells() const
{require(bool(data_),"consumed stress owner");return data_->cells;}
const ElectricVofMaxwellReport& ElectricVofMaxwellStress::report() const
{require(bool(data_),"consumed stress owner");return data_->report;}
void ElectricVofMaxwellStress::validateFor(const QuadtreeMesh& mesh,const ElectricVofField& source,Clock clock) const
{
    require(bool(data_),"consumed stress owner");clock.time();require(sameClock(clock,timeLayer()),"different complete Maxwell clock");
    data_->source.validateFor(mesh,data_->source.geometry(),data_->source.sourceCharge(),clock);
    // The first original validation has just checked this actual live mesh,
    // every retained owner and the full clock. Exact immutable Data equality
    // makes the second call identical, with no callback/mutation in between.
    // A different Field, including equal values, keeps its original full check.
    if(source.data_!=data_->source.data_)
        source.validateFor(mesh,source.geometry(),source.sourceCharge(),clock);
    // Field has no public sameInstance port. Its readonly potential vector is
    // owned by the immutable Field Data: copies share it, fresh solves do not.
    require(&source.potential()==&data_->source.potential(),"foreign same-value electric field owner");
    require(sameClock(source.geometryTime(),clock),"foreign noncoeval electric geometry");
}
ElectricVofMaxwellSample ElectricVofMaxwellStress::sample(int nativeCell,Vector2D point) const
{
    require(bool(data_),"consumed stress owner");require(nativeCell>=0&&static_cast<std::size_t>(nativeCell)<data_->cells.size(),"native stress cell out of range");
    require(field().geometry().domain().cellRowByIndex()[nativeCell]>=0,"inactive native cell has no Maxwell reconstruction");
    const auto& c=data_->cells[nativeCell];require(std::isfinite(point.x)&&std::isfinite(point.y)&&
        point.x>=c.lower.x&&point.x<=c.upper.x&&point.y>=c.lower.y&&point.y<=c.upper.y,"stress sample outside actual closed cell rectangle");
    const long double x0=c.lower.x,x1=c.upper.x,z0=c.lower.y,z1=c.upper.y;
    const auto side=[&](int i){const auto value=c.sideDisplacement[i]/c.sidePhysicalMeasure[i];
        require(std::isfinite(value),"unrepresentable whole side displacement density");return value;};
    ElectricVofMaxwellSample result;result.nativeCell=nativeCell;result.point=point;
    if(!data_->rz)result.displacement.x=interpolate(side(0),side(1),(static_cast<long double>(point.x)-x0)/(x1-x0));
    else if(x0==0){
        // r*Dr=b*r^2, so Dr=b*r including the actual r=0 limit.
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
        if constexpr(FV::NativeSmallExact::supported)
            result.displacement.x=point.x==0?0:rounded(FV::NativeSmallExact::product({side(1),static_cast<long double>(point.x)/x1}));
        else
#endif
        result.displacement.x=point.x==0?0:rounded(FV::Dyadic::product({side(1),static_cast<long double>(point.x)/x1}));
    }else if(point.x==c.lower.x)result.displacement.x=rounded(FV::Dyadic(side(0)));
    else if(point.x==c.upper.x)result.displacement.x=rounded(FV::Dyadic(side(1)));
    else{
        const long double r=point.x,t=((r-x0)/(x1-x0))*((r+x0)/(x1+x0));
        require(std::isfinite(t)&&t>=0&&t<=1,"invalid RZ squared-radius interpolation");
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
        if constexpr(FV::NativeSmallExact::supported){
            const auto first=FV::NativeSmallExact::product({1-t,x0/r,side(0)});
            const auto second=FV::NativeSmallExact::product({t,x1/r,side(1)});
            result.displacement.x=rounded(FV::NativeSmallExact::combine(first,1,second,1));
        }else
#endif
        {
        auto value=FV::Dyadic::product({1-t,x0/r,side(0)});value.add(FV::Dyadic::product({t,x1/r,side(1)}));
        result.displacement.x=rounded(value);
        }
    }
    result.displacement.y=interpolate(side(2),side(3),(static_cast<long double>(point.y)-z0)/(z1-z0));
    result.electricField={result.displacement.x/c.epsilon,result.displacement.y/c.epsilon};
    require(std::isfinite(result.electricField.x)&&std::isfinite(result.electricField.y),"unrepresentable actual cell electric field");
    // Complete contractions use the represented epsilon/E. Final RNE zero is
    // legal for a derived tensor, and never deletes a Q or phase inventory.
    const auto ex=result.electricField.x,ey=result.electricField.y;
    // Differentiate the retained cell RT0 stress analytically. Since curl(E)
    // vanishes within the cell, div(T)=div(D)*E. In RZ this identity includes
    // (Trr-Thoop)/r and Trz/r; regular Dr=0 at the axis gives radial limit 0.
    // The represented field and tensor above retain their original rounding;
    // this derivative port evaluates that analytic identity at represented E.
    const auto& divD=data_->representedDivergence[nativeCell];
    result.tensorDivergence={rounded(divD.multiplied(ex)),rounded(divD.multiplied(ey))};
#if defined(_MSC_VER) && defined(_M_X64) && !defined(_M_CEE) && !defined(__FAST_MATH__) && (!defined(_M_FP_FAST) || !_M_FP_FAST)
    if constexpr(FV::NativeSmallExact::supported){
        const auto x=FV::NativeSmallExact::product({.5L,c.epsilon,ex,ex});
        const auto y=FV::NativeSmallExact::product({.5L,c.epsilon,ey,ey});
        const auto xx=FV::NativeSmallExact::combine(x,1,y,-1);
        const auto yy=FV::NativeSmallExact::combine(y,1,x,-1);
        const auto hoop=FV::NativeSmallExact::combine(x,-1,y,-1);
        result.tensor={rounded(xx),rounded(yy),rounded(FV::NativeSmallExact::product({c.epsilon,ex,ey})),rounded(hoop)};
        return result;
    }
#endif
    const auto x=FV::Dyadic::product({.5L,c.epsilon,ex,ex}),y=FV::Dyadic::product({.5L,c.epsilon,ey,ey});
    auto xx=x;xx.add(y,-1);auto yy=y;yy.add(x,-1);auto hoop=FV::Dyadic{};hoop.add(x,-1);hoop.add(y,-1);
    result.tensor={rounded(xx),rounded(yy),rounded(FV::Dyadic::product({c.epsilon,ex,ey})),rounded(hoop)};
    return result;
}
}
