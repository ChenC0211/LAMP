#include "Modules/Interface/VofPhaseFraction.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace LAMP::Interface
{
namespace
{
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[VOF PhaseFraction] ")+reason);}
void environment(){Geometry::Detail::ExactExpansion::validateActualEnvironment();}
}
VofPhaseFraction VofPhaseFraction::fromAlpha(double alpha)
{
    environment();require(std::isfinite(alpha)&&alpha>=0&&alpha<=1,"alpha is outside [0,1]; clipping prohibited");
    // Sterbenz: for alpha>=.5 the subtraction retains its exact double meaning.
    return alpha<=.5?fromMinor(alpha,Phase::Liquid):fromMinor(1-alpha,Phase::Gas);
}
VofPhaseFraction VofPhaseFraction::fromMinor(double minor,Phase phase)
{
    environment();require(std::isfinite(minor)&&minor>=0&&minor<=.5,"minor fraction is outside [0,.5]; clipping prohibited");
    require(phase==Phase::Liquid||phase==Phase::Gas,"invalid minor phase tag");
    if(minor==.5)phase=Phase::Liquid; // Unique representation of the same value.
    return VofPhaseFraction(minor==0?0:minor,phase);
}
VofPhaseFraction VofPhaseFraction::fromInventories(long double liquid,long double gas,
    long double capacity,long double tolerance)
{
    environment();require(std::isfinite(capacity)&&capacity>0&&std::isfinite(liquid)&&liquid>=0&&std::isfinite(gas)&&gas>=0,
        "invalid independent physical phase inventories");
    require(std::isfinite(tolerance)&&tolerance>0&&tolerance<=1e-12L,"invalid own-phase tolerance");
    require(liquid>0||gas>0,"positive capacity has two empty phase inventories");
    const bool liquidMinor=liquid<=gas;const auto retained=liquidMinor?liquid:gas;
    const auto minor=static_cast<double>(retained/capacity);
    require(retained==0||minor>0,"positive minor fraction is not representable");
    const auto result=fromMinor(minor,liquidMinor?Phase::Liquid:Phase::Gas);const auto actual=result.phaseVolumes(capacity);
    const auto error=[&](long double value,long double expected){return expected>0?std::abs(value-expected)/expected:(value==0?0.L:1.L);};
    const auto residual=std::max(error(actual[0],liquid),error(actual[1],gas));
    if(residual>tolerance){std::ostringstream why;why.precision(18);why<<"[VOF PhaseFraction] independent inventories do not admit this primary fraction"
        <<std::scientific<<" liquid="<<liquid<<" gas="<<gas<<" capacity="<<capacity<<" minor="<<minor<<" ownResidual="<<residual;
        throw std::invalid_argument(why.str());}
    return result;
}
double VofPhaseFraction::alphaRounded()const noexcept
{return phase_==Phase::Liquid?minor_:1-minor_;}
std::array<double,2> VofPhaseFraction::alphaComponents()const
{
    environment();if(phase_==Phase::Liquid)return {{minor_,0}};
    const double high=1-minor_,difference=1-high,low=difference-minor_;
    // For m in [0,.5], Sterbenz plus the retained small correction yields
    // (high+low)=1-m exactly, even when high rounds all the way to one.
    return {{high,low}};
}
bool VofPhaseFraction::hasSingleDoubleAlpha()const{return alphaComponents()[1]==0;}
std::array<long double,2> VofPhaseFraction::phaseVolumes(long double capacity)const
{
    environment();require(std::isfinite(capacity)&&capacity>0,"invalid physical phase capacity");
    const auto minor=static_cast<long double>(minor_)*capacity,dominant=capacity-minor;
    require(std::isfinite(minor)&&std::isfinite(dominant)&&minor>=0&&dominant>0,"phase volume overflow or invalid complement");
    require(minor_==0||minor>0,"positive physical minor volume underflow; phase loss prohibited");
    return phase_==Phase::Liquid?std::array<long double,2>{{minor,dominant}}:std::array<long double,2>{{dominant,minor}};
}
long double VofPhaseFraction::mixture(long double liquid,long double gas)const
{
    environment();require(std::isfinite(liquid)&&liquid>=0&&std::isfinite(gas)&&gas>=0,"invalid nonnegative phase material");
    const auto minor=phase_==Phase::Liquid?liquid:gas,dominant=phase_==Phase::Liquid?gas:liquid;
    const auto value=std::fma(static_cast<long double>(minor_),minor-dominant,dominant);
    require(std::isfinite(value)&&value>=0,"nonfinite physical P0 mixture");
    require(value>0||(dominant==0&&(minor_==0||minor==0)),"positive physical P0 mixture underflow");return value;
}
long double VofPhaseFraction::difference(const VofPhaseFraction& other)const
{
    environment();const long double a=minor_,b=other.minor_;
    if(phase_==other.phase_)return phase_==Phase::Liquid?a-b:b-a;
    const auto combined=(a-.5L)+(b-.5L);return phase_==Phase::Liquid?combined:-combined;
}
}
