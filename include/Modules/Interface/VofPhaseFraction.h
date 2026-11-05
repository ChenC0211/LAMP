#pragma once
#include <array>

namespace LAMP::Interface
{
// ONE primary mathematical fraction: alpha=m or alpha=1-m. The minor phase
// is retained directly; a rounded dominant alpha is never its authority.
class VofPhaseFraction
{
public:
    enum class Phase { Liquid, Gas };
    VofPhaseFraction()=default; // Pure gas.
    static VofPhaseFraction fromAlpha(double); // Exact import of that double.
    static VofPhaseFraction fromMinor(double,Phase);
    static VofPhaseFraction fromInventories(long double liquid,long double gas,
        long double capacity,long double relativeTolerance=1e-12L);
    double minorFraction()const noexcept{return minor_;}
    Phase minorPhase()const noexcept{return phase_;}
    bool isMixed()const noexcept{return minor_>0;}
    bool isPureLiquid()const noexcept{return minor_==0&&phase_==Phase::Gas;}
    bool isPureGas()const noexcept{return minor_==0&&phase_==Phase::Liquid;}
    // Derived display/legacy view only. It can equal 1 for a real mixed cell.
    double alphaRounded()const noexcept;
    // Exact two-native-double expansion of the primary alpha, not its view.
    std::array<double,2> alphaComponents()const;
    bool hasSingleDoubleAlpha()const;
    // Independently stable physical phase volumes: liquid, gas. The dominant
    // complement may round; the directly retained minor can never disappear.
    std::array<long double,2> phaseVolumes(long double capacity)const;
    // Same primary P0 weighted nonnegative material and stable alpha difference.
    long double mixture(long double liquidValue,long double gasValue)const;
    long double difference(const VofPhaseFraction& other)const;
    friend bool operator==(VofPhaseFraction a,VofPhaseFraction b)noexcept
    {return a.minor_==b.minor_&&a.phase_==b.phase_;}
    friend bool operator!=(VofPhaseFraction a,VofPhaseFraction b)noexcept{return !(a==b);}
private:
    VofPhaseFraction(double value,Phase phase):minor_(value),phase_(phase){}
    double minor_=0;
    Phase phase_=Phase::Liquid;
};
}
