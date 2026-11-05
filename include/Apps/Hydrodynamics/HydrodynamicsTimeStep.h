#pragma once
#include "RunControls.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

// Passive App policy result. Infinite bounds are unconstrained/non-finite;
// their *_bound_available is false. A measured zero flow_rate remains available.
struct HydrodynamicsTimeStepReport
{
    double proposed_dt=0,actual_dt=0,end=0;
    long double flow_rate=0,capillary_rate=0;
    double interface_minimum_spacing=0;
    long double flow_bound=std::numeric_limits<long double>::infinity();
    long double capillary_bound=std::numeric_limits<long double>::infinity();
    long double growth_bound=0,max_bound=0;
    bool flow_rate_available=true,capillary_rate_available=false;
    bool flow_bound_available=false,capillary_bound_available=false;
    bool growth_bound_available=false,max_bound_available=true;
    bool interface_spacing_available=false,capillary_available=false;
    bool endpoint_aligned=false,minimum_clock_rounded=false,final=false,below_min_final=false;
    const char* limiter="maxDt";
};

// No mesh, Module, solver, output or state owner is retained here. Rates must
// come from the actual accepted-state native stability query, never base-grid h.
inline HydrodynamicsTimeStepReport selectHydrodynamicsTimeStep(
    const RunControls& run,double begin,double previousActualDt,int acceptedSteps,
    long double flowRate,bool capillaryAvailable,long double capillaryRate,
    double interfaceMinimumSpacing,double maximumCapillaryCourant)
{
    const auto require=[](bool ok,const char* why){
        if(!ok)throw std::invalid_argument(std::string("[HydrodynamicsTimeStep] ")+why);
    };
    require(run.adaptive.has_value(),"adaptive policy is absent");
    const auto& policy=*run.adaptive;
    require(std::isfinite(begin)&&std::isfinite(run.endTime)&&run.endTime>begin&&acceptedSteps>=0,
        "invalid accepted clock");
    require(std::isfinite(policy.minDt)&&policy.minDt>0&&std::isfinite(policy.maxDt)
        &&policy.maxDt>=policy.minDt&&std::isfinite(run.dt)&&run.dt>=policy.minDt&&run.dt<=policy.maxDt,
        "invalid adaptive min/seed/max interval");
    require(std::isfinite(policy.targetCourant)&&policy.targetCourant>0
        &&std::isfinite(policy.safety)&&policy.safety>0&&policy.safety<=1
        &&std::isfinite(policy.maxGrowth)&&policy.maxGrowth>=1,"invalid adaptive Courant/safety/growth");
    require(std::isfinite(flowRate)&&flowRate>=0,"invalid available native flow rate");
    require(std::isfinite(maximumCapillaryCourant)&&maximumCapillaryCourant>0,
        "invalid existing capillary Courant limit");
    require(std::isfinite(capillaryRate)&&std::isfinite(interfaceMinimumSpacing)
        &&(capillaryAvailable?(capillaryRate>0&&interfaceMinimumSpacing>0):
            (capillaryRate==0&&interfaceMinimumSpacing==0)),"invalid capillary availability/rate/spacing");
    if(acceptedSteps>0)require(std::isfinite(previousActualDt)&&previousActualDt>0,"invalid previous accepted dt");
    if(run.stepRetry)require(run.stepRetry->recoverAfter==0&&run.stepRetry->minDt==policy.minDt,
        "adaptive retry must share minDt and leave growth to the adaptive policy");

    HydrodynamicsTimeStepReport report;
    report.flow_rate=flowRate;report.capillary_rate=capillaryRate;
    report.interface_minimum_spacing=interfaceMinimumSpacing;
    report.capillary_available=report.capillary_rate_available=report.interface_spacing_available=capillaryAvailable;
    report.max_bound=policy.maxDt;
    report.growth_bound=acceptedSteps==0?static_cast<long double>(run.dt):
        static_cast<long double>(previousActualDt)*policy.maxGrowth;
    report.growth_bound_available=std::isfinite(report.growth_bound);
    if(flowRate>0){
        report.flow_bound=static_cast<long double>(policy.safety)*policy.targetCourant/flowRate;
        require(report.flow_bound>0,"positive native flow bound is not representable");
        report.flow_bound_available=std::isfinite(report.flow_bound);
    }
    if(capillaryAvailable){
        report.capillary_bound=static_cast<long double>(policy.safety)*maximumCapillaryCourant/capillaryRate;
        require(report.capillary_bound>0,"positive capillary bound is not representable");
        report.capillary_bound_available=std::isfinite(report.capillary_bound);
    }
    long double bound=report.max_bound;
    if(report.growth_bound<bound){bound=report.growth_bound;report.limiter=acceptedSteps==0?"initialDt":"growth";}
    if(report.flow_bound<bound){bound=report.flow_bound;report.limiter="flow";}
    if(report.capillary_bound<bound){bound=report.capillary_bound;report.limiter="capillary";}
    const auto stability=std::min(report.flow_bound,report.capillary_bound);
    const double remaining=run.endTime-begin;
    require(std::isfinite(remaining)&&remaining>0,"invalid remaining physical interval");
    const long double rounding=8.L*std::numeric_limits<double>::epsilon()
        *(static_cast<long double>(acceptedSteps)+1)*std::max(std::abs(begin),std::abs(run.endTime));
    const bool smallFinal=remaining<policy.minDt;
    if(smallFinal){
        // A genuinely final short interval is allowed only when already stable.
        // Never clamp an unstable interval upward to minDt or split a tiny tail.
        require(static_cast<long double>(remaining)<=stability,"remaining final interval is unstable below minDt");
        const auto soft=std::min(report.max_bound,report.growth_bound);
        if(static_cast<long double>(remaining)>soft){
            require(static_cast<long double>(remaining)-soft<=rounding,
                "remaining final interval exceeds nominal max/growth beyond clock rounding");
            report.endpoint_aligned=true;
        }
        report.proposed_dt=remaining;report.end=run.endTime;report.actual_dt=remaining;
        report.final=report.below_min_final=true;report.limiter="endTime";
        return report;
    }
    require(std::isfinite(bound)&&bound>=policy.minDt,"stability/growth bound is below adaptive minDt");
    const auto selected=std::min(bound,static_cast<long double>(remaining));
    double proposed=static_cast<double>(selected);
    if(static_cast<long double>(proposed)>selected)proposed=std::nextafter(proposed,0.);
    require(std::isfinite(proposed)&&proposed>=policy.minDt,"selected adaptive dt is below minDt or unrepresentable");
    if(static_cast<long double>(remaining)<bound)report.limiter="endTime";
    report.proposed_dt=proposed;
    const auto unaligned=std::min(run.endTime,begin+proposed);
    double end=run.stepEnd(begin,proposed,acceptedSteps),actual=end-begin;
    bool aligned=end==run.endTime&&(end>unaligned||static_cast<long double>(actual)>bound);
    if(aligned&&static_cast<long double>(actual)>bound)
        require(static_cast<long double>(actual)-bound<=rounding,"endpoint extension exceeds the existing clock rounding budget");
    // Clock alignment may extend only nominal max/growth/seed bounds, never
    // accepted-state flow/capillary stability. Ordinary endpoint rounding also
    // cannot increase any active bound. One predecessor removes that round-up.
    if(static_cast<long double>(actual)>stability){
        end=unaligned;actual=end-begin;aligned=false;
    }
    if((!aligned&&static_cast<long double>(actual)>bound)||static_cast<long double>(actual)>stability){
        end=std::nextafter(end,begin);actual=end-begin;aligned=false;
    }
    require(end>begin&&end<=run.endTime&&std::isfinite(actual)&&actual>0,"adaptive endpoint is not representable");
    require(static_cast<long double>(actual)<=stability,"actual adaptive interval exceeds flow/capillary stability");
    require(aligned||static_cast<long double>(actual)<=bound,"actual adaptive interval exceeds a nominal bound");
    if(actual<policy.minDt&&!(end==run.endTime&&remaining<policy.minDt)){
        // The nominal policy bound is still >= minDt. Subtraction at a
        // representable clock (and a predecessor used to respect an upper
        // bound) may make its ACTUAL interval a few local clock ulps shorter.
        // This is not permission to choose a genuinely below-min proposal.
        const long double localRounding=4.L*std::numeric_limits<double>::epsilon()
            *std::max(std::abs(begin),std::abs(end));
        require(proposed>=policy.minDt&&static_cast<long double>(policy.minDt)-actual
            <=std::min(rounding,localRounding),"actual adaptive interval is genuinely below minDt");
        report.minimum_clock_rounded=true;
    }
    report.end=end;report.actual_dt=actual;report.endpoint_aligned=aligned;report.final=end==run.endTime;
    report.below_min_final=report.final&&actual<policy.minDt;
    return report;
}
