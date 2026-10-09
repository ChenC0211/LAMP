#include "Modules/Interface/Compatibility/InterfaceColumnMode.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include <cmath>
#include <stdexcept>

InterfaceColumnMode InterfaceColumnMode::evaluate(const InterfaceState& state,double lo,double hi,double radius,int halfWaves)
{
    if(!std::isfinite(lo) || !std::isfinite(hi) || !(hi>lo) || !std::isfinite(radius) || !(radius>0) || (halfWaves!=1 && halfWaves!=2))
        throw std::invalid_argument("column mode requires finite axial bounds and positive reference radius");
    InterfaceColumnMode result;result.source=state.geometry();
    const double length=hi-lo,k=halfWaves*std::acos(-1.)/length;
    const auto& surfaces=state.topology().surfaces();
    if(surfaces.size()!=1)return result;
    const auto& surface=surfaces.front();const auto& chain=surface.geometry.nativeChain();
    if(!chain.axisymmetric() || chain.closed())return result;
    const auto* first=chain.endpointSymmetryPlane(false);const auto* last=chain.endpointSymmetryPlane(true);
    if(!first || !last || first->normalAxis!=1 || last->normalAxis!=1
        || !((first->coordinate==lo && last->coordinate==hi) || (first->coordinate==hi && last->coordinate==lo)))return result;
    bool graph=true;long double mean=0,mode=0,coverage=0;
    for(const auto& g:chain.graphs())graph&=!g.coordinateIsX;
    for(const auto& panel:surface.geometry.panels())for(const auto& q:panel.quadratureNodes) {
        graph&=q.frame.normal.x>0;
        const long double dz=q.weight*q.frame.normal.x;
        coverage+=dz;mean+=dz*q.frame.point.x;
        mode+=dz*q.frame.point.x*std::cos(k*(q.frame.point.y-lo));
    }
    result.axialLength=static_cast<double>(coverage);
    result.meanRadius=static_cast<double>(mean/length);
    result.cosineAmplitude=static_cast<double>(2*mode/length);
    result.relativeAmplitude=static_cast<double>(2*mode/(length*radius));
    result.singleRadiusGraph=graph && std::abs(coverage-length)<1e-10*length
        && std::isfinite(result.relativeAmplitude) && std::isfinite(result.meanRadius);
    return result;
}
