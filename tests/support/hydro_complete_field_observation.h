#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <stdexcept>

// Independent physical observations for complete-field flow controls.
// No production mass matrix, affine projection or minimum-energy rebuild.
namespace HydroCompleteObservation {
using Velocity=Quantum::Hydrodynamics::HydroCompatibleMomentSpline::NativeVelocity;
template<class F>void rectangle(Vector2D lo,Vector2D hi,F&& f){
    if(!(hi.x>lo.x&&hi.y>lo.y))return;
    static constexpr double x[6]={.12523340851146891547,.36783149899818019375,.58731795428661744730,
        .76990267419430468704,.90411725637047485668,.98156063424671925069};
    static constexpr double w[6]={.24914704581340278500,.23349253653835480876,.20316742672306592175,
        .16007832854334622633,.10693932599531843096,.047175336386511827195};
    const auto h=(hi-lo)*.5,center=(hi+lo)*.5;
    for(int i=0;i<6;i++)for(int j=0;j<6;j++)for(int sx:{-1,1})for(int sy:{-1,1})
        f(Vector2D{center.x+sx*h.x*x[i],center.y+sy*h.y*x[j]},h.x*h.y*w[i]*w[j]);
}
template<class F>void volume(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,F&& f){
    // Native Q8 fields are polynomial on a leaf; radial division off-axis
    // makes RZ rational. Four radial panels resolve that smooth factor.
    const bool rz=domain.axisymmetric();const int nx=rz?4:2;
    for(int index:domain.cellIndexByRow()){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);const auto h=cell.size(),lo=cell.centerCoord()-h*.5;
        for(int i=0;i<nx;i++)for(int j=0;j<2;j++)
            rectangle({lo.x+h.x*i/nx,lo.y+h.y*j/2},{lo.x+h.x*(i+1)/nx,lo.y+h.y*(j+1)/2},
                [&](Vector2D p,double da){f(p,da*(rz?2*std::acos(-1.)*p.x:1));});
    }
}
template<class Reference>double errorSquared(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const Velocity& field,Reference&& reference){
    long double sum=0;
    volume(mesh,domain,[&](Vector2D p,double dv){const auto a=field.firstJet(p).value,b=reference(p);
        const long double x=a.x-b.x,y=a.y-b.y;sum+=dv*(x*x+y*y);});
    if(!(sum>=0&&std::isfinite(sum)))throw std::runtime_error("nonfinite complete-field observation");
    return double(sum);
}
inline double distanceSquared(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const Velocity& a,const Velocity& b){return errorSquared(mesh,domain,a,[&](Vector2D p){return b.firstJet(p).value;});}
inline double normSquared(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Velocity& a){
    return errorSquared(mesh,domain,a,[](Vector2D){return Vector2D{0,0};});
}
}
