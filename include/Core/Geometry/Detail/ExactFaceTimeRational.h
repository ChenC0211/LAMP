#pragma once
#include "Core/Geometry/Detail/ExactExpansion.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Quantum::Geometry::Detail::FaceTime
{
namespace Exact=ExactExpansion;
// Shared original face-time rational algebra. Only power-of-two rescaling;
// neither the scalar sign nor the clipped event is a rounded approximation.
struct FaceTimeRational
{
    Exact::Expansion n,d;
    explicit FaceTimeRational(double value=0):n(Exact::fromDouble(value)),d(Exact::fromDouble(1)){}
    FaceTimeRational(Exact::Expansion numerator,Exact::Expansion denominator):n(std::move(numerator)),d(std::move(denominator))
    {
        const int sign=Exact::sign(d);
        if(!sign)throw std::runtime_error("[LevelSetHermiteGeometry] zero face-time rational denominator");
        if(sign<0){n=Exact::subtract(Exact::Expansion{},n);d=Exact::subtract(Exact::Expansion{},d);}
        std::array<Exact::Expansion,1> a{{n}};Exact::centerHomogeneous(a,d);n=std::move(a[0]);
    }
};
inline FaceTimeRational operator+(const FaceTimeRational& a,const FaceTimeRational& b)
{return {Exact::add(Exact::multiply(a.n,b.d),Exact::multiply(b.n,a.d)),Exact::multiply(a.d,b.d)};}
inline FaceTimeRational operator-(const FaceTimeRational& a,const FaceTimeRational& b)
{return {Exact::subtract(Exact::multiply(a.n,b.d),Exact::multiply(b.n,a.d)),Exact::multiply(a.d,b.d)};}
inline FaceTimeRational operator*(const FaceTimeRational& a,const FaceTimeRational& b)
{return {Exact::multiply(a.n,b.n),Exact::multiply(a.d,b.d)};}
inline FaceTimeRational operator/(const FaceTimeRational& a,const FaceTimeRational& b)
{return {Exact::multiply(a.n,b.d),Exact::multiply(a.d,b.n)};}
inline long double faceTimeValue(const FaceTimeRational& x)
{
    long double n=0,d=0;
    for(double c:x.n.components())n+=c;
    for(double c:x.d.components())d+=c;
    const long double value=n/d;
    if(!std::isfinite(value) || (value==0 && Exact::sign(x.n)!=0))
        throw std::runtime_error("[LevelSetHermiteGeometry] unrepresentable face-time measure");
    return value;
}
struct HalfPlaneMoments {FaceTimeRational area{0},first{0};};
// Integral of 1 and u over f0+fs*u+ft*v<0 inside [0,1]^2. A nonconstant
// equality line has zero measure; an identically zero form certifies nothing.
inline HalfPlaneMoments negativeHalfPlane(const Exact::Expansion& f0,
    const Exact::Expansion& fs,const Exact::Expansion& ft)
{
    if(Exact::sign(f0)==0 && Exact::sign(fs)==0 && Exact::sign(ft)==0)return {};
    using R=FaceTimeRational;
    struct Point{R s,t;};
    const std::array<Point,4> square{{{R(0),R(0)},{R(1),R(0)},{R(1),R(1)},{R(0),R(1)}}};
    const std::array<Exact::Expansion,4> values{{f0,Exact::add(f0,fs),Exact::add(Exact::add(f0,fs),ft),Exact::add(f0,ft)}};
    std::vector<Point> polygon;
    for(int i=0;i<4;++i){
        const int j=(i+1)%4,si=Exact::sign(values[i]),sj=Exact::sign(values[j]);
        if(si<=0)polygon.push_back(square[i]);
        if(si*sj<0){const R q(values[i],Exact::subtract(values[i],values[j]));
            polygon.push_back({square[i].s+q*(square[j].s-square[i].s),square[i].t+q*(square[j].t-square[i].t)});}
    }
    R twiceArea,sixMoment;
    for(std::size_t i=0;i<polygon.size();++i){const auto& p=polygon[i];const auto& q=polygon[(i+1)%polygon.size()];
        const auto cross=p.s*q.t-q.s*p.t;twiceArea=twiceArea+cross;sixMoment=sixMoment+(p.s+q.s)*cross;}
    return {twiceArea/R(2),sixMoment/R(6)};
}
}
