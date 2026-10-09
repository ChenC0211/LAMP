#include "Core/Base/ScalarExpression.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using E=Quantum::ScalarExpression;
static void near(double value,double reference){if(std::abs(value-reference)>1e-12*std::max(1.,std::abs(reference)))throw std::runtime_error("expression value/chain-rule mismatch");}
int main()
{
    try{
        const E spatial("sin(x)*exp(y)+x*x/y",{},E::Coordinates::Planar);
        const auto a=spatial.sample(.3,.7);near(a.value,std::sin(.3)*std::exp(.7)+.09/.7);
        near(a.derivative,std::cos(.3)*std::exp(.7)+.6/.7);near(a.partialSecond,std::sin(.3)*std::exp(.7)-.09/.49);
        const E moving("(x-speed*t)^2+sin(y+t)",{{"speed",.2}},E::Coordinates::PlanarTime);
        const auto b=moving.sample(.8,.4,1.5);near(b.value,.25+std::sin(1.9));near(b.derivative,1);near(b.partialSecond,std::cos(1.9));
        near(moving.sample(.8,.4).value,moving.sample(.8,.4,0).value);
        const E rz("r*r+(z-t)^2",{},E::Coordinates::MeridionalTime);
        const auto c=rz.sample(.2,.7,.3);near(c.value,.2);near(c.derivative,.4);near(c.partialSecond,.8);
        const E temporal("3*cos(pi*t/2)",{},E::Coordinates::Time);
        near(temporal.sample(1).value,0);near(temporal.sample(1).derivative,0);near(temporal.sample(1).partialSecond,0);
        const E parameter("s^3-2*s");near(parameter.sample(.5).value,-.875);near(parameter.sample(.5).derivative,-1.25);
        int rejected=0;for(const auto& text:{"x+t","unknown(x)","sin(x);","x^y"})try{E invalid(text,{},E::Coordinates::Planar);}catch(const std::invalid_argument&){++rejected;}
        if(rejected!=4)throw std::runtime_error("invalid/unsupported expression accepted");
        try{E invalid("1/(x-y)",{},E::Coordinates::Planar);invalid.sample(.3,.3);throw std::runtime_error("singular expression accepted");}catch(const std::domain_error&){}
        std::cout<<"scalar expression spatial/time contract PASS\n";return 0;
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
