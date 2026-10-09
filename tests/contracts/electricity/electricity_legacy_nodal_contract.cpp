#include "Modules/Electricity/LegacyNodalMaxwell.h"
#include <cmath>
#include <iostream>

int main() try
{
    std::vector<InterfaceNode> nodes(2);
    nodes[1].setNormal({1.,0.}); nodes[0].setNormal({0.,1.});
    Field<double> normal(2),tangent(2),fx(2),fy(2);
    normal.value[0]=2.; normal.value[1]=-1.;
    tangent.value[0]=.5; tangent.value[1]=1.5;
    fx.value[0]=.25; fx.value[1]=0.; fy.value[0]=0.; fy.value[1]=-.375;
    auto check=[&](double x0,double x1,double y0,double y1) {
        const double actual[]={fx.value[0],fx.value[1],fy.value[0],fy.value[1]};
        const double expected[]={x0,x1,y0,y1};
        for (int i=0;i<4;++i)
            if (!std::isfinite(actual[i]) || std::abs(actual[i]-expected[i])>1e-12)
                throw std::runtime_error("Maxwell normal/tangent orientation or additive load");
    };
    LegacyNodalMaxwell::add({1,0},nodes,normal,tangent,fx,fy);
    check(2.25,-1.5,.5,-1.375);
    // A second contribution must accumulate, preserving pre-existing capillarity.
    LegacyNodalMaxwell::add({1,0},nodes,normal,tangent,fx,fy);
    check(4.25,-3.,1.,-2.375);
    LegacyNodalMaxwell::add({},nodes,normal,tangent,fx,fy);
    check(4.25,-3.,1.,-2.375);
    std::cout<<"Legacy Maxwell rotation and accumulation checks passed\n";
    return 0;
}
catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
