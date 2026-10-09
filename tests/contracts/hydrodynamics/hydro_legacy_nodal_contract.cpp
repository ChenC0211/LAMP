#include "Modules/Hydrodynamics/LegacyNodalForcing.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <iostream>

namespace {
void near(double actual, double expected, const char* reason)
{
    if (!std::isfinite(actual) || std::abs(actual-expected)>1e-12)
        throw std::runtime_error(reason);
}
void field(const Field<double>& actual, std::initializer_list<double> expected, const char* reason)
{
    if (actual.size()!=static_cast<int>(expected.size())) throw std::runtime_error(reason);
    int i=0; for (double value:expected) near(actual.value[i++],value,reason);
}
void cells(const CellField<double>& actual, double initial,
           std::initializer_list<double> increments, const char* reason)
{
    int i=0; for (double value:increments) near(actual.value[i++],initial+value,reason);
    for (;i<actual.size();++i) near(actual.value[i],initial,reason);
}
}

int main() try
{
    DictionaryMesh input; input.Nx=input.Ny=4; input.maxCellNum=64; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);
    std::vector<InterfaceNode> nodes(2);
    const std::vector<int> labels={1,0};
    for (auto& node:nodes)
        for (int j=0;j<25;++j) { node.cellStencil[j]=0; node.cellStencilCoeff[j]=0.; }
    nodes[1].setIndex(0); nodes[1].setNormal({1.,0.}); nodes[1].setCurvature(2.);
    nodes[1].setDs(.25);
    nodes[1].cellStencil[0]=0; nodes[1].cellStencilCoeff[0]=.25;
    nodes[1].cellStencil[1]=1; nodes[1].cellStencilCoeff[1]=.75;
    nodes[0].setIndex(1); nodes[0].setNormal({0.,1.}); nodes[0].setCurvature(-3.);
    nodes[0].setDs(.5);
    nodes[0].cellStencil[0]=1; nodes[0].cellStencilCoeff[0]=.4;
    nodes[0].cellStencil[1]=2; nodes[0].cellStencilCoeff[1]=.6;

    Field<double> fx(7,99.),fy(7,99.),ux,uy;
    LegacyNodalForcing::initializeCapillary(labels,nodes,.125,fx,fy);
    field(fx,{.25,0.},"capillary x and reset"); field(fy,{0.,-.375},"capillary y and reset");
    // Distinct, overlapping stencils check label order, normalization and addition.
    CellField<double> cellUx(mesh,3.),cellUy(mesh,-2.);
    cellUx.value[1]=7.; cellUx.value[2]=12.;
    cellUy.value[1]=2.; cellUy.value[2]=7.;
    LegacyNodalForcing::interpolate(labels,nodes,cellUx,cellUy,ux,uy);
    field(ux,{6.,10.},"mapped interpolation x"); field(uy,{1.,5.},"mapped interpolation y");
    const std::vector<InterfaceNode> noSlip={nodes[1],nodes[0]};
    LegacyNodalForcing::interpolateNoSlip(noSlip,2,cellUx,cellUy,ux,uy);
    field(ux,{6.,10.},"no-slip interpolation x"); field(uy,{1.,5.},"no-slip interpolation y");
    cellUx.initialize(mesh,3.); cellUy.initialize(mesh,-2.);
    LegacyNodalForcing::interpolateNoSlip(noSlip,2,cellUx,cellUy,ux,uy);
    field(ux,{3.,3.},"constant preservation x"); field(uy,{-2.,-2.},"constant preservation y");
    near(LegacyNodalForcing::fixedVelocityForce(2,ux,uy,{.5,-.25},4.,fx,fy),
         std::sqrt(9.3125),"maximum velocity residual");
    field(fx,{-10.,-10.},"fixed-velocity x"); field(fy,{7.,7.},"fixed-velocity y");
    CellField<double> volumeX(mesh,.125),volumeY(mesh,-.25);
    LegacyNodalForcing::spreadVolume(noSlip,2,fx,fy,volumeX,volumeY);
    cells(volumeX,.125,{-2.5,-11.5,-6.},"volume x addition");
    cells(volumeY,-.25,{1.75,8.05,4.2},"volume y addition");
    fx.value[0]=2.25; fx.value[1]=-1.5; fy.value[0]=.5; fy.value[1]=-1.375;
    CellField<double> surfaceX(mesh,.125),surfaceY(mesh,-.25);
    LegacyNodalForcing::spreadSurface(mesh,labels,nodes,fx,fy,surfaceX,surfaceY);
    cells(surfaceX,.125,{2.25,1.95,-7.2},"surface x area/quadrature");
    cells(surfaceY,-.25,{.5,-2.9,-6.6},"surface y area/quadrature");
    double totalX=0.,totalY=0.;
    for (int i=0;i<mesh.cellIndexNum();++i)
    {
        const double area=mesh.cells(mesh.cellIndexToLabel[i]).area();
        totalX+=(surfaceX.value[i]-.125)*area;
        totalY+=(surfaceY.value[i]+.25)*area;
    }
    near(totalX,-.1875,"surface force integral x"); near(totalY,-.5625,"surface force integral y");
    CellField<double> evapX(mesh,.125),evapY(mesh,-.25);
    LegacyNodalForcing::spreadEvaporation(labels,nodes,{3.,2.},evapX,evapY);
    cells(evapX,.125,{-.5,-1.5,0.},"evaporation sign x");
    cells(evapY,-.25,{0.,-1.2,-1.8},"evaporation sign y");
    LegacyNodalForcing::initializeCapillary({},nodes,.125,fx,fy);
    field(fx,{},"empty capillary"); field(fy,{},"empty capillary");
    near(LegacyNodalForcing::fixedVelocityForce(0,ux,uy,{0.,0.},4.,fx,fy),0.,"empty residual");
    std::cout<<"Legacy capillary, interpolation, direct force and spreading checks passed\n";
    return 0;
}
catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
