#include "Modules/Hydrodynamics/HydroNodalTraction.h"
#include "Modules/Interface/Compatibility/UniformInterfaceContour.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <type_traits>

namespace {
void check(bool condition, const char* reason)
{ if (!condition) throw std::runtime_error(reason); }
template<class F> void rejects(F action)
{
    try { action(); } catch (const std::exception&) { return; }
    throw std::runtime_error("invalid contour/load was accepted");
}
DictionaryGeometry circle()
{
    DictionaryGeometry input;
    input.points = {{.73,.5}, {.5,.73}, {.27,.5}, {.5,.27}};
    for (int i=0; i<4; ++i)
    {
        Line line;
        line.points[0]=input.points[i]; line.points[1]=input.points[(i+1)%4];
        line.setType(LineType::arc); line.setArcCenter({.5,.5});
        line.setRadius(.23); line.setRadian(.5);
        line.setGeoBoundary(GeometricBoundary::moving); line.setImmersed(true);
        input.lines.push_back(line);
    }
    input.levelSetName={"liquid"}; input.innerPoints={{.5,.5}}; input.moving={true};
    input.linesOnInterface={{{0,false},{1,false},{2,false},{3,false}}};
    return input;
}
void exercise(bool rz)
{
    DictionaryMesh input; input.Nx=input.Ny=16; input.maxCellNum=1024; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input), foreignMesh(input);
    auto geometry=std::make_unique<ImmersedGeometry<QuadtreeMesh>>(circle(),mesh);
    geometry->initializeLevelSets(mesh,0.,false);
    geometry->refreshLevelSetDerivedData(mesh,0.,false);
    auto contour=UniformInterfaceContour::build(mesh,*geometry,0,rz);
    const auto report=contour.report();
    check(report.nodeCount>3 && std::isfinite(report.maxAbsCurvature),"closed contour missing");
    const double analyticMeasure=2.*Pi*.23*(rz?2.*Pi*.5:1.);
    check(std::abs(report.surfaceMeasure/analyticMeasure-1.)<.08,"circle/torus surface measure");
    const auto* storage=contour.nodes().data();
    const auto first=contour.nodes().front().coord();
    // Owning old geometry must not turn into a view of updated source fields.
    for(int i=0;i<mesh.vertIndexNum();++i)
        geometry->levelSets[0].value[i]=mesh.verts(mesh.vertIndexToLabel[i]).coord().x-.5;
    geometry->refreshLevelSetDerivatives(mesh);
    rejects([&]{UniformInterfaceContour::build(mesh,*geometry,0,rz);});
    geometry.reset();
    check(contour.nodes().data()==storage && contour.nodes().front().coord().x==first.x
        && contour.nodes().front().coord().y==first.y,"source mutation/destruction changed snapshot");
    auto retained=std::move(contour);
    check(retained.nodes().data()==storage,"moving contour copied nodes");
    rejects([&]{contour.validateMesh(mesh);});
    Field<double> tx(report.nodeCount,2.),ty(report.nodeCount,-.75);
    CellField<double> fx(mesh,.125),fy(mesh,-.25);
    auto unchanged=[&] {
        for(int i=0;i<mesh.cellIndexNum();++i)
            check(fx.value[i]==.125 && fy.value[i]==-.25,"failed load partially changed output");
    };
    rejects([&]{HydroNodalTraction::spread(foreignMesh,retained,tx,ty,fx,fy);}); unchanged();
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]); const auto saved=vertex.coord();
    vertex.setCoord(saved+Vector2D(.001,0.));
    rejects([&]{HydroNodalTraction::spread(mesh,retained,tx,ty,fx,fy);}); unchanged();
    vertex.setCoord(saved);
    tx.value[report.nodeCount-1]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{HydroNodalTraction::spread(mesh,retained,tx,ty,fx,fy);}); unchanged();
    tx.value[report.nodeCount-1]=2.;
    Field<double> shortField(report.nodeCount-1,0.);
    rejects([&]{HydroNodalTraction::spread(mesh,retained,shortField,ty,fx,fy);}); unchanged();
    rejects([&]{HydroNodalTraction::spread(mesh,retained,tx,ty,fx,fx);}); unchanged();
    const auto load=HydroNodalTraction::spread(mesh,retained,tx,ty,fx,fy);
    check(std::abs(load.forceX-2.*report.surfaceMeasure)<1e-12
        && std::abs(load.forceY+.75*report.surfaceMeasure)<1e-12,"constant traction surface integral");
    long double totalX=0.,totalY=0.; double maximum=0.;
    for(int i=0;i<mesh.cellIndexNum();++i)
    {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);
        const double dx=cell.size().x,dy=cell.size().y,r=cell.centerCoord().x;
        const double volume=rz?Pi*((r+.5*dx)*(r+.5*dx)-(r-.5*dx)*(r-.5*dx))*dy:dx*dy;
        totalX+=(fx.value[i]-.125)*volume;
        totalY+=(fy.value[i]+.25)*volume;
        maximum=std::max(maximum,std::abs(fx.value[i]-.125));
    }
    check(maximum>0. && std::abs(static_cast<double>(totalX)-load.forceX)<1e-12
        && std::abs(static_cast<double>(totalY)-load.forceY)<1e-12,"Cartesian/annular-volume force conservation");
    check(std::abs(load.spreadForceX-load.forceX)<1e-12 && std::abs(load.spreadForceY-load.forceY)<1e-12,
        "spread report conservation");
}
}
int main() try
{
    static_assert(!std::is_copy_constructible_v<UniformInterfaceContour>);
    static_assert(std::is_same_v<decltype(std::declval<const UniformInterfaceContour&>().nodes()[0]),const InterfaceNode&>);
    exercise(false); exercise(true);
    std::cout<<"Uniform Interface snapshot and XY/RZ nodal traction checks passed\n";
    return 0;
}
catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
