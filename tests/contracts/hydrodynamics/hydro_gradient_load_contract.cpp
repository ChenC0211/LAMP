#include "Modules/Hydrodynamics/HydroFittedScalarGradientLoad.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/IO/DictionaryMesh.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action) {
    try { action(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("invalid gradient load accepted");
}
void run(bool axisymmetric) {
    DictionaryMesh input; input.Nx=input.Ny=8; input.maxCellNum=256; input.maxInitialIterations=0;
    QuadtreeMesh mesh(input); FaceSegmentGraph graph; graph.build(mesh);
    std::vector<HydroFittedDomain::Boundary> boundary(graph.segments().size());
    for(std::size_t s=0;s<boundary.size();++s) if(graph.segments()[s].L<0)
        boundary[s].type=graph.segments()[s].side==FaceSide::North
            ?HydroFittedDomain::BoundaryType::PressureOutlet:HydroFittedDomain::BoundaryType::Symmetry;
    const auto domain=HydroFittedDomain::build(mesh,graph,boundary,axisymmetric);
    std::vector<double> scalar(mesh.cellIndexNum()), coefficient(domain.faces().size(),0.);
    for(int k=0;k<mesh.cellIndexNum();++k)
        scalar[k]=mesh.cells(mesh.cellIndexToLabel[k]).centerCoord().y<.5?1.:2.;
    // Analytic planar Maxwell jump for Dn=1, Et=1/4, epsilon=1 -> 2.
    // This receiver test has no Electricity dependency.
    const double stressJump=-.28125, scale=3.5;
    for(const auto& face:domain.faces())
        if(face.negativeCell>=0 && face.positiveCell>=0 && scalar[face.negativeCell]!=scalar[face.positiveCell])
            coefficient[face.nativeRow]=stressJump;
    const auto result=buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,coefficient,scale);
    check(result.faces.size()==8,"wrong active load face count");
    for(const auto& face:result.faces)
        check(std::abs(face.forceDensity-stressJump*8.*scale)<1e-12,"incorrect force density or normalization");
    const double area=axisymmetric?std::acos(-1.):1.;
    check(result.integralX==0. && std::abs(result.integralY-stressJump*scale*area)<1e-12,
        "integrated load does not equal analytic stress jump times XY/RZ area");
    const auto saved=scalar;
    const auto wrongGeometry=NativeLeafDomain::build(mesh,graph,!axisymmetric);
    rejects([&]{buildFittedScalarGradientLoad(mesh,domain,wrongGeometry,scalar,coefficient,scale);});
    FaceSegmentGraph partial; partial.build(mesh,[](const QuadtreeMesh&,const Cell& cell){return cell.centerCoord().y<.5;});
    const auto wrongSelection=NativeLeafDomain::build(mesh,partial,axisymmetric);
    rejects([&]{buildFittedScalarGradientLoad(mesh,domain,wrongSelection,scalar,coefficient,scale);});
    auto bad=coefficient; bad[0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,bad,scale);});
    bad=coefficient;
    for(const auto& face:domain.faces()) if(face.negativeCell<0 || face.positiveCell<0) {bad[face.nativeRow]=1.;break;}
    rejects([&]{buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,bad,scale);});
    bad.pop_back();
    rejects([&]{buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,bad,scale);});
    check(scalar==saved,"rejected load changed caller state");
    scalar.assign(scalar.size(),1.);
    const auto uniform=buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,coefficient,scale);
    check(uniform.faces.empty() && uniform.maxAbsComponent==0.,"constant scalar creates a load");
    auto& vertex=mesh.verts(mesh.vertIndexToLabel[0]);
    vertex.setCoord(vertex.coord()+Vector2D{.01,0.});
    rejects([&]{buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,coefficient,scale);});
}
void nativeAffine(bool axisymmetric) {
    DictionaryMesh input;input.Nx=input.Ny=8;input.maxCellNum=256;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);DictionaryMesh::RefineBox box;box.min={0.,.5};box.max={1.,1.};box.level=1;
    mesh.applyRefineBoxes({box});check(mesh.cellIndexNum()==160,"native load test did not refine the mesh");
    FaceSegmentGraph graph;graph.build(mesh);
    std::vector<HydroFittedDomain::Boundary> boundary(graph.segments().size());
    for(std::size_t s=0;s<boundary.size();++s)if(graph.segments()[s].L<0)
        boundary[s].type=graph.segments()[s].side==FaceSide::North
            ?HydroFittedDomain::BoundaryType::PressureOutlet:HydroFittedDomain::BoundaryType::Symmetry;
    const auto domain=HydroFittedDomain::build(mesh,graph,boundary,axisymmetric);
    std::vector<double> scalar(mesh.cellIndexNum()),coefficient(domain.faces().size(),0.);
    int transitions=0;
    for(const auto& face:domain.faces())if(face.negativeCell>=0 && face.positiveCell>=0) {
        coefficient[face.nativeRow]=2.;
        transitions+=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).level()
            !=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).level();
    }
    check(transitions==16,"wrong number of real coarse/fine faces");
    for(bool affine:{false,true}) {
        for(int k=0;k<mesh.cellIndexNum();++k) {
            const auto p=mesh.cells(mesh.cellIndexToLabel[k]).centerCoord();
            scalar[k]=affine?10.+.25*p.x+p.y:10.;
        }
        const auto load=buildFittedScalarGradientLoad(mesh,domain,domain.geometry(),scalar,coefficient,3.5);
        long double expectedX=0.,expectedY=0.;
        for(const auto& face:domain.faces())if(face.velocityRow>=0 && face.negativeCell>=0 && face.positiveCell>=0) {
            const double expected=affine?7.*(face.horizontal?1.:.25):0.;
            const double actual=face.horizontal?load.forceY.value[face.nativeRow-mesh.vEdgeIndexNum()]:load.forceX.value[face.nativeRow];
            check(std::abs(actual-expected)<1e-11,"native load differs from analytic constant/affine gradient");
            (face.horizontal?expectedY:expectedX)+=static_cast<long double>(expected)*face.dualVolume;
        }
        check(std::abs(load.integralX-static_cast<double>(expectedX))<1e-11
            &&std::abs(load.integralY-static_cast<double>(expectedY))<1e-11,"native XY/RZ load integral disagrees");
    }
    std::cout<<"Native "<<(axisymmetric?"RZ":"XY")<<" scalar loads: 160 cells, 16 coarse/fine faces; constant/affine reference passed\n";
}
}
int main() try {
    run(false);run(true);
    nativeAffine(false);nativeAffine(true);
    std::cout<<"Hydro-only pressure-compatible XY/RZ load, integral and source rejection passed\n";
    return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
