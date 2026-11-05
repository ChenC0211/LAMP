#include "Modules/Hydrodynamics/HydroDofLayout.h" // Public header must stand alone.
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <type_traits>

namespace
{
using Layout = LAMP::Hydrodynamics::HydroDofLayout;
using Values = LAMP::Hydrodynamics::HydroNativeValues;
int scenarios = 0, rejected = 0;
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action)
{
    try { action(); } catch (const std::invalid_argument&) { ++rejected; return; }
    throw std::runtime_error("invalid HYDRO layout/value input accepted");
}
DictionaryMesh input(int nx, int ny, double xmax = 1)
{
    DictionaryMesh data; data.Nx = nx; data.Ny = ny; data.maxCoord = {xmax, 1};
    data.maxInitialIterations = 0; data.maxCellNum = 20000; return data;
}
NativeLeafDomain domain(const QuadtreeMesh& mesh, bool rz)
{
    FaceSegmentGraph graph; graph.build(mesh); return NativeLeafDomain::build(mesh, graph, rz);
}
void compare(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, bool rz, bool mixed)
{
    const auto d = NativeLeafDomain::build(mesh, graph, rz);
    const auto layout = Layout::build(mesh, d);
    std::map<int, int> incidence;
    std::set<int> selectedCells;
    for (const auto& face : graph.segments())
    {
        const int row = face.isHorizontal ? mesh.vEdgeIndexNum()+mesh.hEdges(face.edgeLabel).index()
                                         : mesh.vEdges(face.edgeLabel).index();
        ++incidence[row]; selectedCells.insert(mesh.cells(face.K).index());
    }
    check(layout.velocityCount() == static_cast<int>(incidence.size()), "not one shared velocity per real face");
    check(layout.pressureCount() == static_cast<int>(selectedCells.size()), "not one pressure per selected cell");
    check(layout.primaryScalarCount() == incidence.size()+selectedCells.size(), "extra primary scalars");
    check(layout.velocityCount() <= 8*layout.pressureCount(), "nonlinear selected 2:1 DOF growth");
    check(layout.domain().meshIdentity() == d.meshIdentity(), "duplicated Core mesh identity");
    check(layout.domain().axisymmetric() == rz, "lost XY/RZ measure");
    int axisFaces = 0;
    std::map<std::tuple<int, bool, int>, int> coarseSides;
    for (int row = 0; row < d.nativeFaceCount(); ++row)
    {
        const auto f = d.face(mesh, row); const int u = layout.velocityDof(row);
        if (!incidence.count(row)) { check(u == -1, "unselected native face gained a DOF"); continue; }
        check(u >= 0 && layout.nativeFace(u) == row, "velocity numbering is not bijective");
        check(layout.normalAxis(u) == (f.horizontal ? Layout::NormalAxis::Second : Layout::NormalAxis::First),
            "native positive-axis component mismatch");
        if (f.area == 0) { ++axisFaces; check(rz && !f.horizontal && f.center.x == 0, "spurious zero face measure"); }
        const bool interior = f.negativeCell >= 0 && f.positiveCell >= 0;
        check(incidence[row] == (interior ? 2 : 1), "internal velocity DOF is not shared by both cells");
        if (interior)
        {
            const auto& a = mesh.cells(mesh.cellIndexToLabel[f.negativeCell]);
            const auto& b = mesh.cells(mesh.cellIndexToLabel[f.positiveCell]);
            if (a.level() != b.level())
            {
                const auto& coarse = a.level() < b.level() ? a : b;
                const double tangent = f.horizontal ? coarse.size().x : coarse.size().y;
                check(f.tangentSpacing == .5*tangent, "coarse face is not represented by fine subface");
                const auto c = coarse.centerCoord();
                const int side = (f.horizontal ? f.center.y > c.y : f.center.x > c.x) ? 1 : -1;
                ++coarseSides[{coarse.index(), f.horizontal, side}];
            }
        }
    }
    if (mixed) { check(!coarseSides.empty(), "test did not exercise coarse/fine sharing");
        for (const auto& side : coarseSides) check(side.second == 2, "extra/missing coarse-side velocity unknown"); }
    if (rz) check(axisFaces > 0, "RZ axis native slots were dropped");
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        const int p = layout.pressureDof(ci);
        if (!selectedCells.count(ci)) check(p == -1, "unselected cell gained pressure");
        else check(p >= 0 && layout.nativeCell(p) == ci, "pressure numbering is not bijective");
    }
    std::cout << "HYDRO_LAYOUT," << rz << ',' << mixed << ',' << mesh.cellIndexNum() << ','
              << layout.velocityCount() << ',' << layout.pressureCount() << ',' << layout.primaryScalarCount() << '\n';
    ++scenarios;
}
void uniform()
{
    for (bool rz : {false, true}) for (int n : {2, 4, 8, 16, 32})
    {
        QuadtreeMesh mesh(input(n, n)); FaceSegmentGraph graph; graph.build(mesh);
        const auto layout = Layout::build(mesh, NativeLeafDomain::build(mesh, graph, rz));
        check(layout.pressureCount() == n*n && layout.velocityCount() == 2*n*(n+1), "uniform analytic DOF count");
        compare(mesh, graph, rz, false);
    }
    QuadtreeMesh rectangular(input(4, 8, 2));
    const auto layout = Layout::build(rectangular, domain(rectangular, false));
    check(layout.pressureCount() == 32 && layout.velocityCount() == 4*9+5*8, "rectangular analytic DOF count");
    ++scenarios;
}
void amr()
{
    for (bool rz : {false, true}) for (int n : {4, 8, 16})
    {
        QuadtreeMesh mesh(input(n, n));
        DictionaryMesh::RefineBox box; box.min = {.25, .25}; box.max = {.75, .75}; box.level = 1;
        mesh.applyRefineBoxes({box}); FaceSegmentGraph graph; graph.build(mesh);
        compare(mesh, graph, rz, true);
        FaceSegmentGraph selected;
        selected.build(mesh, [](const QuadtreeMesh&, const Cell& cell) { return cell.centerCoord().x < .5; });
        compare(mesh, selected, rz, true);
    }
}
void guards()
{
    static_assert(!std::is_default_constructible_v<Layout>);
    static_assert(std::is_same_v<decltype(std::declval<const Values&>().normalVelocity()), const std::vector<double>&>);
    static_assert(std::is_same_v<decltype(std::declval<const Values&>().pressure()), const std::vector<double>&>);
    QuadtreeMesh mesh(input(4, 4)); auto d = domain(mesh, false); const auto layout = Layout::build(mesh, d);
    auto u = std::vector<double>(layout.velocityCount(), 2), p = std::vector<double>(layout.pressureCount(), 3);
    const Values values(layout, u, p); u[0] = 19; p[0] = -7;
    check(values.normalVelocity()[0] == 2 && values.pressure()[0] == 3, "values alias caller input");
    const auto copied = layout; check(copied.sameIdentity(layout), "copy changed layout identity");
    values.validateFor(copied, mesh);
    const auto rebuilt = Layout::build(mesh, d);
    check(!rebuilt.sameIdentity(layout), "fresh layout implicitly reused old field identity");
    rejects([&] { values.validateFor(rebuilt, mesh); });
    const auto rz = Layout::build(mesh, domain(mesh, true)); rejects([&] { values.validateFor(rz, mesh); });
    rejects([&] { Layout::build(mesh, NativeLeafDomain{}); });
    for (int row : {-1, d.nativeFaceCount()}) rejects([&] { layout.velocityDof(row); });
    for (int row : {-1, layout.velocityCount()}) {
        rejects([&] { layout.nativeFace(row); }); rejects([&] { layout.normalAxis(row); }); }
    for (int row : {-1, mesh.cellIndexNum()}) rejects([&] { layout.pressureDof(row); });
    for (int row : {-1, layout.pressureCount()}) rejects([&] { layout.nativeCell(row); });
    auto wrongU = u; wrongU.pop_back(); rejects([&] { Values bad(layout, wrongU, p); });
    wrongU = u; wrongU.push_back(0); rejects([&] { Values bad(layout, wrongU, p); });
    auto wrongP = p; wrongP.pop_back(); rejects([&] { Values bad(layout, u, wrongP); });
    wrongP = p; wrongP.push_back(0); rejects([&] { Values bad(layout, u, wrongP); });
    for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()}) {
        wrongU = u; wrongU[0] = bad; rejects([&] { Values invalid(layout, wrongU, p); });
        wrongP = p; wrongP[0] = bad; rejects([&] { Values invalid(layout, u, wrongP); }); }
    FaceSegmentGraph left, right;
    left.build(mesh, [](const QuadtreeMesh&, const Cell& c) { return c.centerCoord().x < .5; });
    right.build(mesh, [](const QuadtreeMesh&, const Cell& c) { return c.centerCoord().x > .5; });
    const auto l = Layout::build(mesh, NativeLeafDomain::build(mesh, left, false));
    const auto r = Layout::build(mesh, NativeLeafDomain::build(mesh, right, false));
    check(l.velocityCount() == r.velocityCount() && l.pressureCount() == r.pressureCount(), "selection fixture dimensions differ");
    const Values selected(l, std::vector<double>(l.velocityCount()), std::vector<double>(l.pressureCount()));
    rejects([&] { selected.validateFor(r, mesh); });
    QuadtreeMesh stretched(input(4, 4, 2));
    check(stretched.cellIndexNum() == mesh.cellIndexNum(), "same-count stale mesh fixture failed");
    rejects([&] { layout.validateMesh(stretched); }); rejects([&] { values.validateFor(layout, stretched); });
    auto& vertex = mesh.verts(mesh.vertIndexToLabel[0]); const auto saved = vertex.coord();
    vertex.setCoord(saved+Vector2D{.01, 0});
    rejects([&] { layout.validateMesh(mesh); }); rejects([&] { Layout::build(mesh, d); });
    vertex.setCoord(saved); values.validateFor(layout, mesh);
    DictionaryMesh::RefineBox box; box.min = {.25, .25}; box.max = {.75, .75}; box.level = 1;
    mesh.applyRefineBoxes({box}); rejects([&] { values.validateFor(layout, mesh); });
}
void lifetime()
{
    std::optional<Values> retained;
    {
        QuadtreeMesh mesh(input(4, 4)); auto d = domain(mesh, false); const auto layout = Layout::build(mesh, d);
        retained.emplace(layout, std::vector<double>(layout.velocityCount(), .7), std::vector<double>(layout.pressureCount(), -2));
        d = NativeLeafDomain{};
        check(retained->layout().pressureCount() == 16, "layout borrowed caller domain");
    }
    check(retained->normalVelocity().size() == 40 && retained->pressure().size() == 16 &&
          retained->layout().domain().cellVolumes().size() == 16 && retained->layout().nativeCell(15) >= 0,
          "retained native values/layout lost lifetime");
    QuadtreeMesh replacement(input(4, 4)); retained->validateFor(retained->layout(), replacement);
}
}
int main() try
{
    uniform(); amr(); guards(); lifetime();
    std::cout << "H01 HYDRO native layout PASS scenarios=" << scenarios << " rejected=" << rejected << '\n';
    return 0;
}
catch (const std::exception& error) { std::cerr << "H01 FAIL: " << error.what() << '\n'; return 1; }
