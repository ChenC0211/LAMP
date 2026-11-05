#include "Modules/Hydrodynamics/HydroIntegralDivergence.h" // Public header stands alone.
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <type_traits>

namespace
{
using Layout = LAMP::Hydrodynamics::HydroDofLayout;
using Values = LAMP::Hydrodynamics::HydroNativeValues;
using Divergence = LAMP::Hydrodynamics::HydroIntegralDivergence;
constexpr long double tolerance = 1e-13L;
const long double pi = std::acos(-1.0L);
int scenarios = 0, randomFields = 0, rejected = 0, coarseSides = 0;
long double metricDefect = 0, localDefect = 0, globalDefect = 0, coarseDefect = 0;
long double analyticDefect = 0, streamDefect = 0;
void check(bool ok, const char* reason) { if (!ok) throw std::runtime_error(reason); }
template<class F> void rejects(F action)
{
    try { action(); } catch (const std::invalid_argument&) { ++rejected; return; }
    throw std::runtime_error("invalid divergence source/output accepted");
}
void observe(long double error, long double scale, long double& peak)
{
    check(std::isfinite(error) && std::isfinite(scale) && scale >= 0, "invalid independent error scale");
    if (scale == 0) { check(error == 0, "zero-scale flux is not exactly zero"); return; }
    const auto defect = std::abs(error)/scale;
    peak = std::max(peak, defect);
    if (defect > tolerance)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "closure exceeds frozen 1e-13 gate: defect=" << defect
                << " error=" << error << " scale=" << scale << " completed_scenarios=" << scenarios;
        throw std::runtime_error(message.str());
    }
}
DictionaryMesh input(int nx, int ny, int metric)
{
    DictionaryMesh data; data.Nx = nx; data.Ny = ny;
    data.minCoord = {metric == 0 ? -.5 : metric == 1 ? 0.0 : .25, -.25};
    data.maxCoord = {metric == 0 ? 1.5 : metric == 1 ? 1.0 : 1.25, 1.25};
    data.maxInitialIterations = 0; data.maxCellNum = 20000; return data;
}
void refine(QuadtreeMesh& mesh, int depth)
{
    std::vector<DictionaryMesh::RefineBox> boxes;
    for (int level = 1; level <= depth; ++level)
    {
        const double lo = level == 1 ? .25 : .375, hi = 1-lo;
        const auto width = mesh.maxCoord-mesh.minCoord;
        DictionaryMesh::RefineBox box;
        box.min = mesh.minCoord+width*lo; box.max = mesh.minCoord+width*hi; box.level = level;
        boxes.push_back(box);
    }
    if (!boxes.empty()) mesh.applyRefineBoxes(boxes);
}
FaceSegmentGraph selected(const QuadtreeMesh& mesh, int selection)
{
    FaceSegmentGraph graph;
    graph.build(mesh, [selection](const QuadtreeMesh& m, const Cell& c) {
        const auto p = c.centerCoord(), w = m.maxCoord-m.minCoord;
        const double x = (p.x-m.minCoord.x)/w.x, y = (p.y-m.minCoord.y)/w.y;
        if (selection == 1) return x < .5;
        if (selection == 2) return x < .25 || x > .75;
        if (selection == 3) return !(x > .375 && x < .625 && y > .375 && y < .625);
        return true;
    });
    return graph;
}
struct Face
{
    long double x0, x1, y0, y1, area;
    bool horizontal;
};
Face face(const QuadtreeMesh& mesh, const FaceSegment& segment, bool rz)
{
    const auto& edge = segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
    const auto a = mesh.verts(edge.edgeToVert[0]).coord(), b = mesh.verts(edge.edgeToVert[1]).coord();
    Face f{std::min(a.x,b.x),std::max(a.x,b.x),std::min(a.y,b.y),std::max(a.y,b.y),0,segment.isHorizontal};
    f.area = f.horizontal ? f.x1-f.x0 : f.y1-f.y0;
    // Independent endpoint formulas, not CartesianMeasure or domain.face().
    if (rz) f.area = f.horizontal ? pi*(f.x1*f.x1-f.x0*f.x0) : 2*pi*f.x0*(f.y1-f.y0);
    return f;
}
int velocityDof(const QuadtreeMesh& mesh, const Layout& layout, const FaceSegment& segment)
{
    const int native = segment.isHorizontal ? mesh.vEdgeIndexNum()+mesh.hEdges(segment.edgeLabel).index()
                                          : mesh.vEdges(segment.edgeLabel).index();
    const int u = layout.velocityDof(native); check(u >= 0, "selected graph face missing from layout"); return u;
}
long double sign(const FaceSegment& segment) { return segment.isHorizontal ? segment.n.y : segment.n.x; }
long double radialMean(long double a, long double b)
{ return (2.0L/3)*(a*a+a*b+b*b)/(a+b); }
long double psi(long double x, long double y, bool rz)
{ return x*x*y*y+(rz ? x*x*x*x : x*x*x)*y; }
enum class Field { Linear, Quadratic, Stream };
double normalAverage(const Face& f, bool rz, Field field)
{
    const long double x = (f.x0+f.x1)/2, y = (f.y0+f.y1)/2;
    constexpr long double a=.7L, b=-.3L, d=.2L, e=.45L, k=-.25L;
    if (field == Field::Stream)
    {
        if (f.area == 0) return 0; // Analytic regular radial-axis limit.
        const auto flow = f.horizontal ? -(psi(f.x1,y,rz)-psi(f.x0,y,rz))
                                       : psi(x,f.y1,rz)-psi(x,f.y0,rz);
        return static_cast<double>((rz ? 2*pi : 1)*flow/f.area);
    }
    if (!f.horizontal)
        return static_cast<double>(field == Field::Linear ? (rz ? a*x : a*x+b*y+.4L) : a*x+b*x*y);
    const auto mean = rz ? radialMean(f.x0,f.x1) : x;
    if (field == Field::Linear) return static_cast<double>(d*mean+e*y-.2L);
    const auto second = rz ? (f.x0*f.x0+f.x1*f.x1)/2
                           : (f.x0*f.x0+f.x0*f.x1+f.x1*f.x1)/3;
    return static_cast<double>(d*second+e*y*y+k*mean*y);
}
std::vector<double> fieldValues(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const Layout& layout, bool rz, Field field)
{
    std::vector<double> u(layout.velocityCount(),std::numeric_limits<double>::quiet_NaN());
    for (const auto& segment : graph.segments())
        u[velocityDof(mesh,layout,segment)] = normalAverage(face(mesh,segment,rz),rz,field);
    for (double value : u) check(std::isfinite(value), "independent field left an unfilled face");
    return u;
}
struct Reference
{
    std::vector<long double> integral, scale;
    long double boundary=0, totalScale=0;
};
Reference reference(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const Layout& layout, bool rz, const std::vector<double>& u)
{
    Reference result;
    result.integral.assign(layout.pressureCount(),0); result.scale.assign(layout.pressureCount(),0);
    for (const auto& segment : graph.segments())
    {
        const int p = layout.pressureDof(mesh.cells(segment.K).index());
        const auto q = face(mesh,segment,rz).area*u[velocityDof(mesh,layout,segment)]*sign(segment);
        result.integral[p] += q; result.scale[p] += std::abs(q); result.totalScale += std::abs(q);
        // L retains the physical neighbor even if it was not selected. A
        // selected-domain boundary also includes those selection interfaces.
        const auto neighbor = segment.L < 0 ? std::pair<std::size_t,std::size_t>{0,0}
            : graph.segmentRangeForCellIndex(mesh.cells(segment.L).index());
        if (neighbor.first == neighbor.second) result.boundary += q;
    }
    return result;
}
void coefficients(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, const Divergence& op, bool rz)
{
    const auto& layout = op.layout(); const auto& matrix = op.matrix();
    check(matrix.rows() == layout.pressureCount() && matrix.cols() == layout.velocityCount() && matrix.isCompressed(),
        "wrong/uncompressed native integral matrix dimensions");
    check(matrix.nonZeros() <= 2*static_cast<Eigen::Index>(layout.velocityCount()) &&
          matrix.nonZeros() <= 8*static_cast<Eigen::Index>(layout.pressureCount()), "nonlinear integral storage growth");
    std::vector<int> incidence(layout.velocityCount(),0), actual(layout.velocityCount(),0);
    std::vector<long double> columnSum(layout.velocityCount(),0), expected(layout.velocityCount(),0), area(layout.velocityCount(),0);
    Eigen::Index nonzeros=0;
    for (const auto& segment : graph.segments())
    {
        const int p = layout.pressureDof(mesh.cells(segment.K).index()), u = velocityDof(mesh,layout,segment);
        const auto f = face(mesh,segment,rz); ++incidence[u]; area[u]=f.area;
        expected[u] += sign(segment)*f.area;
        observe(matrix.coeff(p,u)-sign(segment)*f.area,f.area,metricDefect);
        nonzeros += f.area != 0;
    }
    check(matrix.nonZeros() == nonzeros, "missing/extra coarse or zero-axis matrix entry");
    for (Eigen::Index p=0; p<matrix.outerSize(); ++p)
    {
        int rowCount=0;
        for (Divergence::Matrix::InnerIterator it(matrix,p); it; ++it)
        {
            check(std::isfinite(it.value()) && it.value() != 0, "nonfinite or explicit-zero coefficient");
            ++rowCount; ++actual[it.col()]; columnSum[it.col()] += it.value();
        }
        check(rowCount <= 8, "extra cell-side entries");
    }
    for (int u=0; u<layout.velocityCount(); ++u)
    {
        check(incidence[u] == 1 || incidence[u] == 2, "invalid independent graph incidence");
        check(actual[u] == (area[u] == 0 ? 0 : incidence[u]), "wrong face incidence or lost axis column");
        if (incidence[u] == 2) check(columnSum[u] == 0, "internal coefficients do not cancel exactly");
        else observe(columnSum[u]-expected[u],area[u],metricDefect);
    }
}
int coarseClosure(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, const Divergence& op,
    bool rz, const std::vector<double>& u)
{
    std::map<std::pair<int,int>,std::vector<const FaceSegment*>> groups;
    for (const auto& segment : graph.segments())
    {
        const auto& cell = mesh.cells(segment.K);
        const auto f = face(mesh,segment,rz);
        const long double length = f.horizontal ? f.x1-f.x0 : f.y1-f.y0;
        const long double whole = f.horizontal ? cell.size().x : cell.size().y;
        // Topology determines splitting. Endpoint subtraction on non-dyadic
        // bounds can round a full side below cell.size(), even without AMR.
        const bool split = segment.L >= 0 && mesh.cells(segment.L).level() > cell.level();
        const long double endpointScale = f.horizontal ? std::abs(f.x0)+std::abs(f.x1)
                                                      : std::abs(f.y0)+std::abs(f.y1);
        const long double geometryRoundoff = 64*std::numeric_limits<double>::epsilon()*std::max(whole,endpointScale);
        check(std::abs((split ? 2*length : length)-whole) <= geometryRoundoff,
            "actual endpoint length disagrees with face topology");
        if (split)
        {
            check(mesh.cells(segment.L).level() == cell.level()+1, "AMR fixture is not actual 2:1 subfaces");
            groups[{cell.index(),static_cast<int>(segment.side)}].push_back(&segment);
        }
    }
    for (const auto& group : groups)
    {
        check(group.second.size() == 2, "coarse side is not exactly two real subfaces");
        const int p = op.layout().pressureDof(group.first.first);
        long double actual=0, exact=0, scale=0;
        for (const auto* segment : group.second)
        {
            const int col = velocityDof(mesh,op.layout(),*segment);
            const auto q = sign(*segment)*face(mesh,*segment,rz).area*u[col];
            actual += static_cast<long double>(op.matrix().coeff(p,col))*u[col]; exact += q; scale += std::abs(q);
        }
        observe(actual-exact,scale,coarseDefect);
    }
    coarseSides += static_cast<int>(groups.size()); return static_cast<int>(groups.size());
}
struct Retained { Divergence op; Values values; };
Retained run(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph, bool rz, bool requireMixed)
{
    const auto layout = Layout::build(mesh,NativeLeafDomain::build(mesh,graph,rz));
    const auto op = Divergence::build(mesh,layout); const auto copy = op;
    check(op.layout().sameIdentity(layout) && &copy.matrix() == &op.matrix(), "operator lost immutable source/ownership");
    coefficients(mesh,graph,op,rz);
    std::mt19937_64 rng(194731); std::uniform_real_distribution<double> random(-1,1);
    std::vector<double> u(layout.velocityCount()), pressure(layout.pressureCount(),0);
    for (int trial=0; trial<32; ++trial)
    {
        const double scale = trial%3 == 0 ? 1e-100 : trial%3 == 1 ? 1 : 1e100;
        for (auto& value : u) value = trial == 0 ? 0 : trial == 1 ? 1 : scale*random(rng);
        const Values values(layout,u,pressure); const auto result = op.apply(mesh,values);
        const auto ref = reference(mesh,graph,layout,rz,u);
        long double total=0;
        for (int p=0; p<layout.pressureCount(); ++p)
        { observe(result[p]-ref.integral[p],ref.scale[p],localDefect); total += result[p]; }
        observe(total-ref.boundary,ref.totalScale,globalDefect);
        if (trial == 1)
        {
            const Values changedPressure(layout,u,std::vector<double>(pressure.size(),42));
            check(result == op.apply(mesh,changedPressure), "integral divergence depends on pressure");
        }
        ++randomFields;
    }
    const int mixed = coarseClosure(mesh,graph,op,rz,u);
    if (requireMixed) check(mixed > 0, "fixture missed AMR coarse/fine closure");
    for (Field field : {Field::Linear,Field::Quadratic,Field::Stream})
    {
        const auto analyticVelocity = fieldValues(mesh,graph,layout,rz,field);
        const auto ref = reference(mesh,graph,layout,rz,analyticVelocity);
        const auto result = op.apply(mesh,Values(layout,analyticVelocity,pressure));
        for (int p=0; p<layout.pressureCount(); ++p)
        {
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[layout.nativeCell(p)]);
            const auto c = cell.centerCoord(), w = cell.size();
            const long double r0=c.x-.5L*w.x, r1=c.x+.5L*w.x;
            const long double volume = rz ? pi*(r1*r1-r0*r0)*w.y : static_cast<long double>(w.x)*w.y;
            long double exact = 0;
            if (field == Field::Linear) exact = volume*(rz ? 1.85L : 1.15L);
            if (field == Field::Quadratic)
                exact = volume*(rz ? 1.4L+.3L*c.y-.25L*radialMean(r0,r1) : .7L+.6L*c.y-.25L*c.x);
            observe(result[p]-exact,ref.scale[p],field == Field::Stream ? streamDefect : analyticDefect);
        }
    }
    ++scenarios;
    std::cout << "HYDRO_DIV," << rz << ',' << mesh.cellIndexNum() << ',' << layout.pressureCount() << ','
              << layout.velocityCount() << ',' << op.matrix().nonZeros() << ',' << mixed << '\n';
    return {op,Values(layout,fieldValues(mesh,graph,layout,rz,Field::Stream),pressure)};
}
void grids()
{
    for (int metric : {0,1,2})
    {
        for (int n : {2,4,8,16,32})
        { QuadtreeMesh mesh(input(n,n,metric)); run(mesh,selected(mesh,0),metric!=0,false); }
        QuadtreeMesh rectangular(input(4,8,metric)); run(rectangular,selected(rectangular,0),metric!=0,false);
        for (int n : {4,8,16})
        {
            QuadtreeMesh mesh(input(n,n,metric)); refine(mesh,2);
            for (int selection : {0,1,2,3}) run(mesh,selected(mesh,selection),metric!=0,selection==0);
        }
        // Non-dyadic bounds exercise true endpoint measures and roundoff,
        // beyond the exactly representable grids used for topology checks.
        auto nonDyadic=input(8,6,metric);
        nonDyadic.minCoord={metric==0 ? -.37 : metric==1 ? 0.0 : .17,-.23};
        nonDyadic.maxCoord={metric==0 ? 1.19 : 1.31,1.07};
        QuadtreeMesh irregularBounds(nonDyadic); refine(irregularBounds,2);
        for (int selection : {0,3}) run(irregularBounds,selected(irregularBounds,selection),metric!=0,selection==0);
    }
}
void rebuilds()
{
    for (int metric : {0,1,2})
    {
        std::optional<Retained> previous;
        for (int depth : {0,1,2,0,1})
        {
            QuadtreeMesh mesh(input(8,8,metric)); refine(mesh,depth);
            auto current = run(mesh,selected(mesh,0),metric!=0,depth>0);
            if (previous)
            {
                rejects([&] { previous->op.apply(mesh,previous->values); });
                rejects([&] { current.op.apply(mesh,previous->values); });
            }
            previous = std::move(current);
        }
    }
}
void guardsAndLifetime()
{
    static_assert(!std::is_default_constructible_v<Divergence>);
    static_assert(std::is_same_v<decltype(std::declval<const Divergence&>().matrix()),const Divergence::Matrix&>);
    QuadtreeMesh mesh(input(4,4,1)); auto graph = selected(mesh,0);
    const auto layout = Layout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    const auto op = Divergence::build(mesh,layout);
    const Values values(layout,std::vector<double>(layout.velocityCount(),1),std::vector<double>(layout.pressureCount()));
    const auto sameSize = Layout::build(mesh,NativeLeafDomain::build(mesh,graph,false));
    rejects([&] { op.apply(mesh,Values(sameSize,values.normalVelocity(),values.pressure())); });
    DictionaryMesh stretchedInput = input(4,4,1); stretchedInput.maxCoord.x += .25;
    QuadtreeMesh stretched(stretchedInput);
    check(stretched.cellIndexNum()==mesh.cellIndexNum(), "stale same-count fixture differs in dimensions");
    rejects([&] { Divergence::build(stretched,layout); }); rejects([&] { op.apply(stretched,values); });
    auto& vertex = mesh.verts(mesh.vertIndexToLabel[0]); const auto saved = vertex.coord();
    vertex.setCoord(saved+Vector2D{.01,0}); rejects([&] { op.apply(mesh,values); });
    vertex.setCoord(saved); op.apply(mesh,values);
    const auto rz = Layout::build(mesh,NativeLeafDomain::build(mesh,graph,true));
    rejects([&] { op.apply(mesh,Values(rz,values.normalVelocity(),values.pressure())); });
    const auto left = Layout::build(mesh,NativeLeafDomain::build(mesh,selected(mesh,1),false));
    FaceSegmentGraph rightGraph; rightGraph.build(mesh,[](const QuadtreeMesh&,const Cell& c) {return c.centerCoord().x>.5;});
    const auto right = Layout::build(mesh,NativeLeafDomain::build(mesh,rightGraph,false));
    check(left.velocityCount()==right.velocityCount() && left.pressureCount()==right.pressureCount(), "selection dimensions differ");
    rejects([&] { Divergence::build(mesh,left).apply(mesh,Values(right,std::vector<double>(right.velocityCount()),std::vector<double>(right.pressureCount()))); });
    refine(mesh,1); rejects([&] { op.apply(mesh,values); });
    DictionaryMesh large = input(2,2,0); large.minCoord={0,0}; large.maxCoord={4,4};
    QuadtreeMesh wide(large); const auto wideLayout=Layout::build(wide,NativeLeafDomain::build(wide,selected(wide,0),false));
    const auto wideOp=Divergence::build(wide,wideLayout); auto huge=std::vector<double>(wideLayout.velocityCount(),0);
    huge[0]=std::numeric_limits<double>::max();
    rejects([&] { wideOp.apply(wide,Values(wideLayout,huge,std::vector<double>(wideLayout.pressureCount()))); });
    const auto zero=wideOp.apply(wide,Values(wideLayout,std::vector<double>(huge.size()),std::vector<double>(wideLayout.pressureCount())));
    check(std::all_of(zero.begin(),zero.end(),[](double value) {return value==0;}), "overflow poisoned immutable operator");
    std::optional<Retained> retained;
    {
        QuadtreeMesh temporary(input(4,4,1)); const auto d=NativeLeafDomain::build(temporary,selected(temporary,0),true);
        const auto l=Layout::build(temporary,d);
        retained.emplace(Retained{Divergence::build(temporary,l),Values(l,std::vector<double>(l.velocityCount()),std::vector<double>(l.pressureCount()))});
    }
    QuadtreeMesh replacement(input(4,4,1));
    check(retained->op.matrix().rows()==16 && retained->op.apply(replacement,retained->values)==std::vector<double>(16,0),
        "prepared operator borrowed mesh/domain lifetime");
}
}
int main() try
{
    grids(); rebuilds(); guardsAndLifetime();
    check(scenarios==75 && randomFields==2400 && coarseSides>0, "incomplete configured divergence coverage");
    std::cout << std::scientific << std::setprecision(16)
              << "H02b HYDRO divergence PASS scenarios=" << scenarios << " random=" << randomFields
              << " rejected=" << rejected << " coarse_sides=" << coarseSides
              << " metric=" << metricDefect << " local=" << localDefect << " global=" << globalDefect
              << " coarse=" << coarseDefect << " analytic=" << analyticDefect << " stream=" << streamDefect << '\n';
    return 0;
}
catch (const std::exception& error) { std::cerr << "H02b FAIL: " << error.what() << '\n'; return 1; }
