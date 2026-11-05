#include "Modules/Interface/HermiteGalsState.h"

#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace
{
using namespace LAMP::Interface;
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Patch = LAMP::Geometry::ActiveQuarticHermitePatch;
using Position = HermiteGalsTimeLayer::Position;

template<class T, class = void> struct hasMixed : std::false_type {};
template<class T> struct hasMixed<T, std::void_t<decltype(std::declval<T>().dxy)>> : std::true_type {};
static_assert(!hasMixed<HermiteGalsVertexJet>::value, "a first jet cannot carry transported dxy");
static_assert(!std::is_copy_constructible_v<HermiteGalsCandidate>, "mutable candidates cannot share storage");
static_assert(!std::is_default_constructible_v<HermiteGalsState>, "no empty completed state");
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsState&>().primary()),
    const HermiteGalsPrimaryState&>, "completed primary must be read-only");
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsCandidate&>().primary()),
    const HermiteGalsPrimaryState&>, "edits must go through identity invalidation");
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsState&>().closure()),
    const HermiteGalsClosure&>, "closure is derived and read-only");
static_assert(std::is_same_v<decltype(std::declval<HermiteGalsGeometryView&>().snapshot()),
    const Snapshot&>, "geometry has no writable view");
static_assert(!std::is_constructible_v<HermiteGalsGeometryView, std::shared_ptr<Snapshot>>,
    "external shared_ptr must not create a mutable alias");

void check(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
template<class F> void rejects(F action, const char* fragment)
{
    try { action(); }
    catch (const std::exception& error) {
        check(std::string(error.what()).find(fragment) != std::string::npos,
            "rejection did not identify the expected failure");
        return;
    }
    throw std::runtime_error(std::string("missing rejection: ") + fragment);
}
DictionaryMesh meshInput(double shift = 0.0)
{
    DictionaryMesh input;
    input.minCoord = {shift, 0}; input.maxCoord = {shift + 1, 1};
    input.Nx = input.Ny = 1; input.maxCellNum = 16; input.maxInitialIterations = 0;
    return input;
}

// Analytic ownership fixture on ONE unit cell: phi=x-(xmin+1/2), M=1/2.
// This is not a V2 initializer, closure solver or transport benchmark.
struct Fixture
{
    HermiteGalsPrimaryState primary;
    HermiteGalsClosure closure;
    Patch patch;
    explicit Fixture(const QuadtreeMesh& mesh)
    {
        const double middle = mesh.minCoord.x + .5;
        patch.activeCells.assign(1, 1);
        patch.cellCenters.assign(1, 0.0);
        for (int v = 0; v < mesh.vertIndexNum(); ++v) {
            const double phi = mesh.verts(mesh.vertIndexToLabel[v]).coord().x - middle;
            primary.vertices.push_back({v, phi, 1, 0});
            patch.nodes.push_back({phi, 1, 0, 0});
        }
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[0]);
        for (int side = 0; side < 4; ++side) {
            const auto& a = mesh.verts(cell.cellToVert[side]);
            const auto& b = mesh.verts(cell.cellToVert[(side + 1) % 4]);
            const bool horizontal = a.coord().y == b.coord().y;
            std::array<int, 2> ids{{a.index(), b.index()}};
            if (horizontal ? a.coord().x > b.coord().x : a.coord().y > b.coord().y)
                std::swap(ids[0], ids[1]);
            const double phi = (a.coord().x + b.coord().x) * .5 - middle;
            primary.edges.push_back({ids, phi, horizontal ? 0.0 : 1.0});
            patch.edgeMidpoints.push_back({ids, phi, horizontal ? 0.0 : 1.0});
        }
        primary.geometryCells = {0}; primary.phaseVolumes = {{0, .5}};
        closure.vertexMixedDerivatives.assign(primary.vertices.size(), 0.0);
        closure.cellBubbleCoefficients = {0.0};
    }
    Snapshot geometry(const QuadtreeMesh& mesh) const { return Snapshot::fromNativeLeaves(mesh, patch); }
};

void ownershipAndLifetime()
{
    std::optional<HermiteGalsGeometryView> saved;
    HermiteGalsIdentity initialIdentity;
    {
        QuadtreeMesh mesh(meshInput()); Fixture f(mesh);
        HermiteGalsCandidate initial(mesh, {}, f.primary);
        initialIdentity = initial.identity();
        check(initialIdentity.value() != 0 && !initial.hasGeometry(), "initial candidate readiness");
        f.primary.vertices[0].phi = 999; f.primary.phaseVolumes[0].liquid = .25;
        check(initial.primary().vertices[0].phi != 999 && initial.primary().phaseVolumes[0].liquid == .5,
            "caller mutation changed primary storage");
        rejects([&] { std::move(initial).finish(); }, "has not been bound");
        auto source = f.geometry(mesh);
        initial.bindGeometry(mesh, initialIdentity, f.closure, std::move(source));
        f.closure.vertexMixedDerivatives[0] = 91;
        auto accepted = std::move(initial).finish(); // test App explicitly chooses this handle
        rejects([&] { initial.identity(); }, "consumed handle");
        rejects([&] { source.nativeCellCount(); }, "moved");
        check(accepted.closure().vertexMixedDerivatives[0] == 0, "closure aliases its caller");
        accepted.validateMesh(mesh);
        auto sameMesh = mesh.clone(); accepted.validateMesh(*sameMesh);
        saved = accepted.geometry();
        check(saved->identity() == accepted.identity() && saved->timeLayer().time() == 0,
            "geometry lost its state identity or time");
        check(&accepted.geometry().snapshot() == &saved->snapshot(), "second snapshot authority created");

        const HermiteGalsTimeLayer midpoint{1, 0, .2, Position::Midpoint};
        {
            HermiteGalsCandidate rejected(mesh, midpoint, accepted.primary());
            auto changed = rejected.primary(); changed.phaseVolumes[0].liquid = .25;
            rejected.replacePrimary(mesh, changed);
            check(rejected.primary().phaseVolumes[0].liquid == .25, "M was recomputed from geometry");
        }
        check(accepted.identity() == initialIdentity && accepted.primary().phaseVolumes[0].liquid == .5,
            "rejected candidate polluted accepted state");

        // A completed module candidate remains separate until the App replaces accepted.
        Fixture fresh(mesh); HermiteGalsCandidate candidate(mesh, midpoint, fresh.primary);
        candidate.bindGeometry(mesh, candidate.identity(), fresh.closure, fresh.geometry(mesh));
        auto completed = std::move(candidate).finish();
        check(accepted.timeLayer().time() == 0 && completed.timeLayer().time() == .1,
            "finishing an Interface candidate advanced the App clock");
        auto copied = accepted;
        accepted = completed;
        check(copied.identity() == initialIdentity && saved->identity() == initialIdentity,
            "replacing App handle changed retained old authority");
        auto moved = std::move(completed);
        rejects([&] { completed.geometry(); }, "consumed handle");
        check(moved.timeLayer().time() == .1, "moved state lost time");
        mesh.verts(mesh.vertIndexToLabel[0]).setCoord({-.125, 0});
        rejects([&] { copied.validateMesh(mesh); }, "mesh binding");
    }
    check(saved->identity() == initialIdentity &&
        std::abs(saved->snapshot().scalarJetInNativeCell(0, {.25, .75}).value + .25) < 1e-14,
        "geometry view did not retain state/mesh lifetime");
    auto movedView = std::move(*saved);
    rejects([&] { saved->snapshot(); }, "consumed handle");
    check(movedView.identity() == initialIdentity, "moved geometry view lost authority");
    std::cout << "ownership: input copies, rejection isolation, explicit App replacement and retained lifetime PASS\n";
}

void revisionAndBinding()
{
    QuadtreeMesh mesh(meshInput()); Fixture f(mesh);
    HermiteGalsCandidate a(mesh, {}, f.primary), b(mesh, {}, f.primary);
    check(a.identity() != b.identity(), "equal time/data reused candidate identity");
    auto geometry = f.geometry(mesh);
    rejects([&] { a.bindGeometry(mesh, b.identity(), f.closure, std::move(geometry)); }, "foreign primary");
    check(geometry.nativeCellCount() == 1 && !a.hasGeometry(), "failed binding consumed input");
    const auto old = a.identity();
    a.bindGeometry(mesh, old, f.closure, std::move(geometry));
    auto another = f.geometry(mesh);
    rejects([&] { a.bindGeometry(mesh, old, f.closure, std::move(another)); }, "unique snapshot");
    auto invalid = a.primary(); invalid.phaseVolumes[0].liquid = -1;
    rejects([&] { a.replacePrimary(mesh, invalid); }, "negative phase volume");
    check(a.identity() == old && a.hasGeometry(), "failed replacement damaged candidate");
    a.replacePrimary(mesh, a.primary()); // also safe when input aliases this candidate
    check(a.identity() != old && !a.hasGeometry(), "primary replacement retained stale geometry");
    rejects([&] { a.bindGeometry(mesh, old, f.closure, std::move(another)); }, "stale/foreign");
    rejects([&] { std::move(a).finish(); }, "has not been bound");

    auto badClosure = f.closure; badClosure.cellBubbleCoefficients.clear();
    rejects([&] { a.bindGeometry(mesh, a.identity(), badClosure, std::move(another)); }, "bubble closure size");
    badClosure = f.closure; badClosure.vertexMixedDerivatives.pop_back();
    rejects([&] { a.bindGeometry(mesh, a.identity(), badClosure, std::move(another)); }, "mixed closure size");
    badClosure = f.closure; badClosure.cellBubbleCoefficients[0] = std::numeric_limits<double>::infinity();
    rejects([&] { a.bindGeometry(mesh, a.identity(), badClosure, std::move(another)); }, "non-finite bubble");
    badClosure = f.closure; badClosure.vertexMixedDerivatives[0] = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { a.bindGeometry(mesh, a.identity(), badClosure, std::move(another)); }, "non-finite mixed");

    LAMP::Geometry::ActiveHermitePatch cubicPatch;
    cubicPatch.activeCells = f.patch.activeCells; cubicPatch.nodes = f.patch.nodes;
    auto cubic = Snapshot::fromNativeLeaves(mesh, cubicPatch);
    rejects([&] { a.bindGeometry(mesh, a.identity(), f.closure, std::move(cubic)); }, "enriched-quartic");
    QuadtreeMesh foreignMesh(meshInput(2)); Fixture foreign(foreignMesh);
    auto foreignGeometry = foreign.geometry(foreignMesh);
    rejects([&] { a.bindGeometry(foreignMesh, a.identity(), f.closure, std::move(foreignGeometry)); }, "mesh binding");
    rejects([&] { a.bindGeometry(mesh, a.identity(), f.closure, std::move(foreignGeometry)); }, "bounds mismatch");
    auto inactive = f.primary; inactive.geometryCells.clear();
    HermiteGalsCandidate wrongActivity(mesh, {}, inactive);
    auto inactiveClosure = f.closure; inactiveClosure.cellBubbleCoefficients.clear();
    rejects([&] { wrongActivity.bindGeometry(mesh, wrongActivity.identity(), inactiveClosure, std::move(another)); }, "activity mismatch");
    a.bindGeometry(mesh, a.identity(), f.closure, std::move(another));
    auto moved = std::move(a);
    rejects([&] { a.hasGeometry(); }, "consumed handle");
    check(std::move(moved).finish().primary().phaseVolumes[0].liquid == .5, "recovery after failed binding");
    QuadtreeMesh decimalMesh(meshInput(.1)); Fixture decimal(decimalMesh);
    HermiteGalsCandidate decimalCandidate(decimalMesh, {}, decimal.primary);
    decimalCandidate.bindGeometry(decimalMesh, decimalCandidate.identity(), decimal.closure,
        decimal.geometry(decimalMesh));
    std::move(decimalCandidate).finish().validateMesh(decimalMesh);
    std::cout << "binding: revision invalidation, stale/foreign identity, unique snapshot and failure atomicity PASS\n";
}

void inputValidation()
{
    QuadtreeMesh mesh(meshInput()); Fixture f(mesh);
    const auto invalidPrimary = [&](auto mutate, const char* message) {
        auto primary = f.primary; mutate(primary);
        rejects([&] { HermiteGalsCandidate candidate(mesh, {}, primary); }, message);
    };
    invalidPrimary([](auto& p) { p.vertices[0].vertex = -1; }, "master vertex indices");
    invalidPrimary([](auto& p) { p.vertices[1].vertex = p.vertices[0].vertex; }, "master vertex indices");
    invalidPrimary([](auto& p) { p.vertices.back().vertex = 999; }, "master vertex indices");
    invalidPrimary([](auto& p) { p.vertices[0].dx = std::numeric_limits<double>::quiet_NaN(); }, "non-finite master vertex");
    invalidPrimary([](auto& p) { p.edges.push_back(p.edges[0]); }, "duplicate master edge");
    invalidPrimary([](auto& p) { std::swap(p.edges[0].vertices[0], p.edges[0].vertices[1]); }, "positive Cartesian tangent");
    invalidPrimary([](auto& p) { p.edges[0].vertices[0] = -1; }, "not a native vertex");
    invalidPrimary([](auto& p) { p.edges[0].dn = std::numeric_limits<double>::infinity(); }, "non-finite master edge");
    invalidPrimary([](auto& p) { p.geometryCells = {0, 0}; }, "geometry cell indices");
    invalidPrimary([](auto& p) { p.geometryCells = {1}; }, "geometry cell indices");
    invalidPrimary([](auto& p) { p.phaseVolumes.clear(); }, "lacks conserved inventory");
    invalidPrimary([](auto& p) { p.phaseVolumes.push_back(p.phaseVolumes[0]); }, "inventory cell indices");
    invalidPrimary([](auto& p) { p.phaseVolumes[0].liquid = -std::numeric_limits<double>::denorm_min(); }, "negative phase volume");
    invalidPrimary([](auto& p) { p.phaseVolumes[0].liquid = std::numeric_limits<double>::quiet_NaN(); }, "non-finite or negative");
    invalidPrimary([](auto& p) { p.measure = static_cast<HermiteGalsMeasure>(99); }, "unknown measure");
    for (auto time : {HermiteGalsTimeLayer{-1, 0, 0, Position::Initial},
            {0, 0, 1, Position::Initial}, {0, 0, 1, Position::End},
            {1, 1, 0, Position::End}, {1, 0, 0, Position::Midpoint},
            {1, 0, std::numeric_limits<double>::infinity(), Position::End},
            {1, 0, 1, static_cast<Position>(99)}})
        rejects([&] { HermiteGalsCandidate candidate(mesh, time, f.primary); }, "time:");
    check(HermiteGalsTimeLayer{2, 1, 3, Position::Midpoint}.time() == 2 &&
        HermiteGalsTimeLayer{2, 1, 3, Position::End}.time() == 3 &&
        HermiteGalsTimeLayer{12, 9, 9, Position::Initial}.time() == 9, "time layer/restart semantics");
    auto rz = f.primary; rz.measure = HermiteGalsMeasure::AxisymmetricVolume;
    rz.phaseVolumes[0].liquid = std::acos(-1.0) * .25; // true revolved volume r<1/2
    HermiteGalsCandidate radial(mesh, {}, rz);
    check(radial.primary().phaseVolumes[0].liquid == rz.phaseVolumes[0].liquid,
        "RZ inventory was silently converted into planar area");
    mesh.verts(mesh.vertIndexToLabel[0]).setType(TVertType::west); // explicit hanging-role fixture
    rejects([&] { HermiteGalsCandidate candidate(mesh, {}, f.primary); }, "hanging vertex");
    std::cout << "input: native indices, finite data, unclipped inventory, explicit measure and time layers PASS\n";
}
}

int main()
{
    try { ownershipAndLifetime(); revisionAndBinding(); inputValidation(); }
    catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
    std::cout << "I01 state/API contract PASS (no V2 numerical publication/transport claim)\n";
    return 0;
}
