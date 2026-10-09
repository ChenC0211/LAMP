#include "Modules/Interface/HermiteGalsTemporalPhaseFlux.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace
{
using namespace Quantum::Interface;
using Velocity = HermiteGalsTemporalVelocitySource;
using Characteristic = HermiteGalsTemporalCharacteristicSource;
using Interval = HermiteGalsPhaseFluxIntervalSource;
using Transfer = HermiteGalsIntegratedFaceTransfer;
using Volumes = HermiteGalsIntegratedFaceVolumes;
using Orientation = HermiteGalsFaceTransferOrientation;
using Position = HermiteGalsTimeLayer::Position;
using Query = HermiteGalsVelocityView::Query;
using Jet = LevelSetBoundaryTrace::VelocityJet;
static_assert(!std::is_default_constructible_v<Interval> && !std::is_default_constructible_v<Transfer>);
static_assert(!std::is_default_constructible_v<Velocity> && !std::is_default_constructible_v<Characteristic>);
static_assert(std::is_same_v<decltype(std::declval<Transfer&>().integratedVolumes()), const Volumes&>);
static_assert(std::is_same_v<decltype(std::declval<Interval&>().previousAccepted()), const HermiteGalsState&>);
static_assert(std::is_same_v<decltype(std::declval<Velocity&>().previous()), const HermiteGalsVelocityView&>);
static_assert(std::is_same_v<decltype(std::declval<Characteristic&>().domain()), const LevelSetBoundaryTrace::Domain&>);
int rejections = 0;
void check(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
template<class F> void rejects(F action, const char* why)
{
    try { action(); }
    catch (const std::exception& e) {
        check(std::string(e.what()).find(why) != std::string::npos, e.what());
        ++rejections; return;
    }
    throw std::runtime_error(std::string("missing rejection: ") + why);
}
std::unique_ptr<QuadtreeMesh> meshFor(int n = 1)
{
    DictionaryMesh d; d.Nx = d.Ny = n; d.maxCellNum = 64; d.maxInitialIterations = 0;
    return std::make_unique<QuadtreeMesh>(d);
}
// Structural state fixture, not a production initializer or temporal sample.
HermiteGalsState stateFor(const QuadtreeMesh& mesh, bool radial, HermiteGalsTimeLayer time = {})
{
    Quantum::Geometry::ActiveQuarticHermitePatch patch;
    HermiteGalsPrimaryState p; HermiteGalsClosure closure;
    p.measure = radial ? HermiteGalsMeasure::AxisymmetricVolume : HermiteGalsMeasure::CartesianArea;
    patch.activeCells.assign(mesh.cellIndexNum(), 1);
    patch.cellCenters.resize(mesh.cellIndexNum());
    for (int v = 0; v < mesh.vertIndexNum(); ++v) {
        const auto x = mesh.verts(mesh.vertIndexToLabel[v]).coord();
        p.vertices.push_back({v, x.y - .4, 0, 1}); patch.nodes.push_back({x.y - .4, 0, 1, 0});
    }
    std::set<std::array<int, 2>> sides;
    for (int c = 0; c < mesh.cellIndexNum(); ++c) {
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[c]);
        const auto lo = mesh.verts(cell.cellToVert[0]).coord(), hi = mesh.verts(cell.cellToVert[2]).coord();
        const auto liquidHeight = std::clamp(.4 - lo.y, 0., hi.y - lo.y);
        p.geometryCells.push_back(c);
        p.phaseVolumes.push_back({c, liquidHeight * (radial ? std::acos(-1.) * (hi.x * hi.x - lo.x * lo.x) : hi.x - lo.x)});
        patch.cellCenters[c] = cell.centerCoord().y - .4;
        for (int k = 0; k < 4; ++k) {
            const auto& a = mesh.verts(cell.cellToVert[k]); const auto& b = mesh.verts(cell.cellToVert[(k+1)%4]);
            std::array<int, 2> key{{a.index(), b.index()}};
            if (a.coord().y == b.coord().y ? a.coord().x > b.coord().x : a.coord().y > b.coord().y) std::swap(key[0], key[1]);
            sides.insert(key);
        }
    }
    for (const auto& key : sides) {
        const auto a = mesh.verts(mesh.vertIndexToLabel[key[0]]).coord(), b = mesh.verts(mesh.vertIndexToLabel[key[1]]).coord();
        const double value = .5 * (a.y + b.y) - .4, dn = a.y == b.y ? 1. : 0.;
        p.edges.push_back({key, value, dn}); patch.edgeMidpoints.push_back({key, value, dn});
    }
    closure.vertexMixedDerivatives.assign(p.vertices.size(), 0); closure.cellBubbleCoefficients = patch.cellCenters;
    HermiteGalsCandidate candidate(mesh, time, p);
    candidate.bindGeometry(mesh, candidate.identity(), closure,
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh, patch, radial));
    return std::move(candidate).finish();
}
Query constant(double speed)
{ return [speed](const Vector2D&) { return Jet{{0, speed}, {0, 0, 0, 0}}; }; }
LevelSetBoundaryTrace::Domain domainFor(const QuadtreeMesh& mesh, bool radial)
{ return {mesh.minCoord, mesh.maxCoord, {1. / mesh.Nx, 1. / mesh.Ny}, radial}; }
Velocity velocityFor(const HermiteGalsState& old, HermiteGalsIdentity arrival, double speed = .1)
{ return {old, arrival, {1, 0, 1, Position::End}, constant(speed), constant(speed)}; }
int nonAxisVerticalFace(const QuadtreeMesh& mesh, const Interval& interval)
{
    for (int row = 0; row < mesh.vEdgeIndexNum(); ++row) if (interval.face(mesh, row).area > 0) return row;
    throw std::runtime_error("no positive-measure vertical face");
}
void contracts(bool radial, int n)
{
    auto mesh = meshFor(n); auto accepted = stateFor(*mesh, radial);
    const auto oldId = accepted.identity(); const auto* geometry = &accepted.geometry().snapshot();
    const auto oldM = accepted.primary().phaseVolumes;
    auto next = stateFor(*mesh, radial, {1, 0, 1, Position::End});
    const auto velocity = velocityFor(accepted, next.identity());
    const auto domain = domainFor(*mesh, radial); const Characteristic characteristic(*mesh, domain);
    const Interval half(*mesh, accepted, velocity, characteristic, 0, .5), full(*mesh, accepted, velocity, characteristic, 0, 1);
    const Interval second(*mesh, accepted, velocity, characteristic, .5, 1);
    auto copy = half; auto meshCopy = mesh->clone(); copy.validateMesh(*meshCopy);
    check(half.dt() == .5 && full.dt() == 1 && !half.sameInterval(full) && half.sameInterval(copy), "extent identity");
    check(half.previousGeometryIdentity() == oldId && full.previousGeometryIdentity() == oldId
        && &half.previousAccepted().geometry().snapshot() == geometry && &full.previousAccepted().geometry().snapshot() == geometry,
        "old geometry was rebuilt");
    check(half.velocitySource().sameSource(velocity) && &half.velocitySource().previous() == &velocity.previous()
        && half.characteristicSource().sameSource(characteristic) && &half.meshIdentity() == &characteristic.meshIdentity(),
        "source owner not retained");
    check(half.measure() == accepted.primary().measure
        && characteristic.method() == Characteristic::Method::CoherentExplicitMidpointRK2, "measure/method lost");
    const int face = nonAxisVerticalFace(*mesh, half);
    Volumes values{5, 2, 3, 1}; const Transfer a(*mesh, half, face, Orientation::PositiveX, values);
    values.integratedLiquidPositiveVolume = 4;
    check(a.signedLiquidVolume() == 2 && a.signedTotalVolume() == 3
        && a.averageLiquidFlux() == 4 && a.averageTotalFlux() == 6, "volume/rate or copied input contract");
    const Transfer b(*mesh, second, face, Orientation::PositiveX, {7, 4, 2, 3});
    const auto joined = a.appended(*mesh, b); const auto& sum = joined.integratedVolumes();
    check(sum.integratedTotalPositiveVolume == 12 && sum.integratedTotalNegativeVolume == 6
        && sum.integratedLiquidPositiveVolume == 5 && sum.integratedLiquidNegativeVolume == 4
        && joined.interval().sameInterval(full) && joined.averageLiquidFlux() == 1, "component-wise interval sum");
    check(a.signedLiquidVolume() == 2 && b.signedLiquidVolume() == -1, "composition mutated inputs");
    for (auto pair : {std::array<double,2>{0,0}, {1,0}, {0,std::numeric_limits<double>::infinity()},
            {std::numeric_limits<double>::quiet_NaN(),1}})
        rejects([&] { Interval bad(*mesh, accepted, velocity, characteristic, pair[0], pair[1]); }, "invalid finite positive interval");
    rejects([&] { Interval bad(*mesh, accepted, velocity, characteristic, -.1, .5); }, "outside retained velocity");
    rejects([&] { Interval bad(*mesh, accepted, velocity, characteristic, .5, 1.1); }, "outside retained velocity");
    auto foreign = stateFor(*mesh, radial);
    rejects([&] { Interval bad(*mesh, foreign, velocity, characteristic, 0, .5); }, "previous geometry/time");
    auto midpoint = stateFor(*mesh, radial, {1, 0, 1, Position::Midpoint});
    rejects([&] { Interval bad(*mesh, midpoint, velocity, characteristic, .5, 1); }, "midpoint candidate");
    auto wrongMeasure = domain; wrongMeasure.radialAxisSymmetry = !radial;
    const Characteristic wrongMetric(*mesh, wrongMeasure);
    rejects([&] { Interval bad(*mesh, accepted, velocity, wrongMetric, 0, .5); }, "measure mismatch");
    for (int k = 0; k < 4; ++k) for (long double bad : {-1.L, -std::numeric_limits<long double>::denorm_min(),
            std::numeric_limits<long double>::infinity(), std::numeric_limits<long double>::quiet_NaN()}) {
        Volumes v{5,2,3,1};
        long double* components[] = {&v.integratedTotalPositiveVolume, &v.integratedTotalNegativeVolume,
            &v.integratedLiquidPositiveVolume, &v.integratedLiquidNegativeVolume};
        *components[k] = bad;
        rejects([&] { Transfer invalid(*mesh, half, face, Orientation::PositiveX, v); }, "finite and nonnegative");
    }
    const Transfer zero(*mesh, half, face, Orientation::PositiveX, {-0.L,0,-0.L,0});
    check(zero.signedLiquidVolume() == 0, "signed zero is not zero");
    rejects([&] { Transfer bad(*mesh, half, face, Orientation::PositiveX, {5,2,6,1}); }, "exceeds total");
    rejects([&] { Transfer bad(*mesh, half, face, Orientation::PositiveX, {5,2,3,3}); }, "exceeds total");
    rejects([&] { Transfer bad(*mesh, half, -1, Orientation::PositiveX, {}); }, "native face identity");
    rejects([&] { Transfer bad(*mesh, half, face, Orientation::PositiveY, {}); }, "face orientation");
    const Transfer differentFace(*mesh, second, mesh->vEdgeIndexNum(), Orientation::PositiveY, {});
    rejects([&] { a.appended(*mesh, differentFace); }, "face identity/orientation");
    rejects([&] { b.appended(*mesh, a); }, "chronologically contiguous");
    rejects([&] { a.appended(*mesh, a); }, "chronologically contiguous");
    const Interval gap(*mesh, accepted, velocity, characteristic, .75, 1);
    const Transfer gapTransfer(*mesh, gap, face, Orientation::PositiveX, {});
    rejects([&] { a.appended(*mesh, gapTransfer); }, "chronologically contiguous");
    const auto otherVelocity = velocityFor(accepted, next.identity(), .2);
    check(otherVelocity.previous().geometryIdentity() == velocity.previous().geometryIdentity()
        && otherVelocity.arrival().geometryIdentity() == velocity.arrival().geometryIdentity()
        && !otherVelocity.sameSource(velocity), "same geometry/time conflated different velocity providers");
    const Interval otherV(*mesh, accepted, otherVelocity, characteristic, .5, 1);
    const Transfer otherVT(*mesh, otherV, face, Orientation::PositiveX, {});
    rejects([&] { a.appended(*mesh, otherVT); }, "different velocity source");
    const Characteristic otherCharacteristic(*mesh, domain);
    const Interval otherC(*mesh, accepted, velocity, otherCharacteristic, .5, 1);
    const Transfer otherCT(*mesh, otherC, face, Orientation::PositiveX, {});
    rejects([&] { a.appended(*mesh, otherCT); }, "different characteristic source");
    const auto foreignVelocity = velocityFor(foreign, next.identity());
    const Interval foreignI(*mesh, foreign, foreignVelocity, characteristic, .5, 1);
    const Transfer foreignT(*mesh, foreignI, face, Orientation::PositiveX, {});
    rejects([&] { a.appended(*mesh, foreignT); }, "different previous geometry");
    if (radial) for (int row = 0; row < mesh->vEdgeIndexNum(); ++row) if (half.face(*mesh,row).area == 0) {
        const Transfer axis(*mesh, half, row, Orientation::PositiveX, {});
        check(axis.signedTotalVolume() == 0, "axis zero volume");
        rejects([&] { Transfer bad(*mesh, half, row, Orientation::PositiveX, {1,0,0,0}); }, "zero-measure face");
    }
    if (n > 1) {
        int internal = 0;
        for (int row = 0; row < mesh->vEdgeIndexNum() + mesh->hEdgeIndexNum(); ++row) {
            const auto f = half.face(*mesh, row); if (f.negativeCell < 0 || f.positiveCell < 0) continue;
            const Transfer t(*mesh, half, row, f.horizontal ? Orientation::PositiveY : Orientation::PositiveX, {5,2,3,1});
            check(t.face().negativeCell == f.negativeCell && t.face().positiveCell == f.positiveCell,
                "shared actual face owners lost"); ++internal;
        }
        check(internal == 4, "missing uniform internal face coverage");
    }
    // Mutating even one coordinate invalidates the exact retained revision.
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({-.125, 0});
    rejects([&] { half.validateMesh(*mesh); }, "mesh");
    rejects([&] { a.appended(*mesh, b); }, "mesh");
    mesh->verts(mesh->vertIndexToLabel[0]).setCoord({0, 0});
    half.validateMesh(*mesh);
    check(accepted.identity() == oldId && &accepted.geometry().snapshot() == geometry, "failure changed old authority");
    for (std::size_t k = 0; k < oldM.size(); ++k) check(accepted.primary().phaseVolumes[k].liquid == oldM[k].liquid, "failure changed M");
    std::cout << "INTERVAL_CONTRACT," << (radial ? "RZ" : "XY") << ',' << n << ",PASS\n";
}
void sourceLifetimeAndGuards()
{
    auto mesh = meshFor(); auto state = stateFor(*mesh, false); auto next = stateFor(*mesh, false, {1,0,1,Position::End});
    std::optional<Interval> retained; std::weak_ptr<const Jet> weak;
    int calls = 0;
    {
        auto field = std::make_shared<const Jet>(Jet{{0,.125},{0,0,0,0}}); weak = field;
        Query query = [field,&calls](const Vector2D&) { ++calls; return *field; };
        Velocity velocity(state, next.identity(), {1,0,1,Position::End}, query, query);
        Characteristic characteristic(*mesh, domainFor(*mesh, false));
        retained.emplace(*mesh, state, velocity, characteristic, 0, .5);
        auto movedV = std::move(velocity); auto movedC = std::move(characteristic);
        rejects([&] { velocity.previous(); }, "moved-from velocity");
        rejects([&] { characteristic.domain(); }, "moved-from characteristic");
        check(retained->velocitySource().sameSource(movedV) && retained->characteristicSource().sameSource(movedC), "move lost retained identity");
    }
    check(!weak.expired() && calls == 0, "contract evaluated query or lost callback ownership");
    const auto traced = traceHermiteGalsWithJacobian(retained->characteristicSource().domain(), {.5,.5},
        retained->velocitySource().previous(), retained->velocitySource().arrival());
    check(std::abs(traced.departurePoint().y - .375) < 1e-14 && calls > 0, "retained actual tracer/velocity source unusable");
    retained.reset(); check(weak.expired(), "callback ownership leaked");
    rejects([&] { Velocity bad(state, next.identity(), {1,0,1,Position::End}, {}, constant(0)); }, "missing paired");
    rejects([&] { Velocity bad(state, {}, {1,0,1,Position::End}, constant(0), constant(0)); }, "unbound geometry");
    rejects([&] { Velocity bad(state, next.identity(), {1,.1,1,Position::End}, constant(0), constant(0)); }, "compatible previous/arrival");
    auto domain = domainFor(*mesh, false); domain.minCoord.x = -.1;
    rejects([&] { Characteristic bad(*mesh, domain); }, "mesh bounds");
    domain = domainFor(*mesh, false); domain.cellSize.x = 0;
    rejects([&] { Characteristic bad(*mesh, domain); }, "cell size");
    const auto velocity = velocityFor(state, next.identity()); const Characteristic characteristic(*mesh, domainFor(*mesh, false));
    Interval half(*mesh, state, velocity, characteristic, 0, .5), second(*mesh, state, velocity, characteristic, .5, 1);
    const int face = nonAxisVerticalFace(*mesh, half); const auto maximum = std::numeric_limits<long double>::max();
    const Transfer a(*mesh, half, face, Orientation::PositiveX, {maximum,0,maximum,0});
    const Transfer b(*mesh, second, face, Orientation::PositiveX, {maximum,0,maximum,0});
    rejects([&] { a.appended(*mesh, b); }, "finite and nonnegative");
    rejects([&] { a.averageLiquidFlux(); }, "not representable");
    check(a.signedLiquidVolume() == maximum, "failed derived rate changed integrated authority");
    auto differentMesh = meshFor(2); auto differentState = stateFor(*differentMesh, false);
    auto differentEnd = stateFor(*differentMesh, false, {1,0,1,Position::End});
    auto differentVelocity = velocityFor(differentState, differentEnd.identity());
    const Characteristic differentCharacteristic(*differentMesh, domainFor(*differentMesh, false));
    const Interval differentInterval(*differentMesh, differentState, differentVelocity, differentCharacteristic, .5, 1);
    const Transfer differentTransfer(*differentMesh, differentInterval, 0, Orientation::PositiveX, {});
    rejects([&] { a.appended(*differentMesh, differentTransfer); }, "mesh");
    rejects([&] { a.appended(*mesh, differentTransfer); }, "mesh");
    rejects([&] { Interval bad(*mesh, state, velocity, characteristic,
        -std::numeric_limits<double>::max(), std::numeric_limits<double>::max()); }, "finite positive interval");
    const Interval quarter(*mesh, state, velocity, characteristic, 0, .25);
    const Interval rest(*mesh, state, velocity, characteristic, .25, 1);
    const Transfer q(*mesh, quarter, face, Orientation::PositiveX, {5,2,3,1});
    const Transfer r(*mesh, rest, face, Orientation::PositiveX, {7,4,2,3});
    const auto unequal = q.appended(*mesh, r);
    check(unequal.signedLiquidVolume() == 1 && unequal.averageLiquidFlux() == 1
        && unequal.signedTotalVolume() == 6, "unequal intervals averaged rates instead of adding volumes");
    auto moved = std::move(half); rejects([&] { half.dt(); }, "moved-from interval");
    Transfer t(*mesh, moved, face, Orientation::PositiveX, {}); auto retainedT = std::move(t);
    rejects([&] { t.signedLiquidVolume(); }, "moved-from integrated");
    check(retainedT.interval().dt() == .5, "moved transfer lost interval");
    std::cout << "SOURCE_LIFETIME,NO_SAMPLING,FAILURE_ISOLATION,PASS\n";
}
void authorityLifetime()
{
    std::optional<Interval> saved; HermiteGalsIdentity id;
    {
        auto mesh = meshFor(); const auto old = stateFor(*mesh, false, {0,3,3,Position::Initial});
        const auto end = stateFor(*mesh, false, {1,3,4,Position::End}); id = old.identity();
        const Velocity velocity(old, end.identity(), end.timeLayer(), constant(.1), constant(.1));
        const Characteristic characteristic(*mesh, domainFor(*mesh, false));
        saved.emplace(*mesh, old, velocity, characteristic, 3.25, 3.75);
    }
    check(saved->previousGeometryIdentity() == id && saved->timeBegin() == 3.25 && saved->dt() == .5
        && saved->previousAccepted().primary().phaseVolumes[0].liquid == .4
        && std::abs(saved->previousAccepted().geometry().snapshot().scalarJetNative({.5,.5}).value - .1) < 1e-14,
        "interval failed to retain authority after inputs and mesh died");
    auto equivalentMesh = meshFor(); saved->validateMesh(*equivalentMesh);
    std::cout << "AUTHORITY_LIFETIME,RESTART_SUBINTERVAL,PASS\n";
}
}
int main() try
{
    for (bool radial : {false,true}) for (int n : {1,2}) contracts(radial,n);
    sourceLifetimeAndGuards();
    authorityLifetime();
    std::cout << "I12a1 temporal interval contract PASS rejections=" << rejections << "\n";
    return 0;
}
catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
