#include "Core/Geometry/ActiveLevelSetTransport.h"
#include "Core/Geometry/ActiveQuarticHermiteReconstruction.h"
#include "Core/Geometry/ExactPointKey.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Field/VertField.h"
#include "Core/Geometry/LevelSetCurvature.h"
#include "Core/Geometry/LevelSetQueryReconstruction.h"
#include "Core/Geometry/LevelSetPhaseFraction.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <exception>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace ActiveLevelSetTransport
{
namespace
{
void require(bool ok, const char *reason)
{
    if (!ok)
        throw std::runtime_error(std::string("[ActiveLevelSetTransport] ") + reason);
}
bool same(Vector2D a, Vector2D b)
{
    return a.x == b.x && a.y == b.y;
}
double extension(double phi, const Options &options)
{
    if (options.physicalVelocity) return 1.0;
    const double a = std::abs(phi);
    if (a <= options.fullWidth)
        return 1.0;
    if (a >= options.zeroWidth)
        return 0.0;
    // Reuse the shared tail-factored quintic. Direct evaluation near q=1
    // can round above one; neither clipping nor a looser guard is needed.
    const double halfWidth = 0.5 * (options.zeroWidth - options.fullWidth);
    const double centered = 0.5 * ((options.zeroWidth - a) - (a - options.fullWidth));
    const double chi = LevelSetPhaseFraction::regularizedPositiveIndicator(centered, halfWidth);
    require(std::isfinite(chi) && chi >= 0.0 && chi <= 1.0, "invalid C2 velocity extension");
    return chi;
}

// One retained map for scalar transport and material endpoints. The inverse
// below reuses getNextCoord's bounded fixed-point algorithm, applied to the
// actual fitted-velocity backtrace, not to an interpolated vertex field.
class Characteristic
{
  public:
    Characteristic(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
                   const LevelSetBoundaryTrace::Domain &domain, double dt,
                   LevelSetBoundaryTrace::Velocity current, LevelSetBoundaryTrace::Velocity previous,
                   const Options &options)
        : source_(std::move(source)), domain_(domain), dt_(dt), current_(std::move(current)),
          previous_(std::move(previous))
    {
        widths_.fullWidth = options.fullWidth;
        widths_.zeroWidth = options.zeroWidth;
        widths_.physicalVelocity = options.physicalVelocity;
        inflow_=options.inflowScalar;
        bool selected=source_->restrictedDomain();
        if (source_->usesNativeLeaves()) {
            selected=false;
            for (std::size_t c=0;c<source_->nativeCellCount();++c)
                if (!source_->nativeCellActive(static_cast<int>(c))) {selected=true;break;}
        }
        if (selected)
            domain_.activeSegment=[source=source_](Vector2D a,Vector2D b){return source->segmentInsideActiveDomain(a,b);};
    }

    const LevelSetHermiteGeometry::Snapshot &source() const
    {
        return *source_;
    }

    Vector2D backward(Vector2D point) const
    {
        const auto sample=backwardSample(point);
        require(!sample.inflowValue,"material map cannot cross a prescribed inflow");
        return sample.point;
    }

    struct BackwardSample { Vector2D point; double sourceValue; std::optional<double> inflowValue; };
    BackwardSample backwardSample(Vector2D point) const
    {
        // One RK2 backtrace revisits its starting point and queries the same
        // midpoint with both velocity time levels. Only this call owns the
        // two scalar samples; velocity callbacks are still evaluated in their
        // original order. Exact coordinate bits, no near-point snapping.
        std::array<LAMP::Geometry::ExactPointKey, 2> keys;
        std::array<double, 2> values{};
        std::size_t count = 0;
        const auto scalarAt = [&](Vector2D p) {
            // A velocity callback can change the current floating mode.
            // Retained point validity is not permission to ignore that change.
            LAMP::Geometry::Detail::ExactExpansion::validateEnvironment();
            const auto key = LAMP::Geometry::exactPointKey(p);
            for (std::size_t index = 0; index < count; ++index)
                if (keys[index] == key) return values[index];
            const double value = source_->value(p);
            if (count < keys.size())
            {
                keys[count] = key;
                values[count++] = value;
            }
            return value;
        };
        // Even dt==0 and the zero-speed tail require an active scalar point.
        // This never makes a solid point legal or invents an inflow value.
        const double sourceValue=scalarAt(point);
        const LevelSetBoundaryTrace::Velocity current = [this, &scalarAt](const Vector2D &p) {
            return extended(current_, p, scalarAt(p));
        };
        const LevelSetBoundaryTrace::Velocity previous = [this, &scalarAt](const Vector2D &p) {
            return extended(previous_, p, scalarAt(p));
        };
        const auto trace = LevelSetBoundaryTrace::trace(domain_, point, dt_, current, previous);
        require(!trace.boundaryHit || inflow_,
                "transport extension reached a physical inflow; explicit scalar closure is required");
        // Activity coverage is not phase sampling. In particular, a tiny
        // valid trace need not contain a representable phase sampling point.
        require(same(point, trace.point) || source_->segmentInsideActiveDomain(point, trace.point),
                "backtrace crosses inactive scalar cells");
        if (trace.boundaryHit) {
            const double value=inflow_(*trace.side,trace.point,trace.elapsed);
            require(std::isfinite(value),"nonfinite prescribed inflow scalar");
            return {trace.point,sourceValue,value};
        }
        return {trace.point,sourceValue,std::nullopt};
    }

    Vector2D forward(Vector2D point) const
    {
        Vector2D next = point;
        Vector2D traced = backward(next);
        const double h = std::min(domain_.cellSize.x, domain_.cellSize.y);
        const double scale = std::max({h, std::abs(domain_.minCoord.x), std::abs(domain_.minCoord.y),
                                       std::abs(domain_.maxCoord.x), std::abs(domain_.maxCoord.y)});
        const double tolerance = std::max(128.0 * std::numeric_limits<double>::epsilon() * scale, 1e-12 * h);
        require(std::isfinite(tolerance) && tolerance > 0.0 && tolerance <= 1e-8 * h,
                "forward trace coordinate scale is unresolved");
        Vector2D previousNext = next, previousTraced = traced;
        for (int iteration = 0; iteration < 32; ++iteration)
        {
            if (iteration != 0)
                traced = backward(next);
            const Vector2D residual = point - traced;
            if (std::hypot(residual.x, residual.y) <= tolerance)
                return next;
            // Solve B(next)=point using the same B as nodal advection. No
            // independent forward RK2 formula, projection or wall clipping.
            const Vector2D candidate = next + residual;
            require(std::isfinite(candidate.x) && std::isfinite(candidate.y) &&
                        std::hypot(candidate.x - point.x, candidate.y - point.y) <= 0.5 * h,
                    "forward trace exceeds half a cell displacement");
            require(!same(candidate, next), "forward trace inverse stagnated above tolerance");
            if (iteration == 31)
            {
                // Failure-only diagnostics: do not reevaluate either velocity
                // callback, change the inverse map, or relax its acceptance.
                std::ostringstream message;
                message
                    << std::setprecision(17)
                    << "[ActiveLevelSetTransport] forward trace inverse did not converge in 32 evaluations"
                    << "; point=(" << point.x << ',' << point.y << ')' << "; next=(" << next.x << ','
                    << next.y << ')' << "; backward(next)=(" << traced.x << ',' << traced.y << ')'
                    << "; residual=(" << residual.x << ',' << residual.y << ')' << "; previous_next=("
                    << previousNext.x << ',' << previousNext.y << ')' << "; previous_backward=("
                    << previousTraced.x << ',' << previousTraced.y << ')' << "; tolerance=" << tolerance
                    << "; dt=" << dt_ << "; h=" << h;
                throw std::runtime_error(message.str());
            }
            previousNext = next;
            previousTraced = traced;
            next = candidate;
        }
        throw std::runtime_error(
            "[ActiveLevelSetTransport] forward trace inverse did not converge in 32 evaluations");
    }

  private:
    Vector2D extended(const LevelSetBoundaryTrace::Velocity &velocity, Vector2D point, double scalar) const
    {
        const double chi = extension(scalar, widths_);
        if (chi == 0.0)
            return Vector2D(0.0, 0.0);
        const auto u = velocity(point);
        require(std::isfinite(u.x) && std::isfinite(u.y), "nonfinite physical velocity sample");
        return u * chi;
    }

    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
    LevelSetBoundaryTrace::Domain domain_;
    double dt_;
    LevelSetBoundaryTrace::Velocity current_, previous_;
    Options widths_; // Only the two scalar-amplitude widths, no nodal masks.
    InflowScalar inflow_;
};
struct Masks
{
    std::vector<std::uint8_t> active, fixed;
};
template<class Patch>
Masks prepareMasks(const QuadtreeMesh &mesh, const LevelSetHermiteGeometry::Snapshot &source,
                   const Patch &previous, const Options &options)
{
    std::vector<std::uint8_t> active(mesh.vertIndexNum(), 0), fixed(mesh.vertIndexNum(), 0);
    for (int k = 0; k < mesh.cellIndexNum(); ++k)
        if (source.usesNativeLeaves() ? source.nativeCellActive(k) : previous.activeCells[k])
            for (int label : mesh.cells(mesh.cellIndexToLabel[k]).cellToVert)
                active[mesh.verts(label).index()] = 1;
    require(options.oneSidedSamples.empty() || options.oneSidedSamples.size() == active.size(),
            "one-sided sample mask size mismatch");
    require(options.extensionSideSamples.empty() || options.extensionSideSamples.size() == active.size(),
            "extension-side sample mask size mismatch");
    for (std::size_t k = 0; k < options.oneSidedSamples.size(); ++k)
        require(options.oneSidedSamples[k] <= 1 && (!options.oneSidedSamples[k] || active[k]),
                "one-sided derivative support must be an active binary subset");
    for (std::size_t k = 0; k < options.extensionSideSamples.size(); ++k)
        require(options.extensionSideSamples[k] <= 1 && (!options.extensionSideSamples[k] || active[k]),
                "extension derivative support must be an active binary subset");
    for (int k : options.fixedZeroVertices)
    {
        require(k >= 0 && k < mesh.vertIndexNum() && active[k] && !fixed[k] &&
                    (!source.usesNativeLeaves() || mesh.verts(mesh.vertIndexToLabel[k]).supCell < 0) &&
                    previous.nodes[k].value == 0.0,
                "fixed zero must identify a unique active zero DOF");
        fixed[k] = 1;
    }
    for (const auto &oldSector : previous.cornerSectors)
    {
        bool retained = false;
        for (const auto &extension : options.fixedWallCornerExtensions)
            retained |= extension.cellIndex == oldSector.cellIndex &&
                        extension.vertexIndex == oldSector.vertexIndex &&
                        extension.phaseSign == oldSector.phaseSign;
        require(retained, "a retained corner sector cannot be silently dropped during transport");
    }
    return {std::move(active), std::move(fixed)};
}
void validateInputs(double dt, const LevelSetBoundaryTrace::Velocity &current,
                    const LevelSetBoundaryTrace::Velocity &previous, const Options &options)
{
    require(std::isfinite(dt) && dt >= 0.0 && current && previous,
            "finite dt and both velocity callbacks required");
    require(options.physicalVelocity || (std::isfinite(options.fullWidth) && options.fullWidth > 0.0 && std::isfinite(options.zeroWidth) &&
                options.zeroWidth > options.fullWidth),
            "explicit positive nested scalar-amplitude transport widths required");
}
void validatePin(Vector2D point, const LevelSetBoundaryTrace::Velocity &current,
                 const LevelSetBoundaryTrace::Velocity &previous)
{
    const auto a = current(point), b = previous(point);
    require(a.x == 0.0 && a.y == 0.0 && b.x == 0.0 && b.y == 0.0,
            "this pinned zero requires both actual no-slip velocity traces to vanish");
}
template<class MapJet>
void pullbackNodalJet(LAMP::Geometry::ActiveHermitePatch::NodalDof& node,
    const LevelSetHermiteGeometry::Snapshot::ScalarJet2& jet,const MapJet& x,const MapJet& y)
{
    const double xx=1+x.dx,xy=x.dy,yx=y.dx,yy=1+y.dy;
    node.dx=jet.dx*xx+jet.dy*yx;
    node.dy=jet.dx*xy+jet.dy*yy;
    node.dxy=jet.dxx*xx*xy+jet.dxy*(xx*yy+yx*xy)+jet.dyy*yx*yy+jet.dx*x.dxy+jet.dy*y.dxy;
    require(std::isfinite(node.dx) && std::isfinite(node.dy) && std::isfinite(node.dxy),"nonfinite pulled-back nodal jet");
}

template<class CandidateType>
void finishCandidate(const QuadtreeMesh &mesh, const Options &options, const Masks &masks,
                     const std::shared_ptr<const Characteristic> &characteristic, CandidateType &result)
{
    result.patch.cornerSectors.clear();
    for (const auto &extension : options.fixedWallCornerExtensions)
    {
        require(extension.vertexIndex >= 0 && extension.vertexIndex < mesh.vertIndexNum() &&
                    masks.fixed[extension.vertexIndex],
                "corner extension requires a declared fixed zero vertex");
        LAMP::Geometry::ActiveHermitePatch::CornerSector sector;
        sector.cellIndex = extension.cellIndex;
        sector.vertexIndex = extension.vertexIndex;
        sector.phaseSign = extension.phaseSign;
        sector.dof = result.patch.nodes[extension.vertexIndex];
        if (extension.wallHorizontal)
            sector.dof.dx = 0.0;
        else
            sector.dof.dy = 0.0;
        result.patch.cornerSectors.push_back(sector);
    }
    try
    {
        if (characteristic->source().usesNativeLeaves())
            result.geometry = std::make_unique<LevelSetHermiteGeometry::Snapshot>(
                LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh, result.patch,
                                                                    options.radialAxisSymmetry));
        else if constexpr (std::is_same_v<typename CandidateType::Patch, LAMP::Geometry::ActiveHermitePatch>)
            result.geometry = std::make_unique<LevelSetHermiteGeometry::Snapshot>(
                mesh, result.patch, options.radialAxisSymmetry, CartesianGridScope::UniformLeaves);
        else
            throw std::runtime_error("[ActiveLevelSetTransport] quartic candidate requires native leaves");
    }
    catch (const std::exception &error)
    {
        std::ostringstream message;
        message << std::setprecision(17)
                << "[ActiveLevelSetTransport] rejected candidate geometry: changedValues="
                << result.changedValues << ", maxValueIncrement=" << result.maxValueIncrement
                << ", maxTraceRatio=" << result.maxTraceRatio;
        for (int k : options.fixedZeroVertices)
        {
            const auto &node = result.patch.nodes[k];
            message << "; fixed vertex=" << k << " DOF=(" << node.value << ',' << node.dx << ',' << node.dy
                    << ',' << node.dxy << ')';
        }
        throw std::runtime_error(message.str() + ": " + error.what());
    }
    result.forwardMaterialPoint = [characteristic](Vector2D point) { return characteristic->forward(point); };
}

struct NativeSpacing
{
    std::vector<Vector2D> vertices;
    Vector2D minimum;
};
NativeSpacing nativeSpacing(const QuadtreeMesh &mesh, const LevelSetHermiteGeometry::Snapshot &source)
{
    const double infinity = std::numeric_limits<double>::infinity();
    NativeSpacing result{std::vector<Vector2D>(mesh.vertIndexNum(), {infinity, infinity}),
                         {infinity, infinity}};
    for (int cell = 0; cell < mesh.cellIndexNum(); ++cell)
        if (source.nativeCellActive(cell))
        {
            const auto bounds = source.nativeCellBounds(cell);
            const Vector2D h = bounds[1] - bounds[0];
            result.minimum.x = std::min(result.minimum.x, h.x);
            result.minimum.y = std::min(result.minimum.y, h.y);
            for (int label : mesh.cells(mesh.cellIndexToLabel[cell]).cellToVert)
            {
                auto &local = result.vertices[mesh.verts(label).index()];
                local.x = std::min(local.x, h.x);
                local.y = std::min(local.y, h.y);
            }
        }
    return result;
}

// Local query stencils on the retained scalar, NOT a second uniform mesh.
// Each sample transports the same scalar increment by the same characteristic.
// Geometry alone selects the nearest admitted 5x5 box containing the vertex;
// all positive-area leaf intersections must belong to the requested chart.
class NativeStencil
{
  public:
    NativeStencil(const QuadtreeMesh &mesh, const Characteristic &characteristic, bool radial)
        : mesh_(mesh), characteristic_(characteristic), reconstruction_(radial)
    {
    }

    LevelSetPolynomial::Derivatives derivatives(Vector2D point, Vector2D spacing,
                                                const std::vector<std::uint8_t> &mask)
    {
        // This point belongs to every admitted stencil. Keep one public
        // active-scalar / strict-FP check per reconstruction even if ALL
        // of its sample increments have already been evaluated.
        characteristic_.source().value(point);
        long double variation = 0;
        return reconstruction_.evaluateDerivatives(point, spacing,
            [&](Vector2D lo, Vector2D hi) { return admitted(lo, hi, mask); },
            [&](Vector2D sample) { return increment(sample); }, variation);
    }

    std::size_t sampleRequests() const { return sampleRequests_; }
    std::size_t sampleEvaluations() const { return increments_.size(); }

    template<class Sample>
    LevelSetPolynomial::Derivatives derivativesOf(Vector2D point,Vector2D spacing,
        const std::vector<std::uint8_t>& mask,const Sample& sample) const
    {
        characteristic_.source().value(point);
        long double variation=0;
        return reconstruction_.evaluateDerivatives(point,spacing,
            [&](Vector2D lo,Vector2D hi){return admitted(lo,hi,mask);},sample,variation);
    }

  private:
    const QuadtreeMesh &mesh_;
    const Characteristic &characteristic_;
    LevelSetPolynomial::QueryReconstruction reconstruction_;
    // Owned by one serial advanceNative pass with one immutable scalar,
    // dt, extension and pair of physical velocity fields. Chart masks select
    // query locations, but do not redefine the scalar at those locations.
    // Integer keys distinguish every double bit pattern, including +/- zero;
    // no rounding, tolerance, cell ID or reconstructed lattice is used.
    using SampleKey = std::array<std::uint64_t, 2>;
    std::map<SampleKey, double> increments_;
    std::size_t sampleRequests_ = 0;
    struct AdmissionKey {
        LAMP::Geometry::ExactPointKey lower,upper;
        const std::vector<std::uint8_t>* mask;
        bool operator==(const AdmissionKey& other) const noexcept
        {return lower==other.lower && upper==other.upper && mask==other.mask;}
    };
    struct AdmissionHash {
        std::size_t operator()(const AdmissionKey& key) const noexcept {
            const LAMP::Geometry::ExactPointKeyHash hash;
            return hash({static_cast<std::uint64_t>(hash(key.lower)),
                static_cast<std::uint64_t>(hash(key.upper)^std::hash<const void*>{}(key.mask))});
        }
    };
    // One serial advance owns this immutable source and fixed chart masks.
    // x/y pullbacks (and neighbouring nodes) revisit identical support boxes.
    // Cache only exact coordinate bits + mask identity, never sampled values
    // or velocity callbacks. The public scalar/FP check above still runs on
    // every derivative query. Failed geometric queries are never cached.
    mutable std::unordered_map<AdmissionKey,bool,AdmissionHash> admissions_;

    double increment(Vector2D sample)
    {
        static_assert(sizeof(double) == sizeof(std::uint64_t), "sample key requires binary64 storage");
        SampleKey key;
        std::memcpy(&key[0], &sample.x, sizeof(double));
        std::memcpy(&key[1], &sample.y, sizeof(double));
        ++sampleRequests_;
        const auto found = increments_.find(key);
        if (found != increments_.end())
            return found->second;
        const auto traced = characteristic_.backward(sample);
        const double delta = same(traced, sample) ? 0 : characteristic_.source().value(traced) -
            characteristic_.source().value(sample);
        require(std::isfinite(delta), "nonfinite native Q4 scalar increment");
        // Failed trajectory / scalar queries are never cached. A later
        // derivative or candidate failure destroys this private map too.
        increments_.emplace(key, delta);
        return delta;
    }

    bool admitted(Vector2D lo, Vector2D hi, const std::vector<std::uint8_t> &mask) const
    {
        const auto &source = characteristic_.source();
        const auto lower = source.minCoord(), upper = source.maxCoord();
        if (!std::isfinite(lo.x) || !std::isfinite(lo.y) || !std::isfinite(hi.x) || !std::isfinite(hi.y) ||
            lo.x < lower.x || lo.y < lower.y || hi.x > upper.x || hi.y > upper.y ||
            !(lo.x < hi.x && lo.y < hi.y))
            return false;
        const AdmissionKey key{LAMP::Geometry::exactPointKey(lo),LAMP::Geometry::exactPointKey(hi),&mask};
        const auto cached=admissions_.find(key);
        if(cached!=admissions_.end())return cached->second;
        const auto evaluate=[&]() {
        bool covered = false;
        for (int cell : source.nativeCellCandidates(lo, hi))
        {
            const auto bounds = source.nativeCellBounds(cell);
            if (std::max(lo.x, bounds[0].x) >= std::min(hi.x, bounds[1].x) ||
                std::max(lo.y, bounds[0].y) >= std::min(hi.y, bounds[1].y))
                continue;
            if (!source.nativeCellActive(cell))
                return false;
            if (!mask.empty())
                for (int label : mesh_.cells(mesh_.cellIndexToLabel[cell]).cellToVert)
                    if (!mask[mesh_.verts(label).index()])
                        return false;
            covered = true;
        }
        return covered;
        };
        const bool result=evaluate();admissions_.emplace(key,result);return result;
    }
};

} // namespace

Candidate advance(const QuadtreeMesh &mesh, const LAMP::Geometry::ActiveHermitePatch &previous, double dt,
                  const LevelSetBoundaryTrace::Velocity &newVelocity,
                  const LevelSetBoundaryTrace::Velocity &oldVelocity, const Options &options)
{
    const UniformLeafGrid grid(mesh, CartesianGridScope::UniformLeaves);
    validateInputs(dt, newVelocity, oldVelocity, options);
    require(!options.independentVelocityQueries, "concurrent velocity queries require native quartic transport");
    require(!options.inflowScalar,"scalar inflow requires the native jet transport port");
    // Construct from the supplied patch itself; a separately retained scalar
    // cannot silently supply different DOFs or activity to this characteristic.
    const auto characteristic = std::make_shared<const Characteristic>(
        std::make_shared<const LevelSetHermiteGeometry::Snapshot>(mesh, previous, options.radialAxisSymmetry,
                                                                  CartesianGridScope::UniformLeaves),
        LevelSetBoundaryTrace::Domain{grid.minCoord(), grid.maxCoord(), grid.spacing(),
                                      options.radialAxisSymmetry},
        dt, newVelocity, oldVelocity, options);
    const auto &source = characteristic->source();
    const auto masks = prepareMasks(mesh, source, previous, options);
    const auto &active = masks.active;
    const auto &fixed = masks.fixed;
    Candidate result;
    result.patch = previous;
    VertField<double> increment(mesh, 0.0);
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        if (!active[k])
            continue;
        ++result.activeNodes;
        const auto point = mesh.verts(mesh.vertIndexToLabel[k]).coord();
        try
        {
            if (fixed[k])
            {
                validatePin(point, newVelocity, oldVelocity);
                ++result.stationaryNodes;
                continue; // Value only; increment derivatives below can change the angle.
            }
            const double chi = extension(source.value(point), options);
            if (chi > 0.0 && chi < 1.0)
                ++result.taperedNodes;
            const auto traced = characteristic->backward(point);
            if (same(traced, point))
            {
                ++result.stationaryNodes;
                continue;
            }
            // Native coordinates and canonical DOF coordinates may differ by
            // roundoff. Transport their scalar INCREMENT; do not move either
            // coordinate or change phi merely by sampling a stationary node.
            const double delta = source.value(traced) - source.value(point);
            const double next = previous.nodes[k].value + delta;
            require(std::isfinite(delta) && std::isfinite(next), "nonfinite advected nodal scalar");
            increment.value[k] = delta;
            result.patch.nodes[k].value = next;
            result.changedValues += next != previous.nodes[k].value;
            result.maxValueIncrement = std::max(result.maxValueIncrement, std::abs(delta));
            result.maxTraceRatio =
                std::max(result.maxTraceRatio, std::hypot(traced.x - point.x, traced.y - point.y) /
                                                   std::min(grid.spacing().x, grid.spacing().y));
        }
        catch (const std::exception &error)
        {
            std::ostringstream message;
            message << std::setprecision(17) << "[ActiveLevelSetTransport] vertex=" << k << ", point=("
                    << point.x << ',' << point.y << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    }
    using Reconstruction = LevelSetCurvature::Reconstruction<QuadtreeMesh>;
    const Reconstruction derivatives(mesh, increment, 4, LevelSetCurvature::Basis::TensorProduct,
                                     options.radialAxisSymmetry, CartesianGridScope::UniformLeaves, active);
    std::unique_ptr<Reconstruction> oneSided;
    if (!options.oneSidedSamples.empty())
        oneSided = std::make_unique<Reconstruction>(
            mesh, increment, 4, LevelSetCurvature::Basis::TensorProduct, options.radialAxisSymmetry,
            CartesianGridScope::UniformLeaves, options.oneSidedSamples);
    std::unique_ptr<Reconstruction> extensionSide;
    if (!options.extensionSideSamples.empty())
        extensionSide = std::make_unique<Reconstruction>(
            mesh, increment, 4, LevelSetCurvature::Basis::TensorProduct, options.radialAxisSymmetry,
            CartesianGridScope::UniformLeaves, options.extensionSideSamples);
    for (int j = 0; j <= grid.ny(); ++j)
        for (int i = 0; i <= grid.nx(); ++i)
        {
            const int k =
                mesh.verts(grid.vertices()[static_cast<std::size_t>(j) * (grid.nx() + 1) + i]).index();
            if (!active[k])
                continue;
            const Vector2D point(
                i == grid.nx() ? grid.maxCoord().x : grid.minCoord().x + i * grid.spacing().x,
                j == grid.ny() ? grid.maxCoord().y : grid.minCoord().y + j * grid.spacing().y);
            const auto &reconstruction = oneSided && options.oneSidedSamples[k]             ? *oneSided
                                         : extensionSide && options.extensionSideSamples[k] ? *extensionSide
                                                                                            : derivatives;
            const auto jet = reconstruction.evaluateDerivatives(point);
            auto &node = result.patch.nodes[k];
            node.dx += jet.dx;
            node.dy += jet.dy;
            node.dxy += jet.dxy;
            require(std::isfinite(node.dx) && std::isfinite(node.dy) && std::isfinite(node.dxy),
                    "nonfinite derivative increment");
        }
    finishCandidate(mesh, options, masks, characteristic, result);
    return result;
}

Candidate advanceNative(const QuadtreeMesh &mesh, const LAMP::Geometry::ActiveHermitePatch &previous,
                        double dt, const LevelSetBoundaryTrace::Velocity &newVelocity,
                        const LevelSetBoundaryTrace::Velocity &oldVelocity, const Options &options)
{
    validateInputs(dt, newVelocity, oldVelocity, options);
    require(!options.independentVelocityQueries, "concurrent velocity queries require native quartic transport");
    require(!options.inflowScalar,"scalar inflow requires the native jet transport port");
    const auto source = std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh, previous, options.radialAxisSymmetry));
    const auto masks = prepareMasks(mesh, *source, previous, options);
    const auto layout = nativeSpacing(mesh, *source);
    const auto &spacing = layout.vertices;
    const auto minimum = layout.minimum;
    // First native entry deliberately retains a conservative global finest-cell
    // trajectory bound. Local Q4 spacing is separate; no new CFL relaxation.
    const auto characteristic = std::make_shared<const Characteristic>(
        source,
        LevelSetBoundaryTrace::Domain{source->minCoord(), source->maxCoord(), minimum,
                                      options.radialAxisSymmetry},
        dt, newVelocity, oldVelocity, options);
    NativeStencil stencil(mesh, *characteristic, options.radialAxisSymmetry);
    Candidate result;
    result.patch = previous;
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const auto &vertex = mesh.verts(mesh.vertIndexToLabel[k]);
        if (!masks.active[k] || vertex.supCell >= 0)
            continue;
        ++result.activeNodes;
        const auto point = vertex.coord();
        try
        {
            if (masks.fixed[k])
            {
                validatePin(point, newVelocity, oldVelocity);
                ++result.stationaryNodes;
            }
            else
            {
                const double chi = extension(source->value(point), options);
                result.taperedNodes += chi > 0 && chi < 1;
                const auto traced = characteristic->backward(point);
                if (same(traced, point))
                    ++result.stationaryNodes;
                else
                {
                    const double delta = source->value(traced) - source->value(point);
                    auto &node = result.patch.nodes[k];
                    node.value += delta;
                    require(std::isfinite(delta) && std::isfinite(node.value),
                            "nonfinite advected nodal scalar");
                    result.changedValues += node.value != previous.nodes[k].value;
                    result.maxValueIncrement = std::max(result.maxValueIncrement, std::abs(delta));
                    result.maxTraceRatio =
                        std::max(result.maxTraceRatio, std::hypot(traced.x - point.x, traced.y - point.y) /
                                                           std::min(minimum.x, minimum.y));
                }
            }
            const auto &support = !options.oneSidedSamples.empty() && options.oneSidedSamples[k]
                                      ? options.oneSidedSamples
                                  : !options.extensionSideSamples.empty() && options.extensionSideSamples[k]
                                      ? options.extensionSideSamples
                                      : masks.active;
            const auto jet = stencil.derivatives(point, spacing[k], support);
            auto &node = result.patch.nodes[k];
            node.dx += jet.dx;
            node.dy += jet.dy;
            node.dxy += jet.dxy;
            require(std::isfinite(node.dx) && std::isfinite(node.dy) && std::isfinite(node.dxy),
                    "nonfinite derivative increment");
        }
        catch (const std::exception &error)
        {
            std::ostringstream message;
            message << std::setprecision(17) << "[ActiveLevelSetTransport] native vertex=" << k << ", point=("
                    << point.x << ',' << point.y << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    }
    result.nativeStencilSampleRequests = stencil.sampleRequests();
    result.nativeStencilSampleEvaluations = stencil.sampleEvaluations();
    finishCandidate(mesh, options, masks, characteristic, result);
    return result;
}

Candidate advanceNativeJet(const QuadtreeMesh& mesh,const LAMP::Geometry::ActiveHermitePatch& previous,
    double dt,const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity,const Options& options)
{
    validateInputs(dt,newVelocity,oldVelocity,options);
    require(!options.independentVelocityQueries,"concurrent velocity queries require native quartic transport");
    const auto source=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,previous,options.radialAxisSymmetry));
    const auto masks=prepareMasks(mesh,*source,previous,options);
    // A pinned wet corner has two scalar jets. The base jet belongs to the
    // cap cell across the seam, not to the deterministic point-query owner
    // (which may be the pure wet extension). Retain explicit old sectors;
    // introducing a new contact representation is a separate operation.
    require(previous.cornerSectors.size()==options.fixedWallCornerExtensions.size(),
        "jet pullback requires the retained wet-corner representation");
    std::vector<int> capCell(mesh.vertIndexNum(),-1);
    for (const auto& extension:options.fixedWallCornerExtensions) {
        const int v=extension.vertexIndex,c=extension.cellIndex;
        require(v>=0 && v<mesh.vertIndexNum() && masks.fixed[v]
            && !options.oneSidedSamples.empty() && options.oneSidedSamples[v]
            && !options.extensionSideSamples.empty() && options.extensionSideSamples[v],
            "wet-corner jet requires a fixed zero and explicit two-sided seam supports");
        const auto old=std::find_if(previous.cornerSectors.begin(),previous.cornerSectors.end(),[&](const auto& s){
            return s.vertexIndex==v && s.cellIndex==c && s.phaseSign==extension.phaseSign;
        });
        require(old!=previous.cornerSectors.end(),"wet-corner jet changed its retained sector identity");
        const auto& base=previous.nodes[v];
        require(old->dof.value==base.value && old->dof.dxy==base.dxy &&
            (extension.wallHorizontal?(old->dof.dx==0 && old->dof.dy==base.dy):
                (old->dof.dy==0 && old->dof.dx==base.dx)),
            "retained corner is not the declared flat wall-tangent extension");
        const auto point=mesh.verts(mesh.vertIndexToLabel[v]).coord();
        const auto wet=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();
        const int seam=extension.wallHorizontal?0:1;
        const auto coordinate=[](Vector2D p,int axis){return axis?p.y:p.x;};
        const double origin=coordinate(point,seam);const bool wetBelow=coordinate(wet,seam)<origin;
        for (int n=0;n<mesh.vertIndexNum();++n) if (masks.active[n]) {
            const double x=coordinate(mesh.verts(mesh.vertIndexToLabel[n]).coord(),seam);
            const bool onSeam=x==origin,wetSide=(x<origin)==wetBelow;
            require(static_cast<bool>(options.oneSidedSamples[n])==(onSeam || !wetSide)
                && static_cast<bool>(options.extensionSideSamples[n])==(onSeam || wetSide),
                "wet-corner supports do not represent the actual cap/wet half charts");
        }
        for (int candidate:source->nativeCellCandidates(point,point)) {
            if (!source->nativeCellActive(candidate)) continue;
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[candidate]);const auto centre=cell.centerCoord();
            bool incident=false;
            for (int label:cell.cellToVert) incident|=mesh.verts(label).index()==v;
            if (!incident) continue;
            if ((coordinate(centre,seam)<coordinate(point,seam))==(coordinate(wet,seam)<coordinate(point,seam))
                || (coordinate(centre,1-seam)<coordinate(point,1-seam))!=(coordinate(wet,1-seam)<coordinate(point,1-seam))) continue;
            bool admitted=true;
            for (int label:cell.cellToVert) admitted&=options.oneSidedSamples[mesh.verts(label).index()]!=0;
            require(admitted && capCell[v]<0,"wet-corner cap cell is absent from its unique one-sided chart");
            capCell[v]=candidate;
        }
        require(capCell[v]>=0,"wet-corner jet has no actual cap-side cell");
    }
    const auto layout=nativeSpacing(mesh,*source);
    const auto characteristic=std::make_shared<const Characteristic>(source,
        LevelSetBoundaryTrace::Domain{source->minCoord(),source->maxCoord(),layout.minimum,options.radialAxisSymmetry},
        dt,newVelocity,oldVelocity,options);
    // Map components are not even scalar fields at the radial axis.
    NativeStencil stencil(mesh,*characteristic,false);
    struct Sample { Vector2D departure,displacement; double value; bool inflow; };
    std::map<LAMP::Geometry::ExactPointKey,Sample> samples;
    Candidate result;result.patch=previous;
    const auto sample=[&](Vector2D point) -> const Sample& {
        ++result.nativeStencilSampleRequests;
        const auto key=LAMP::Geometry::exactPointKey(point);
        const auto found=samples.find(key);if (found!=samples.end()) return found->second;
        const auto traced=characteristic->backwardSample(point);
        const double value=traced.inflowValue?*traced.inflowValue:
            same(point,traced.point)?traced.sourceValue:source->value(traced.point);
        require(std::isfinite(value),"nonfinite transported cubic scalar");
        result.maxTraceRatio=std::max(result.maxTraceRatio,std::hypot(traced.point.x-point.x,traced.point.y-point.y)
            /std::min(layout.minimum.x,layout.minimum.y));
        return samples.emplace(key,Sample{traced.point,traced.point-point,value,static_cast<bool>(traced.inflowValue)}).first->second;
    };
    for (int k=0;k<mesh.vertIndexNum();++k) {
        const auto& vertex=mesh.verts(mesh.vertIndexToLabel[k]);
        if (!masks.active[k] || vertex.supCell>=0) continue;
        ++result.activeNodes;const auto point=vertex.coord();
        if (masks.fixed[k]) validatePin(point,newVelocity,oldVelocity);
        const auto& traced=sample(point);auto& node=result.patch.nodes[k];
        if (masks.fixed[k]) require(same(traced.departure,point) && !traced.inflow && traced.value==0,
            "fixed zero acquired a displaced or inflow characteristic");
        node.value=traced.value;
        const auto& support=!options.oneSidedSamples.empty() && options.oneSidedSamples[k]?options.oneSidedSamples:
            !options.extensionSideSamples.empty() && options.extensionSideSamples[k]?options.extensionSideSamples:masks.active;
        bool inflow=traced.inflow;
        const auto x=stencil.derivativesOf(point,layout.vertices[k],support,[&](Vector2D p) {
            const auto& s=sample(p);inflow|=s.inflow;return s.displacement.x;
        });
        const auto y=stencil.derivativesOf(point,layout.vertices[k],support,[&](Vector2D p) {
            const auto& s=sample(p);inflow|=s.inflow;return s.displacement.y;
        });
        if (inflow) {
            const auto fit=stencil.derivativesOf(point,layout.vertices[k],support,[&](Vector2D p){return sample(p).value;});
            node.dx=fit.dx;node.dy=fit.dy;node.dxy=fit.dxy;
        } else if (!same(traced.departure,point) || x.dx!=0 || x.dy!=0 || x.dxy!=0 || y.dx!=0 || y.dy!=0 || y.dxy!=0) {
            const auto jet=capCell[k]>=0?source->scalarJetInNativeCell(capCell[k],traced.departure):
                source->scalarJetNative(traced.departure);
            pullbackNodalJet(node,jet,x,y);
        }
        if (options.radialAxisSymmetry && point.x==0) node.dx=node.dxy=0;
        require(std::isfinite(node.dx) && std::isfinite(node.dy) && std::isfinite(node.dxy),"nonfinite pulled-back cubic jet");
        const bool changed=node.value!=previous.nodes[k].value;
        result.changedValues+=changed;result.stationaryNodes+=!changed;
        const double chi=extension(source->value(point),options);
        result.taperedNodes+=chi>0 && chi<1;
        result.maxValueIncrement=std::max(result.maxValueIncrement,std::abs(node.value-previous.nodes[k].value));
    }
    result.nativeStencilSampleEvaluations=samples.size();
    finishCandidate(mesh,options,masks,characteristic,result);
    return result;
}

static QuarticCandidate advanceQuarticImpl(const QuadtreeMesh &mesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch &previous, double dt,
    const LevelSetBoundaryTrace::Velocity &newVelocity,
    const LevelSetBoundaryTrace::Velocity &oldVelocity, const Options &options, bool pullbackJet)
{
    using Clock=std::chrono::steady_clock;
    auto stageStarted=Clock::now();
    const auto elapsed=[&] {
        const auto now=Clock::now();
        const double seconds=std::chrono::duration<double>(now-stageStarted).count();
        stageStarted=now;
        return seconds;
    };
    validateInputs(dt, newVelocity, oldVelocity, options);
    require(!pullbackJet || (previous.cornerSectors.empty() && options.fixedZeroVertices.empty()
        && options.fixedWallCornerExtensions.empty()), "jet pullback does not prescribe contact data");
    require(pullbackJet || !options.inflowScalar,"scalar inflow requires the native jet transport port");
    const auto source = std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh, previous, options.radialAxisSymmetry));
    const auto masks = prepareMasks(mesh, *source, previous, options);
    const auto layout = nativeSpacing(mesh, *source);
    const auto characteristic = std::make_shared<const Characteristic>(source,
        LevelSetBoundaryTrace::Domain{source->minCoord(), source->maxCoord(), layout.minimum,
                                      options.radialAxisSymmetry},
        dt, newVelocity, oldVelocity, options);
    for (int k : options.fixedZeroVertices)
        validatePin(mesh.verts(mesh.vertIndexToLabel[k]).coord(), newVelocity, oldVelocity);

    QuarticCandidate result;
    result.timings.sourceGeometry=elapsed();
    result.patch = previous;
    ActiveQuarticHermiteReconstruction::Options reconstruction;
    reconstruction.radialAxisSymmetry = options.radialAxisSymmetry && !pullbackJet;
    reconstruction.oneSidedSamples = options.oneSidedSamples;
    reconstruction.extensionSideSamples = options.extensionSideSamples;
    const auto prepared = ActiveQuarticHermiteReconstruction::prepare(mesh, previous.activeCells,
        reconstruction, options.quarticStencils.get(), &result.reconstruction);
    const auto& points = prepared.samplePoints();
    result.timings.stencilPreparation=elapsed();
    struct Increment { double value, traceRatio; Vector2D departure,displacement; double transportedValue; bool inflow; };
    const auto evaluateIncrement = [&](Vector2D point) {
        const auto tracedSample = characteristic->backwardSample(point);
        const auto traced = tracedSample.point;
        const double ratio = std::hypot(traced.x - point.x, traced.y - point.y) /
            std::min(layout.minimum.x, layout.minimum.y);
        // The SAME backtrace has already read and validated the old scalar at
        // its starting point. Reuse that value, not a second physical query.
        const double transported=tracedSample.inflowValue?*tracedSample.inflowValue
            :same(traced,point)?tracedSample.sourceValue:source->value(traced);
        const double delta = transported-tracedSample.sourceValue;
        require(std::isfinite(delta), "nonfinite supplied scalar");
        return Increment{delta, ratio, traced, traced-point,transported,static_cast<bool>(tracedSample.inflowValue)};
    };
    // Only this call owns current values. Workers do not mutate the exact-key
    // lookup or diagnostics; those are populated below in first-use order.
    const auto sampleBatch = [&](const std::vector<Vector2D>& queries) {
        std::vector<Increment> sampled(queries.size());
        if (!options.independentVelocityQueries)
        {
            for (std::size_t i=0; i<queries.size(); ++i) sampled[i] = evaluateIncrement(queries[i]);
        }
        else
        {
            std::vector<std::exception_ptr> errors(queries.size());
            if (!queries.empty()) sampled[0] = evaluateIncrement(queries[0]);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) if(queries.size()>64)
#endif
            for (std::ptrdiff_t i=1; i<static_cast<std::ptrdiff_t>(queries.size()); ++i)
            {
                try {sampled[i] = evaluateIncrement(queries[i]);}
                catch (...) {errors[i] = std::current_exception();}
            }
            for (const auto& error:errors) if (error) std::rethrow_exception(error);
        }
        return sampled;
    };
    const auto sampled=sampleBatch(points);
    std::vector<double> values(points.size());
    // The generator and the pointwise edge/center updates use ONE repeatable
    // scalar increment. Reuse bitwise-identical coordinates only during this
    // advance; neither these values nor the callbacks enter the stencil cache
    // or the retained material trajectory. Never snap near-equal coordinates.
    std::unordered_map<LAMP::Geometry::ExactPointKey, Increment, LAMP::Geometry::ExactPointKeyHash> increments;
    increments.reserve(static_cast<std::size_t>(mesh.vertIndexNum())*4);
    for (std::size_t i=0; i<points.size(); ++i)
    {
        values[i] = sampled[i].value;
        increments.emplace(LAMP::Geometry::exactPointKey(points[i]), sampled[i]);
        result.maxTraceRatio = std::max(result.maxTraceRatio, sampled[i].traceRatio);
    }
    result.scalarIncrementRequests = result.scalarIncrementEvaluations = points.size();
    const auto scalarIncrement = [&](Vector2D point) {
        ++result.scalarIncrementRequests;
        const auto key = LAMP::Geometry::exactPointKey(point);
        const auto found = increments.find(key);
        require(found != increments.end(), "quartic point update lost its prepared scalar increment");
        return found->second.value;
    };
    result.timings.scalarSampling=elapsed();
    LAMP::Geometry::ActiveQuarticHermitePatch increment,mapX,mapY;
    std::optional<LAMP::Geometry::ActiveQuarticHermitePatch> boundaryFit;
    const auto boundaryData=[&]() -> const LAMP::Geometry::ActiveQuarticHermitePatch& {
        if (!boundaryFit) {
            std::vector<double> transported(points.size());
            for (std::size_t i=0;i<points.size();++i) transported[i]=sampled[i].transportedValue;
            boundaryFit=ActiveQuarticHermiteReconstruction::reconstruct(prepared,transported);
        }
        return *boundaryFit;
    };
    const bool hasInflowSamples=std::any_of(sampled.begin(),sampled.end(),[](const auto& x){return x.inflow;});
    ActiveQuarticHermiteReconstruction::Prepared::AffectedData affected;
    if (hasInflowSamples) {
        std::vector<std::uint8_t> marked(points.size());
        for (std::size_t i=0;i<points.size();++i) marked[i]=sampled[i].inflow;
        affected=prepared.affectedBy(marked);
        boundaryData();
    }
    if (pullbackJet) {
        std::vector<double> displacementX(points.size()),displacementY(points.size());
        for (std::size_t i=0;i<points.size();++i) {
            displacementX[i]=sampled[i].displacement.x;
            displacementY[i]=sampled[i].displacement.y;
        }
        mapX=ActiveQuarticHermiteReconstruction::reconstruct(prepared,displacementX);
        mapY=ActiveQuarticHermiteReconstruction::reconstruct(prepared,displacementY);
        // Retain the generator's topology-owned independent edge identities.
        increment=mapX;
    } else increment=ActiveQuarticHermiteReconstruction::reconstruct(prepared, values);
    result.timings.dataReconstruction=elapsed();
    const auto add = [](double &old, double delta) {
        require(std::isfinite(old) && std::isfinite(delta), "nonfinite independent quartic increment");
        if (delta != 0.0) old += delta; // Preserve stored data, including signed zero, in a stationary field.
        require(std::isfinite(old), "nonfinite advanced quartic datum");
    };
    for (int k = 0; k < mesh.vertIndexNum(); ++k)
    {
        const auto &vertex = mesh.verts(mesh.vertIndexToLabel[k]);
        if (!masks.active[k] || vertex.supCell >= 0) continue;
        ++result.activeNodes;
        auto &node = result.patch.nodes[k];
        const auto &delta = increment.nodes[k];
        if (!masks.fixed[k])
        {
            const double chi = extension(source->value(vertex.coord()), options);
            result.taperedNodes += chi > 0.0 && chi < 1.0;
            add(node.value, pullbackJet?scalarIncrement(vertex.coord()):delta.value);
        }
        if (pullbackJet) {
            const auto& traced=increments.at(LAMP::Geometry::exactPointKey(vertex.coord()));
            const auto& x=mapX.nodes[k]; const auto& y=mapY.nodes[k];
            if (hasInflowSamples && affected.nodes[k]) {
                const auto& fit=boundaryData().nodes[k];
                node.dx=fit.dx;node.dy=fit.dy;node.dxy=fit.dxy;
            } else if (traced.traceRatio!=0 || x.dx!=0 || x.dy!=0 || x.dxy!=0 || y.dx!=0 || y.dy!=0 || y.dxy!=0) {
                const auto jet=source->scalarJetNative(traced.departure);
                pullbackNodalJet(node,jet,x,y);
            }
            if (options.radialAxisSymmetry && vertex.coord().x==0) node.dx=node.dxy=0;
        } else {
            add(node.dx, delta.dx);
            add(node.dy, delta.dy);
            add(node.dxy, delta.dxy);
        }
        const bool changed = node.value != previous.nodes[k].value;
        result.changedValues += changed;
        result.stationaryNodes += !changed;
        result.maxValueIncrement = std::max(result.maxValueIncrement, std::abs(node.value - previous.nodes[k].value));
    }
    // A valid supplied patch may omit dependent fine sides or reorder full
    // sides. Match explicit identities, never equate vector positions from
    // the generator with the retained scalar's independent column indices.
    std::map<std::array<int, 2>, std::size_t> retainedEdges;
    for (std::size_t i = 0; i < previous.edgeMidpoints.size(); ++i)
        retainedEdges.emplace(previous.edgeMidpoints[i].vertices, i);
    // Q4's complete support uses integer lattice points, whereas independent
    // values also require actual edge midpoints and native cell centers. Batch
    // their missing queries in the ORIGINAL edge/center order. Exact-key
    // admission and publication remain serial; callbacks use the same explicit
    // concurrency contract as the stencil batch. No new physical sample.
    std::vector<Vector2D> additionalPoints;
    std::unordered_set<LAMP::Geometry::ExactPointKey,LAMP::Geometry::ExactPointKeyHash> additionalKeys;
    const auto requestAdditional = [&](Vector2D point) {
        const auto key=LAMP::Geometry::exactPointKey(point);
        if (increments.find(key)==increments.end() && additionalKeys.insert(key).second)
            additionalPoints.push_back(point);
    };
    for (const auto& delta:increment.edgeMidpoints)
    {
        if (!std::isfinite(delta.value)) continue;
        const auto found=retainedEdges.find(delta.vertices);
        require(found!=retainedEdges.end(), "retained quartic patch lost an independent edge");
        const auto& edge=result.patch.edgeMidpoints[found->second];
        const auto a=mesh.verts(mesh.vertIndexToLabel[edge.vertices[0]]).coord();
        const auto b=mesh.verts(mesh.vertIndexToLabel[edge.vertices[1]]).coord();
        requestAdditional((a+b)*0.5);
    }
    for (int c=0;c<mesh.cellIndexNum();++c)
        if (previous.activeCells[c]) requestAdditional(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord());
    const auto additionalValues=sampleBatch(additionalPoints);
    for (std::size_t i=0;i<additionalPoints.size();++i)
    {
        increments.emplace(LAMP::Geometry::exactPointKey(additionalPoints[i]),additionalValues[i]);
        result.maxTraceRatio=std::max(result.maxTraceRatio,additionalValues[i].traceRatio);
    }
    result.scalarIncrementEvaluations+=additionalPoints.size();
    for (std::size_t edgeIndex=0;edgeIndex<increment.edgeMidpoints.size();++edgeIndex)
    {
        const auto& delta=increment.edgeMidpoints[edgeIndex];
        if (!std::isfinite(delta.value)) continue; // Generator's topology-owned dependent side.
        const auto found = retainedEdges.find(delta.vertices);
        require(found != retainedEdges.end(), "retained quartic patch lost an independent edge");
        auto &edge = result.patch.edgeMidpoints[found->second];
        const double oldValue = edge.value;
        // Values, like regular nodal values, follow the characteristic at
        // their OWN point. Q4 supplies derivatives, not a replacement value:
        // its midpoint stencil may cross a stationary wall endpoint or the
        // zero-speed cutoff even when this point has no displacement.
        const auto a = mesh.verts(mesh.vertIndexToLabel[edge.vertices[0]]).coord();
        const auto b = mesh.verts(mesh.vertIndexToLabel[edge.vertices[1]]).coord();
        add(edge.value, scalarIncrement((a + b) * 0.5));
        if (pullbackJet) {
            const auto& traced=increments.at(LAMP::Geometry::exactPointKey((a+b)*0.5));
            const bool horizontal=a.y==b.y;
            const double x=mapX.edgeMidpoints[edgeIndex].normalDerivative;
            const double y=mapY.edgeMidpoints[edgeIndex].normalDerivative;
            if (traced.inflow || (hasInflowSamples && affected.edgeMidpoints[edgeIndex]))
                edge.normalDerivative=boundaryData().edgeMidpoints[edgeIndex].normalDerivative;
            else if (traced.traceRatio!=0 || x!=0 || y!=0) {
                const auto jet=source->scalarJetNative(traced.departure);
                edge.normalDerivative=jet.dx*(x+(horizontal?0:1))+jet.dy*(y+(horizontal?1:0));
                require(std::isfinite(edge.normalDerivative),"nonfinite pulled-back edge jet");
            }
            if (options.radialAxisSymmetry && a.x==0 && b.x==0) edge.normalDerivative=0;
        } else add(edge.normalDerivative, delta.normalDerivative);
        result.changedEdgeValues += edge.value != oldValue;
        result.maxValueIncrement = std::max(result.maxValueIncrement, std::abs(edge.value - oldValue));
    }
    for (int c = 0; c < mesh.cellIndexNum(); ++c)
    {
        if (!previous.activeCells[c]) continue;
        add(result.patch.cellCenters[c], scalarIncrement(mesh.cells(mesh.cellIndexToLabel[c]).centerCoord()));
        result.changedCenterValues += result.patch.cellCenters[c] != previous.cellCenters[c];
        result.maxValueIncrement = std::max(result.maxValueIncrement,
            std::abs(result.patch.cellCenters[c] - previous.cellCenters[c]));
    }
    result.timings.dataUpdate=elapsed();
    finishCandidate(mesh, options, masks, characteristic, result);
    result.timings.candidateGeometry=elapsed();
    return result;
}
QuarticCandidate advanceQuarticNative(const QuadtreeMesh& mesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch& previous,double dt,
    const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity,const Options& options)
{ return advanceQuarticImpl(mesh,previous,dt,newVelocity,oldVelocity,options,false); }
QuarticCandidate advanceQuarticNativeJet(const QuadtreeMesh& mesh,
    const LAMP::Geometry::ActiveQuarticHermitePatch& previous,double dt,
    const LevelSetBoundaryTrace::Velocity& newVelocity,
    const LevelSetBoundaryTrace::Velocity& oldVelocity,const Options& options)
{ return advanceQuarticImpl(mesh,previous,dt,newVelocity,oldVelocity,options,true); }
} // namespace ActiveLevelSetTransport
