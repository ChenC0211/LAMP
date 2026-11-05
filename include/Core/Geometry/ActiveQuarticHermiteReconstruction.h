#pragma once

#include "Core/Base/Vector2D.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include <cstdint>
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

class QuadtreeMesh;

namespace ActiveQuarticHermiteReconstruction
{
struct Options
{
    bool radialAxisSymmetry=false;
    // Optional native-vertex masks for smooth one-sided sampling charts.
    // A datum uses the first chart containing all its defining vertices
    // (vertex, edge endpoints or cell corners), then the second, else the
    // entire active domain. A seam therefore takes the first chart's jet.
    // The ENTIRE chosen support box must lie in that chart; no extrapolation.
    std::vector<std::uint8_t> oneSidedSamples,extensionSideSamples;
};

struct Report
{
    std::size_t stencilRequests = 0, stencilBuilds = 0;
    std::size_t scalarRequests = 0, scalarEvaluations = 0;
    bool reusedTopology = false;
};

class StencilCache;
// Immutable, owning geometry-only batch. Points are unique by their exact
// binary64 coordinates and ordered by first use in the original node/edge/
// center traversal. No mesh borrow, scalar, velocity or time layer is stored.
// An existing batch remains valid when its cache is invalidated or destroyed.
class Prepared
{
  public:
    struct Impl;
    const std::vector<Vector2D>& samplePoints() const;
    struct AffectedData
    {
        std::vector<std::uint8_t> nodes,edgeMidpoints;
    };
    // Independent derivative data whose ORIGINAL support includes a marked
    // sample. Edge order is exactly the order returned by reconstruct().
    AffectedData affectedBy(const std::vector<std::uint8_t>& markedSamples) const;
  private:
    explicit Prepared(std::shared_ptr<const Impl> impl);
    std::shared_ptr<const Impl> impl_;
    friend Prepared prepare(const QuadtreeMesh&, const std::vector<std::uint8_t>&,
        const Options&, StencilCache*, Report*);
    friend LAMP::Geometry::ActiveQuarticHermitePatch reconstruct(
        const Prepared&, const std::vector<double>&);
};

Prepared prepare(const QuadtreeMesh& mesh, const std::vector<std::uint8_t>& activeCells,
    const Options& options={}, StencilCache* cache=nullptr, Report* report=nullptr);
// Current values correspond exactly to samplePoints(). Independent datum fits
// may run concurrently; all use the SAME Q4 arithmetic. Exceptions join before
// return; no partial patch is published. The supplied array is never retained.
LAMP::Geometry::ActiveQuarticHermitePatch reconstruct(
    const Prepared& prepared, const std::vector<double>& values);

// Optional, nonphysical workspace for repeated calls on the SAME native
// topology/activity/charts. Owns exact identity data and a prepared batch,
// never scalar values, callbacks or borrowed mesh arrays. A changed identity
// clears the selections before use. Sequential use only; not a global cache.
class StencilCache
{
  public:
    struct Impl;
    StencilCache();
    ~StencilCache();
    StencilCache(const StencilCache&) = delete;
    StencilCache& operator=(const StencilCache&) = delete;
  private:
    std::unique_ptr<Impl> impl_;
    friend Prepared prepare(const QuadtreeMesh&,
        const std::vector<std::uint8_t>&,
        const Options&, StencilCache*, Report*);
};

// Generate explicit independent quartic DATA from a supplied physical scalar
// using the existing Q4 recursive extension.
// Convenience port: prepare, sequentially sample on the calling thread, then
// evaluate the prepared batch. Successful calls keep first-use callback order.
// The scalar must be repeatable
// throughout this call, with well-defined fluid-side boundary values. Queries
// are cached only for this invocation; no callback/mesh reference is retained.
// Optional cache retains geometry-only support/sample layout; current topology, native
// coordinates, active cells and both chart masks are checked on every call.
// Native nodal values are sampled at the exact stored vertex, never replaced
// by fitted values. Corner spacing is the finest incident active spacing;
// full-edge data use the actual side length and finest incident normal width;
// center data use their own leaf dimensions. No replacement mesh is created.
// Dependent hanging jets/fine-side data are left NaN for the shared scalar
// constraint owner. Returns no corner sector/contact law or moving candidate;
// the caller must construct/validate geometry before committing any state.
LAMP::Geometry::ActiveQuarticHermitePatch reconstruct(const QuadtreeMesh& mesh,
    const std::vector<std::uint8_t>& activeCells,const std::function<double(Vector2D)>& scalar,
    const Options& options={}, StencilCache* cache=nullptr, Report* report=nullptr);
}
