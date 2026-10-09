#include "Core/Geometry/LevelSetInterfacePanels.h"

#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/LevelSetInterfaceQuadrature.h"
#include "Core/Geometry/Detail/LevelSetLocalProjectionSearch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

namespace
{
namespace E = Quantum::Geometry::Detail::ExactExpansion;
namespace Search = Quantum::Geometry::Detail::LocalInterfaceProjection;
using Search::require;
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Frame = LevelSetHermiteGeometry::Frame;
using Chain = LevelSetInterfacePanels::NativeGraphChain;
using Query = LevelSetInterfacePanels::NativeProjectionQuery;
using Result = LevelSetInterfacePanels::LocalProjection;
using Chart = Quantum::Geometry::CellChartPoint;
using Coordinates = std::array<E::Expansion, 2>;

E::Expansion physicalComponent(const Snapshot& source,int cell,const Chart& chart,int axis)
{
    const auto bounds=source.nativeCellBounds(cell);
    const double low=axis==0?bounds[0].x:bounds[0].y;
    const double high=axis==0?bounds[1].x:bounds[1].y;
    const auto local=axis==0?chart.x:chart.y;
    const auto width=E::subtract(E::fromDouble(high),E::fromDouble(low));
    const auto tail=E::multiply(width,E::fromDouble(local.offset));
    return local.fromUpper?E::subtract(E::fromDouble(high),tail):E::add(E::fromDouble(low),tail);
}
struct PreparedGraphCoordinates
{
    std::array<Coordinates,10> fixed;
    Coordinates lower,upper;
};
using PreparedCoordinates=std::vector<PreparedGraphCoordinates>;

struct Sample
{
    Result result;
    double derivative = 0.0;
};

// Physical distances on ACTUAL native rectangles. Low/high chart offsets
// remain separate expansions; no subtraction of rounded display coordinates.
class NativeProjectionContext
{
  public:
    NativeProjectionContext(const Snapshot &source, Vector2D spacing, Vector2D query) : source_(source)
    {
        E::validateEnvironment();
        require(source.usesNativeLeaves() && std::isfinite(source.value(query)),
                "native load query requires an actual active-domain point");
        const auto lo = source.minCoord(), hi = source.maxCoord();
        const long double scale =
            std::max({static_cast<long double>(spacing.x), static_cast<long double>(spacing.y),
                      std::abs(static_cast<long double>(lo.x)), std::abs(static_cast<long double>(lo.y)),
                      std::abs(static_cast<long double>(hi.x)), std::abs(static_cast<long double>(hi.y))});
        const long double bound = 2048.0L * std::numeric_limits<double>::epsilon() * scale;
        require(std::isfinite(bound) && bound > 0 &&
                    bound <= std::sqrt(static_cast<long double>(std::numeric_limits<double>::epsilon())) *
                                 std::min(spacing.x, spacing.y),
                "coordinate precision cannot resolve the requested boundary projection");
        tolerance = static_cast<double>(bound);
        target_ = {{E::fromDouble(query.x), E::fromDouble(query.y)}};
    }
    NativeProjectionContext(const Snapshot &source, Vector2D spacing, const Query &query) : source_(source)
    {
        E::validateEnvironment();
        require(source.usesNativeLeaves() && source.nativeCellActive(query.cellIndex),
                "native projection requires an active unchanged leaf anchor");
        Quantum::Geometry::validateCellChart(query.cellChartPoint);
        const double h = std::min(spacing.x, spacing.y);
        radius = h / 4.0;
        tolerance = 1e-9 * h;
        const auto lo = source.minCoord(), hi = source.maxCoord();
        boxTolerance =
            512.0 * std::numeric_limits<double>::epsilon() *
            std::max({spacing.x, spacing.y, std::abs(lo.x), std::abs(lo.y), std::abs(hi.x), std::abs(hi.y)});
        require(std::isfinite(query.displacement.x) && std::isfinite(query.displacement.y) &&
                    std::isfinite(radius) && radius > 0 && std::isfinite(tolerance) && tolerance > 0 &&
                    std::isfinite(boxTolerance) &&
                    std::hypot(query.displacement.x, query.displacement.y) <= h / 2.0,
                "local material displacement exceeds half a source cell");
        for (int axis = 0; axis < 2; ++axis)
            target_[axis] = E::add(component(query.cellIndex, query.cellChartPoint, axis),
                                   E::fromDouble(axis == 0 ? query.displacement.x : query.displacement.y));
    }

    Sample evaluate(int index, int cell, Frame frame, double fraction, bool endpoint = false,
                    bool begin = false, LevelSetMappedInterface::Endpoint root = {},
                    const Coordinates* coordinates=nullptr) const
    {
        require(frame.hasCellChartPoint && !frame.hasCellLocalPoint && frame.cellI == -1 &&
                    frame.cellJ == -1 && frame.nativeCellIndex == cell,
                "native projection lost its source-jet provenance");
        const auto dx = coordinates?E::subtract((*coordinates)[0],target_[0]):difference(frame,0);
        const auto dy = coordinates?E::subtract((*coordinates)[1],target_[1]):difference(frame,1);
        Sample sample;
        sample.derivative = E::estimate(E::add(E::multiply(dx, E::fromDouble(-frame.normal.y)),
                                               E::multiply(dy, E::fromDouble(frame.normal.x))));
        auto &result = sample.result;
        result.frame = std::move(frame);
        result.frame.stationarityResidual = std::abs(sample.derivative);
        result.panelIndex = index;
        result.endpoint = endpoint;
        result.atPanelBegin = begin;
        result.distance = std::hypot(E::estimate(dx), E::estimate(dy));
        result.inwardDerivative = endpoint ? (begin ? sample.derivative : -sample.derivative) : 0.0;
        result.kktResidual = endpoint ? std::max(0.0, -result.inwardDerivative) : std::abs(sample.derivative);
        result.nativeGraphFraction = fraction;
        if (endpoint)
            result.root = {root.side, root.rootOrdinal};
        require(std::isfinite(result.distance) && std::isfinite(result.frame.levelSetResidual) &&
                    std::isfinite(result.frame.stationarityResidual) && std::isfinite(result.kktResidual),
                "non-finite local projection residual");
        return sample;
    }

    // The load-search upper bound needs only this same exact-coordinate
    // distance, not a stationarity jet for every endpoint of every graph.
    // This never accepts a projection: retained candidates still pass the
    // complete frame, scalar, orientation and KKT checks in evaluate/search.
    double endpointDistance(const Coordinates& coordinates) const
    {
        const auto dx=E::subtract(coordinates[0],target_[0]);
        const auto dy=E::subtract(coordinates[1],target_[1]);
        const double distance=std::hypot(E::estimate(dx),E::estimate(dy));
        require(std::isfinite(distance),"non-finite native endpoint distance bound");
        return distance;
    }

    double boxDistance(int cell,const Frame& first,const Frame& last,bool coordinateIsX,
        Vector2D localLower,Vector2D localUpper,const PreparedGraphCoordinates* prepared=nullptr) const
    {
        const auto lowerChart=Quantum::Geometry::cellChartFromLocal(localLower);
        const auto upperChart=Quantum::Geometry::cellChartFromLocal(localUpper);
        double distance[2]{};
        for (int axis = 0; axis < 2; ++axis)
        {
            auto lower=prepared?prepared->lower[axis]:component(cell,lowerChart,axis);
            auto upper=prepared?prepared->upper[axis]:component(cell,upperChart,axis);
            if (!prepared && axis==(coordinateIsX?0:1)) {
                // The retained graph is monotone in its declared coordinate.
                // Its exact endpoint-relative interval bounds that coordinate;
                // The other coordinate retains the discovered chart's full
                // dyadic subbox, which already bounds this original graph.
                // Whole-cell graphs retain their original bounds unchanged.
                // No sampled/display envelope and no adjacency relaxation.
                lower=component(cell,first.cellChartPoint,axis);
                upper=component(cell,last.cellChartPoint,axis);
                if (E::sign(E::subtract(lower,upper))>0) std::swap(lower,upper);
            }
            const auto a = E::subtract(lower, target_[axis]);
            const auto b = E::subtract(upper, target_[axis]);
            if (E::sign(a) > 0)
                distance[axis] = E::estimate(a);
            else if (E::sign(b) < 0)
                distance[axis] = -E::estimate(b);
        }
        return std::hypot(distance[0], distance[1]);
    }

    double targetFraction(const Frame &first, const Frame &last, bool coordinateIsX,
        const PreparedGraphCoordinates* prepared=nullptr) const
    {
        const int axis = coordinateIsX ? 0 : 1;
        const auto a = prepared?prepared->fixed[0][axis]:component(first.nativeCellIndex,first.cellChartPoint,axis);
        const auto b = prepared?prepared->fixed[1][axis]:component(last.nativeCellIndex,last.cellChartPoint,axis);
        // Only a search seed: endpoints/topological identity are never inferred
        // from this terminal quotient. The graph sampler validates the trial.
        return E::estimate(E::subtract(target_[axis], a)) / E::estimate(E::subtract(b, a));
    }

    int compareDistance(const Frame &a, const Frame &b) const
    {
        const auto squared = [&](const Frame &frame) {
            const auto dx = difference(frame, 0), dy = difference(frame, 1);
            return E::add(E::multiply(dx, dx), E::multiply(dy, dy));
        };
        return E::sign(E::subtract(squared(a), squared(b)));
    }
    bool samePosition(const Frame &a, const Frame &b) const
    {
        for (int axis = 0; axis < 2; ++axis)
            if (E::sign(E::subtract(component(a.nativeCellIndex, a.cellChartPoint, axis),
                                    component(b.nativeCellIndex, b.cellChartPoint, axis))) != 0)
                return false;
        return true;
    }

    void accept(const Result &result) const
    {
        if(result.distance>radius) {
            std::ostringstream message;message<<std::setprecision(17)
                <<"local material projection exceeds a quarter source cell: distance="<<result.distance
                <<", radius="<<radius<<", graph="<<result.panelIndex<<", fraction="<<result.nativeGraphFraction
                <<", point=("<<result.frame.point.x<<','<<result.frame.point.y<<')';
            require(false,message.str().c_str());
        }
        require(result.frame.levelSetResidual <= tolerance && result.kktResidual <= tolerance,
                "local material projection exceeds the original level-set/stationarity or endpoint KKT "
                "tolerance");
    }

    double radius = 0, tolerance = 0, boxTolerance = 0;

  private:
    const Snapshot &source_;
    std::array<E::Expansion, 2> target_;
    E::Expansion component(int cell, const Chart &chart, int axis) const
    {
        return physicalComponent(source_,cell,chart,axis);
    }
    E::Expansion difference(const Frame &frame, int axis) const
    {
        return E::subtract(component(frame.nativeCellIndex, frame.cellChartPoint, axis), target_[axis]);
    }
};

class NativeProjectionLayout
{
  public:
    using Sample = ::Sample;
    NativeProjectionLayout(const Chain &chain, const std::vector<std::array<int, 2>> &roots,
        const PreparedCoordinates* coordinates=nullptr)
        : chain_(chain), roots_(roots),coordinates_(coordinates)
    {
        require(!chain.graphs().empty() && roots.size() == chain.graphs().size(),
                "native projection requires an owned complete root chain");
    }
    int count() const
    {
        return static_cast<int>(chain_.graphs().size());
    }
    bool closed() const
    {
        return chain_.closed();
    }
    int identity(int index, bool begin) const
    {
        const int id = roots_.at(index)[begin ? 0 : 1];
        require(id >= 0, "native projection lost its discovered endpoint identity");
        return id;
    }
    Sample endpoint(const NativeProjectionContext &context, int index, bool begin) const
    {
        identity(index, begin);
        const auto &graph = chain_.graphs().at(index);
        return context.evaluate(index, graph.nativeCellIndex, frame(index, begin ? 0.0 : 1.0),
                                begin ? 0.0 : 1.0, true, begin, begin ? graph.begin : graph.end,
                                coordinates_?&coordinates_->at(index).fixed[begin?0:1]:nullptr);
    }
    double endpointDistance(const NativeProjectionContext& context,int index,bool begin) const
    {
        identity(index,begin);
        return coordinates_?context.endpointDistance(coordinates_->at(index).fixed[begin?0:1])
            :endpoint(context,index,begin).result.distance;
    }
    double boxDistance(const NativeProjectionContext &context, int index) const
    {
        const auto& graph=chain_.graphs().at(index);
        return context.boxDistance(graph.nativeCellIndex,frame(index,0),frame(index,1),graph.coordinateIsX,
            graph.localLower,graph.localUpper,coordinates_?&coordinates_->at(index):nullptr);
    }

    std::vector<Sample> samples(const NativeProjectionContext &context, int index) const
    {
        std::vector<Sample> samples;
        samples.push_back(endpoint(context, index, true));
        for (double x : Quantum::Geometry::Detail::InterfaceQuadrature::gaussX)
            samples.push_back(evaluate(context, index, (1.0 + x) / 2.0));
        samples.push_back(endpoint(context, index, false));
        const double seed = context.targetFraction(samples.front().result.frame, samples.back().result.frame,
            chain_.graphs()[index].coordinateIsX,coordinates_?&coordinates_->at(index):nullptr);
        if (std::isfinite(seed) && seed > 0 && seed < 1)
            samples.push_back(evaluate(context, index, seed));
        std::sort(samples.begin(), samples.end(), [](const auto &a, const auto &b) {
            return a.result.nativeGraphFraction < b.result.nativeGraphFraction;
        });
        return samples;
    }
    Sample minimum(const NativeProjectionContext &context, int index, Sample first, Sample last) const
    {
        return loadMinimum(context, index, std::move(first), std::move(last), context.tolerance);
    }
    Sample loadMinimum(const NativeProjectionContext &context, int index, Sample first, Sample last,
                       double tolerance) const
    {
        return Search::stationaryMinimum(
            std::move(first), std::move(last), tolerance,
            [](const Sample &sample) { return sample.result.nativeGraphFraction; },
            [&](double fraction) { return evaluate(context, index, fraction); });
    }

  private:
    const Chain &chain_;
    const std::vector<std::array<int, 2>> &roots_;
    const PreparedCoordinates* coordinates_;
    // Query-scoped exact-key reuse on one immutable source; no cross-time cache.
    mutable std::map<std::pair<int, double>, Frame> frameCache_;
    Frame frame(int index, double fraction) const
    {
        // Prepared material chains already own these immutable fixed frames.
        // Avoid copying all endpoints/GL8 nodes into a per-query tree map;
        // only query-dependent stationary iterates need that local cache.
        if(coordinates_) {
            if(fraction==0 || fraction==1)return chain_.frameAt(index,fraction);
            for(double x:Quantum::Geometry::Detail::InterfaceQuadrature::gaussX)
                if(fraction==(1.0+x)/2.0)return chain_.frameAt(index,fraction);
        }
        const auto key = std::make_pair(index, fraction);
        auto found = frameCache_.find(key);
        if (found == frameCache_.end())
            found = frameCache_.emplace(key, chain_.frameAt(index, fraction)).first;
        return found->second;
    }
    Sample evaluate(const NativeProjectionContext &context, int index, double fraction) const
    {
        const Coordinates* coordinates=nullptr;
        if(coordinates_) {
            const auto& nodes=Quantum::Geometry::Detail::InterfaceQuadrature::gaussX;
            for(std::size_t k=0;k<nodes.size();++k)
                if(fraction==(1.0+nodes[k])/2.0) {coordinates=&coordinates_->at(index).fixed[k+2];break;}
        }
        return context.evaluate(index, chain_.graphs()[index].nativeCellIndex, frame(index, fraction),
                                fraction,false,false,{},coordinates);
    }
};
} // namespace

struct LevelSetInterfacePanels::NativeGraphChain::ProjectionCoordinates
{
    PreparedCoordinates graphs;
};

void LevelSetInterfacePanels::NativeGraphChain::prepareProjectionCoordinates()
{
    E::validateEnvironment();
    require(source_ && fixedProjectionFrames_.size()==graphs_.size(),
        "native projection coordinates require complete source-bound fixed frames");
    auto result=std::make_shared<ProjectionCoordinates>();
    result->graphs.resize(graphs_.size());
    for(std::size_t g=0;g<graphs_.size();++g) {
        const auto& graph=graphs_[g];auto& prepared=result->graphs[g];
        for(std::size_t k=0;k<prepared.fixed.size();++k) {
            const auto& frame=fixedProjectionFrames_[g][k];
            require(frame.hasCellChartPoint && frame.nativeCellIndex==graph.nativeCellIndex,
                "prepared projection coordinate lost its native source");
            for(int axis=0;axis<2;++axis)
                prepared.fixed[k][axis]=physicalComponent(*source_,graph.nativeCellIndex,frame.cellChartPoint,axis);
        }
        const auto lower=Quantum::Geometry::cellChartFromLocal(graph.localLower);
        const auto upper=Quantum::Geometry::cellChartFromLocal(graph.localUpper);
        for(int axis=0;axis<2;++axis) {
            prepared.lower[axis]=physicalComponent(*source_,graph.nativeCellIndex,lower,axis);
            prepared.upper[axis]=physicalComponent(*source_,graph.nativeCellIndex,upper,axis);
            if(axis==(graph.coordinateIsX?0:1)) {
                prepared.lower[axis]=prepared.fixed[0][axis];prepared.upper[axis]=prepared.fixed[1][axis];
                if(E::sign(E::subtract(prepared.lower[axis],prepared.upper[axis]))>0)
                    std::swap(prepared.lower[axis],prepared.upper[axis]);
            }
        }
    }
    projectionCoordinates_=std::move(result);
}

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::NativeGraphChain::projectLocal(
    const NativeProjectionQuery &query) const
{
    require(static_cast<bool>(source_), "native material projection requires a retained source");
    const NativeProjectionContext context(*source_, minimumSpacing_, query);
    return Search::search(context, NativeProjectionLayout(*this,rootNodes_,
        projectionCoordinates_?&projectionCoordinates_->graphs:nullptr));
}

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::NativeGraphChain::projectLocalEndpoint(
    const NativeProjectionQuery &query, int graphIndex, bool atGraphBegin) const
{
    require(source_ && graphIndex >= 0 && static_cast<std::size_t>(graphIndex) < graphs_.size(),
            "native endpoint projection requires a retained source and valid graph index");
    const NativeProjectionContext context(*source_, minimumSpacing_, query);
    const NativeProjectionLayout layout(*this,rootNodes_,projectionCoordinates_?&projectionCoordinates_->graphs:nullptr);
    auto result = layout.endpoint(context, graphIndex, atGraphBegin).result;
    context.accept(result);
    Search::guardBranches(context, layout, result);
    return result;
}

LevelSetInterfacePanels::LocalProjection LevelSetInterfacePanels::NativeGraphChain::projectedFrameUncached(
    Vector2D query) const
{
    require(source_ && !graphs_.empty(), "native load projection requires an owned retained curve");
    const NativeProjectionContext context(*source_, minimumSpacing_, query);
    LocalProjection best;
    bool found = false;
    const auto consider = [&](LocalProjection candidate) {
        const int order = found ? context.compareDistance(candidate.frame, best.frame) : -1;
        if (!found || order < 0 ||
            (order == 0 && candidate.endpoint && context.samePosition(candidate.frame, best.frame)))
        {
            best = std::move(candidate);
            found = true;
        }
    };
    // Native leaves use the original retained-graph search directly. There
    // is no uniform Newton seed lattice, rounded-point inverse or solid jet.
    for (auto candidate : loadProjectionCandidates(query))
        consider(std::move(candidate));
    require(found,
            "native load projection has no retained-curve stationary or admissible endpoint candidate");
    return best;
}
std::vector<LevelSetInterfacePanels::LocalProjection>
LevelSetInterfacePanels::NativeGraphChain::loadProjectionCandidatesUncached(Vector2D query) const
{
    require(source_ && !graphs_.empty(), "native load projection requires an owned retained curve");
    const NativeProjectionContext context(*source_, minimumSpacing_, query);
    const NativeProjectionLayout layout(*this,rootNodes_,projectionCoordinates_?&projectionCoordinates_->graphs:nullptr);
    auto candidates=Search::loadCandidates(context, layout, context.tolerance);
    if (!closed_)
        for (int end = 0; end < 2; ++end)
        {
            const auto kind = end == 0 ? beginKind_ : endKind_;
            require(kind == EndpointKind::Axis || kind == EndpointKind::DeclaredBoundary || kind==EndpointKind::SymmetryPlane,
                    "native load projection encountered an undeclared endpoint");
            auto candidate = layout.endpoint(context, end == 0 ? 0 : layout.count() - 1, end == 0).result;
            require(candidate.frame.levelSetResidual <= context.tolerance,
                    "original native load endpoint has an unresolved scalar residual");
            if (kind == EndpointKind::Axis)
                require(axisymmetric_ && candidate.frame.point.x == 0 && candidate.frame.normal.x == 0,
                        "original native load pole lost radial symmetry");
            if (kind==EndpointKind::SymmetryPlane) {
                const auto* plane=endpointSymmetryPlane(end!=0);
                require(plane && (plane->normalAxis==0?candidate.frame.point.x:candidate.frame.point.y)==plane->coordinate &&
                    (plane->normalAxis==0?candidate.frame.normal.x:candidate.frame.normal.y)==0,
                    "original native load endpoint lost section symmetry");
            }
            if (candidate.inwardDerivative >= -context.tolerance)
                candidates.push_back(std::move(candidate));
        }
    return candidates;
}
