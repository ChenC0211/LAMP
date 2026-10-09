#pragma once

#include "Core/Geometry/Detail/ExactExpansion.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Quantum::Geometry::Detail
{
namespace
{
namespace LocalInterfaceProjection
{
inline void require(bool valid, const char *reason)
{
    if (!valid)
        throw std::runtime_error(std::string("[LevelSetInterfacePanels] ") + reason);
}
inline int projectionNeighbor(int index, bool begin, int count, bool closed)
{
    const int next = index + (begin ? -1 : 1);
    return next < 0 ? (closed ? count - 1 : -1) : next == count ? (closed ? 0 : -1) : next;
}

// Safeguarded Brent interpolation in the layout's retained curve parameter.
// Uniform layout supplies its original local chosen coordinate; native layout
// supplies the original graph fraction, retaining endpoint-relative charts.
// Only the trial coordinate changes. Keep the original sign bracket, 128-query
// budget, representably-interior termination and final stationarity gate.
// An interpolation is NOT a root/nearest-point certificate; the actual retained
// graph is evaluated at every trial. No residual-only early exit or x tolerance.
template <class Sample, class Coordinate, class Evaluate>
Sample stationaryMinimum(Sample first, Sample last, double stationaryTolerance,
                         const Coordinate &coordinateOf, const Evaluate &evaluate)
{
    require(first.derivative < 0.0 && last.derivative > 0.0, "invalid local minimum bracket");
    Sample best = std::abs(first.derivative) < std::abs(last.derivative) ? first : last;
    bool interior = false;
    double previousCoordinate = coordinateOf(first), previousDerivative = first.derivative;
    double lastStep = coordinateOf(last) - coordinateOf(first), precedingStep = lastStep;
    for (int iteration = 0; iteration < 128; ++iteration)
    {
        const double a = coordinateOf(first), b = coordinateOf(last);
        const auto midpointExpansion = ExactExpansion::scalePow2(
            ExactExpansion::add(ExactExpansion::fromDouble(a), ExactExpansion::fromDouble(b)), -1);
        const double midpoint = ExactExpansion::estimate(midpointExpansion);
        const auto strictlyInterior = [&](double x) { return x > std::min(a, b) && x < std::max(a, b); };
        if (!strictlyInterior(midpoint))
            break;
        const bool firstIsCloser = std::abs(first.derivative) < std::abs(last.derivative);
        const double x = firstIsCloser ? a : b, opposite = firstIsCloser ? b : a;
        const double f = firstIsCloser ? first.derivative : last.derivative;
        const double oppositeF = firstIsCloser ? last.derivative : first.derivative;
        // If the old best survived the last trial, that trial is the new
        // opposite end; use the two distinct bracket samples for a secant.
        if (previousCoordinate == x)
        {
            previousCoordinate = opposite;
            previousDerivative = oppositeF;
        }
        double coordinate = midpoint;
        bool interpolated = false;
        if (std::abs(previousDerivative) > std::abs(f) && std::isfinite(precedingStep))
        {
            // Secant for two distinct samples; inverse quadratic interpolation
            // when the previous best is a third sample. Ratios avoid products
            // of small derivatives. Degenerate/non-finite proposals bisect.
            const double ratio = f / previousDerivative;
            double step;
            if (previousCoordinate == opposite)
                step = (opposite - x) * ratio / (ratio - 1.0);
            else
            {
                const double q = previousDerivative / oppositeF, r = f / oppositeF;
                step = -ratio * ((opposite - x) * q * (q - r) -
                                (x - previousCoordinate) * (r - 1.0)) /
                       ((q - 1.0) * (r - 1.0) * (ratio - 1.0));
            }
            if (std::isfinite(step) && step * (opposite > x ? 1.0 : -1.0) > 0.0 &&
                std::abs(step) < .75 * std::abs(opposite - x) &&
                std::abs(step) < .5 * std::abs(precedingStep))
            {
                double proposed = x + step;
                // Do not loosen convergence when interpolation rounds back to
                // an endpoint. Test its adjacent representable interior point.
                if (proposed == x) proposed = std::nextafter(x, opposite);
                if (strictlyInterior(proposed))
                {
                    coordinate = proposed;
                    interpolated = true;
                }
            }
        }
        precedingStep = interpolated ? lastStep : midpoint - x;
        lastStep = coordinate - x;
        previousCoordinate = x;
        previousDerivative = f;
        auto trial = evaluate(coordinate);
        trial.result.frame.iterations = iteration + 1;
        if (!interior || std::abs(trial.derivative) < std::abs(best.derivative))
            best = trial;
        interior = true;
        if (trial.derivative == 0.0)
            return trial;
        if (trial.derivative < 0.0)
            first = std::move(trial);
        else
            last = std::move(trial);
    }
    require(interior && !best.result.endpoint && best.result.kktResidual <= stationaryTolerance,
            "local stationary bracket has no resolved interior projection");
    return best;
}

// Nonadjacent-branch guard using the layout's enclosing geometric bounds,
// never a display-curve envelope. Native graphs tighten only their certified
// monotone coordinate to the retained endpoint interval; the other coordinate
// keeps its discovered chart bound (the whole cell for unsplit graphs).
// Thresholds and topological exclusions are fixed.
template <class Context, class Layout, class Result>
void guardBranches(const Context &context, const Layout &layout, const Result &result)
{
    const int count = layout.count();
    for (int index = 0; index < count; ++index)
    {
        if (index == result.panelIndex ||
            index == projectionNeighbor(result.panelIndex, true, count, layout.closed()) ||
            index == projectionNeighbor(result.panelIndex, false, count, layout.closed()))
            continue;
        require(layout.boxDistance(context, index) > result.distance + context.boxTolerance,
                "material projection is ambiguous near a nonadjacent interface branch");
    }
}

// Original retained load-curve search. Unlike material remapping this is not
// bounded by h/4: actual endpoint distances bound the candidate-cell search.
// Open endpoint KKT selection remains with the caller. Both layouts retain
// their own chosen coordinate arithmetic, samples and original root identity.
template <class Context, class Layout>
auto loadCandidates(const Context &context, const Layout &layout, double tolerance)
{
    using Sample = typename Layout::Sample;
    using Result = decltype(layout.endpoint(context, 0, true).result);
    std::vector<Result> candidates;
    double upperDistance = std::numeric_limits<double>::infinity();
    for (int index = 0; index < layout.count(); ++index)
        for (bool begin : {true, false})
            upperDistance = std::min(upperDistance, layout.endpointDistance(context, index, begin));
    require(std::isfinite(upperDistance), "retained load search has no actual endpoint distance bound");
    const auto add = [&](Sample sample) {
        require(sample.result.frame.levelSetResidual <= tolerance && sample.result.kktResidual <= tolerance,
                "retained load graph exceeds the original scalar/stationarity tolerance");
        candidates.push_back(std::move(sample.result));
    };
    for (int index = 0; index < layout.count(); ++index)
    {
        if (layout.boxDistance(context, index) > upperDistance + tolerance)
            continue;
        const auto samples = layout.samples(context, index);
        for (const auto &sample : samples)
        {
            if (!sample.result.endpoint)
            {
                if (sample.derivative == 0.0)
                    add(sample);
                continue;
            }
            const bool begin = sample.result.atPanelBegin;
            const int neighbor = projectionNeighbor(index, begin, layout.count(), layout.closed());
            if (neighbor < 0)
                continue;
            require(layout.identity(index, begin) == layout.identity(neighbor, !begin),
                    "retained load graph endpoints have different original root identities");
            const auto adjacent = layout.endpoint(context, neighbor, !begin);
            if (std::abs(sample.derivative) > tolerance || std::abs(adjacent.derivative) > tolerance)
                continue;
            auto crossing = sample;
            crossing.result.endpoint = false;
            crossing.result.inwardDerivative = 0.0;
            crossing.result.kktResidual = std::abs(sample.derivative);
            add(std::move(crossing));
        }
        for (std::size_t k = 1; k < samples.size(); ++k)
            if (samples[k - 1].derivative < 0.0 && samples[k].derivative > 0.0)
                add(layout.loadMinimum(context, index, samples[k - 1], samples[k], tolerance));
    }
    return candidates;
}

// Same material-search policy for both layouts. Geometry adapters provide only
// retained frames, root identities, sample coordinates and physical metrics.
template <class Context, class Layout> auto search(const Context &context, const Layout &layout)
{
    using Sample = typename Layout::Sample;
    std::vector<Sample> candidates;
    const int count = layout.count();
    const auto addEndpoint = [&](Sample sample) {
        if (sample.result.inwardDerivative < 0.0)
            return;
        const int index = sample.result.panelIndex;
        const bool begin = sample.result.atPanelBegin;
        const auto identity = layout.identity(index, begin);
        const int adjacent = projectionNeighbor(index, begin, count, layout.closed());
        if (adjacent >= 0)
        {
            require(identity == layout.identity(adjacent, !begin),
                    "ordered local panel endpoints have different original root identities");
            if (layout.endpoint(context, adjacent, !begin).result.inwardDerivative < 0.0)
                return;
        }
        for (const auto &prior : candidates)
            if (prior.result.endpoint &&
                identity == layout.identity(prior.result.panelIndex, prior.result.atPanelBegin))
                return;
        candidates.push_back(std::move(sample));
    };
    for (int index = 0; index < count; ++index)
    {
        if (layout.boxDistance(context, index) > context.radius + context.boxTolerance)
            continue;
        const auto samples = layout.samples(context, index);
        for (const auto &sample : samples)
        {
            if (sample.result.endpoint)
                addEndpoint(sample);
            else if (sample.derivative == 0.0)
                candidates.push_back(sample);
        }
        for (std::size_t k = 1; k < samples.size(); ++k)
            if (samples[k - 1].derivative < 0.0 && samples[k].derivative > 0.0)
                candidates.push_back(layout.minimum(context, index, samples[k - 1], samples[k]));
    }
    require(!candidates.empty(), "bounded local graph search found no stationary or endpoint KKT projection");
    const auto chosen =
        std::min_element(candidates.begin(), candidates.end(),
                         [](const auto &a, const auto &b) { return a.result.distance < b.result.distance; });
    const auto &result = chosen->result;
    context.accept(result);
    guardBranches(context, layout, result);
    return result;
}
} // namespace LocalInterfaceProjection
} // namespace
} // namespace Quantum::Geometry::Detail
