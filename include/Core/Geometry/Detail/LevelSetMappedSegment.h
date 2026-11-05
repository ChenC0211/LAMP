#pragma once

#include "Core/Geometry/Detail/LevelSetMappedCellLayout.h"
#include "Core/Geometry/Detail/LevelSetMappedCoordinates.h"
#include "Core/Geometry/Detail/LevelSetMappedRootRefinement.h"
#include "Core/Geometry/Detail/PolynomialSquareFree.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include <algorithm>
#include <functional>
#include <iomanip>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace LAMP::Geometry::Detail
{
namespace
{
namespace MappedHermiteSegment
{
struct CardinalQuarticLayout;
struct Phase { Ratio begin, end; bool negative; };
struct ZeroTrace {Chart chart;std::vector<CellAddress> cells;};
struct RootRecord
{
    LevelSetMappedSegment::Root info;
    Ratio parameter, lower, upper; // Original canonical input-segment chart.
};

template <class Layout> class BasicQuery
{
public:
    using Polynomial = typename Layout::Polynomial;
    Vector2D inputBegin, inputEnd;
    std::vector<LevelSetMappedSegment::Range> ranges;
    std::vector<LevelSetMappedSegment::Root> rootInfo;
    std::vector<Chart> charts;
    std::vector<RootRecord> roots;
    std::vector<Phase> phases;
    std::size_t gridPieces = 0;
    bool hasMultipleRootCharts = false;
    std::vector<ZeroTrace> zeroTraces;

    BasicQuery(const typename Layout::Source& source, Vector2D begin, Vector2D end)
        : BasicQuery(source,begin,end,Layout(source),false) {}
    BasicQuery(const typename Layout::Source& source,Vector2D begin,Vector2D end,
        Layout layout,bool permitZero=false)
        : inputBegin(begin),inputEnd(end),layout_(std::move(layout)),a_(begin),b_(end),permitZero_(permitZero)
    {
        Exact::validateEnvironment();
        const auto lower = source.minCoord(), upper = source.maxCoord();
        const auto valid = [&](Vector2D p) {
            return std::isfinite(p.x) && std::isfinite(p.y) && p.x >= lower.x && p.x <= upper.x
                && p.y >= lower.y && p.y <= upper.y;
        };
        require(valid(begin) && valid(end), "segment is outside the finite snapshot rectangle");
        require(begin.x != end.x || begin.y != end.y, "zero-length segment has no phase partition");
        const bool reverse = std::make_pair(b_.x,b_.y) < std::make_pair(a_.x,a_.y);
        if (reverse) std::swap(a_, b_);

        auto cuts = layout_.cuts(a_, b_);
        std::sort(cuts.begin(), cuts.end(), [](const Ratio& x, const Ratio& y) { return compare(x,y) < 0; });
        cuts.erase(std::unique(cuts.begin(), cuts.end(),
            [](const Ratio& x, const Ratio& y) { return compare(x,y) == 0; }), cuts.end());
        gridPieces = cuts.size()-1;
        if (gridPieces == 1 && retainCertifiedSinglePhase(cuts.front(), cuts.back(), reverse)) return;
        for (std::size_t k = 1; k < cuts.size(); ++k)
        {
            try { partition(cuts[k-1], cuts[k]); }
            catch (const std::exception& error)
            {
                throw std::runtime_error(std::string("[LevelSetMappedSegment] phase partition: ") + error.what());
            }
        }

        // All closed-cell endpoint incidences, including a cell meeting ONLY
        // a grid-cut point. Keep identity; do not assign such a cell a phase arc.
        for (const auto& cut : cuts)
            for (const auto& cell : incidentCells(cut))
            {
                if (layout_.pieceSign(cell,cut,cut) != 0) continue;
                const auto p = polynomial(cell,cut,cut);
                if (Exact::sign(p[0]) != 0) continue;
                // Positive-length zero traces have typed side events instead
                // of a fictitious finite inventory of value roots.
                bool zeroIncidence=false;
                for(const auto& z:zeroTraces)if(compare(z.chart.begin,cut)<=0 && compare(cut,z.chart.end)<=0
                    && std::find(z.cells.begin(),z.cells.end(),cell)!=z.cells.end())zeroIncidence=true;
                if(zeroIncidence)continue;
                const std::size_t chart = charts.size();
                charts.push_back({cut,cut});
                retainRoot(cell,chart,0.0,0.0,0.0,true);
            }

        finish();
    }

    void refineRoot(std::size_t publicIndex)
    {
        require(publicIndex < roots.size(), "root refinement index is out of range");
        const bool reverse = std::make_pair(inputEnd.x,inputEnd.y) < std::make_pair(inputBegin.x,inputBegin.y);
        const RootRecord parent = roots[reverse ? roots.size()-1-publicIndex : publicIndex];
        require(!parent.info.pointOnly && compare(parent.lower,parent.upper)<0,
            "an exact point root needs no refinement");
        const auto child=MappedHermiteRootRefinement::refine(parent.lower,parent.upper,
            [&](const Ratio& a,const Ratio& b) {
                auto p=polynomial({parent.info.cellI,parent.info.cellJ,parent.info.nativeCellIndex},a,b);
                if constexpr(std::is_same_v<Layout,NativeQuarticCellLayout> || std::is_same_v<Layout,NativeC2MacroCellLayout> || std::is_same_v<Layout,CardinalQuarticLayout>)
                    if(hasMultipleRootCharts)p=PolynomialSquareFree::squareFree(p).polynomial;
                // The simple-root refinement certificate belongs to the
                // square-free factor here; only its same-root bracket is used.
                return p;
            });
        const std::size_t chart=charts.size();
        charts.push_back(child.chart);
        // A shared-side incidence is not a second geometric root. Update all
        // incidences and both adjacent phase endpoints in this transaction.
        for (auto& root:roots)
            if (compare(root.parameter,parent.parameter)==0)
            {
                require(compare(root.lower,child.originalUpper)<=0 && compare(child.originalLower,root.upper)<=0,
                    "same-point root incidences have disjoint refined enclosures");
                root.info.chart=chart;
                root.info.parameter=child.parameter;
                root.info.lower=child.lower;
                root.info.upper=child.upper;
                root.parameter=child.originalParameter;
                root.lower=child.originalLower;
                root.upper=child.originalUpper;
            }
        for (auto& phase:phases)
        {
            if (compare(phase.begin,parent.parameter)==0) phase.begin=child.originalParameter;
            if (compare(phase.end,parent.parameter)==0) phase.end=child.originalParameter;
        }
        finish();
    }

    CellAddress cellAtCanonicalParameter(const Ratio& parameter) const
    {
        const auto cells=layout_.incidentCells(parameter);
        require(!cells.empty(),"mapped trace sample has no active owner");
        return cells.front();
    }
    Polynomial completeTracePolynomial() const
    {
        require(gridPieces==1,"complete native trace requires one positive-length chart");
        return layout_.polynomial(cellAtCanonicalParameter(constant(.5)),constant(0),constant(1));
    }
    bool proportionalTrace(const BasicQuery& other) const
    {
        const auto p=completeTracePolynomial(),q=other.completeTracePolynomial();
        int pivot=-1;
        for (int k=0;k<static_cast<int>(p.size());++k) if (Exact::sign(p[k])!=0) { pivot=k; break; }
        require(pivot>=0,"zero polynomial has no trace proportionality");
        return Exact::sign(q[pivot])!=0 && layout_.sameTrace(p,q,pivot);
    }
private:
    void finish()
    {
        ranges.clear(); rootInfo.clear();
        const bool reverse = std::make_pair(inputEnd.x,inputEnd.y) < std::make_pair(inputBegin.x,inputBegin.y);
        std::sort(phases.begin(), phases.end(), [](const Phase& x, const Phase& y) {
            return compare(x.begin,y.begin) < 0;
        });
        std::vector<Chart> coverage;
        Ratio covered=constant(0.0);
        for (const auto& phase : phases)
        {
            if(zeroTraces.empty()){
                require(compare(covered,phase.begin)==0 && compare(phase.begin,phase.end)<0,
                    "phase partition has a gap, overlap or collapsed interval");covered=phase.end;
            }else{
                require(compare(phase.begin,phase.end)<0,"phase interval collapsed");
                coverage.push_back({phase.begin,phase.end});
            }
            const bool fromEnd = compare(add(phase.begin,phase.end),constant(1.0)) > 0;
            const Ratio offset = fromEnd ? subtract(constant(1.0),phase.end) : phase.begin;
            const Ratio width = subtract(phase.end,phase.begin);
            const long double o = approximate(offset), w = approximate(width);
            require(o >= 0.0L && o <= 1.0L && w > 0.0L && w <= 1.0L,
                "invalid terminal integration offset or width");
            ranges.push_back({fromEnd != reverse,o,w,phase.negative});
        }
        if(!zeroTraces.empty()){
            for(const auto& z:zeroTraces)coverage.push_back(z.chart);
            std::sort(coverage.begin(),coverage.end(),[](const Chart& x,const Chart& y){return compare(x.begin,y.begin)<0;});
            for(const auto& c:coverage){require(compare(covered,c.begin)==0 && compare(c.begin,c.end)<0,
                "phase partition has a gap, overlap or collapsed interval");covered=c.end;}
        }
        require(compare(covered,constant(1.0)) == 0, "phase partition does not cover the whole segment");
        std::stable_sort(roots.begin(), roots.end(), [](const RootRecord& x, const RootRecord& y) {
            return compare(x.parameter,y.parameter) < 0;
        });
        for (const auto& root : roots) rootInfo.push_back(root.info);
        if (reverse)
        {
            std::reverse(ranges.begin(),ranges.end());
            std::reverse(rootInfo.begin(),rootInfo.end());
        }
    }

    Layout layout_;
    Vector2D a_, b_;
    bool permitZero_=false;
    std::vector<CellAddress> incidentCellsAtInputPoint(Vector2D point) const
    { return layout_.incidentCellsAtInputPoint(point); }
    std::vector<CellAddress> incidentCells(const Ratio &parameter) const
    { return layout_.incidentCells(parameter); }
    Polynomial polynomial(const CellAddress &cell, const Ratio &begin, const Ratio &end) const
    {
        try { return layout_.polynomial(cell,begin,end); }
        catch (const std::exception &error)
        {
            std::ostringstream message;
            message << "[LevelSetMappedSegment] cleared polynomial cell=(" << cell << "): " << error.what();
            throw std::runtime_error(message.str());
        }
    }
    bool retainCertifiedSinglePhase(const Ratio& begin, const Ratio& end, bool reverse)
    {
        // With no interior grid cut, an active closed cell containing BOTH
        // endpoints covers the entire segment. Require strict whole-cell
        // certificates on EVERY endpoint incidence too: a zero in a cell
        // meeting only an endpoint must still use the full root inventory.
        // Never substitute a center/node sample, distance threshold or root
        // tolerance for these existing immutable scalar certificates.
        const auto first = incidentCellsAtInputPoint(a_), last = incidentCellsAtInputPoint(b_);
        int phaseSign = 0;
        for (const auto* cells : {&first, &last})
            for (const auto& cell : *cells)
            {
                const int sign = layout_.strictSign(cell);
                if (sign == 0 || (phaseSign != 0 && sign != phaseSign)) return false;
                phaseSign = sign;
            }
        bool covered = false;
        for (const auto& cell : first)
            covered = covered || std::find(last.begin(),last.end(),cell) != last.end();
        if (phaseSign == 0 || !covered) return false;
        phases.push_back({begin,end,phaseSign < 0});
        // Exactly the finish() result for the unchanged single [0,1] phase.
        // Keep reversal's fromEnd flag, not just an equivalent integral.
        ranges.push_back({reverse,0.0L,1.0L,phaseSign < 0});
        return true;
    }

    void retainRoot(const CellAddress &cell, std::size_t chart, double parameter, double lower, double upper, bool pointOnly)
    {
        const auto& c = charts[chart];
        const Ratio p = between(c.begin,c.end,parameter);
        // Shared cut endpoints can be visited from two grid pieces. Only
        // identical same-cell incidences are deduplicated; other cells remain.
        for (const auto& root : roots)
            if (root.info.cellI == cell.i && root.info.cellJ == cell.j && root.info.nativeCellIndex == cell.native && compare(root.parameter,p) == 0) return;
        auto lo = between(c.begin,c.end,lower), hi = between(c.begin,c.end,upper);
        if (compare(lo,hi) > 0) std::swap(lo,hi);
        roots.push_back({{cell.i,cell.j,chart,parameter,lower,upper,pointOnly,cell.native},p,std::move(lo),std::move(hi)});
    }
    // P0.1 bounded dyadic fallback for a one-shot near-critical failure on ONE
    // native Q4 leaf piece. When an interface arc runs almost collinear with a
    // physical face, the degree-four along-face trace can carry a tiny extremum
    // that does NOT cross zero, yet makes the single-piece interval isolator
    // report an interior near-critical root (Zalesak cell 3422, rotated slot wall
    // on its west face). Re-restricting the SAME retained scalar to dyadic
    // subpieces is exactly a Bernstein de Casteljau subdivision: it only
    // tightens the coefficient hull; it rounds no coefficient and neither invents
    // nor deletes a root. A pure-phase subpiece isolates with no roots; a monotone
    // transverse-root subpiece isolates cleanly. A genuine interior
    // multiple/tangential root not resolved by the native Q4 algebraic
    // certificate stays critical at every depth and is rethrown once
    // the bounded budget is exhausted. Cubic traces keep the original one-shot
    // behaviour (tuple_size 7), so this recovery is quartic-only.
    PolynomialSegmentRoots::Result isolatePieceWithDyadicFallback(
        const CellAddress& owner,const Chart& chart,const Polynomial& p,bool& repeated) const
    {
        using PolynomialSegmentRoots::isolateEndpointRelative;
        try { return isolateEndpointRelative(p); }
        catch (const std::exception& outer)
        {
            if constexpr(std::is_same_v<Layout,NativeQuarticCellLayout> || std::is_same_v<Layout,NativeC2MacroCellLayout> || std::is_same_v<Layout,CardinalQuarticLayout>) {
                // Exact algebraic deflation supplies distinct zeros only.
                // If bounded expansion arithmetic cannot certify a factor,
                // retain the existing same-polynomial subdivision below.
                std::optional<Polynomial> squareFree;
                try {const auto factor=PolynomialSquareFree::squareFree(p);if(factor.deflated)squareFree=factor.polynomial;}
                catch(const std::runtime_error&) {}
                if(squareFree){auto result=isolateEndpointRelative(*squareFree);repeated=true;return result;}
            }
            if constexpr (std::tuple_size<Polynomial>::value != 9) throw;
            constexpr int maxDepth=6;
            PolynomialSegmentRoots::Result merged;
            std::vector<bool> exactSplitRoots;
            bool failed=false;
            std::function<void(double,double,int)> recurse=[&](double lo,double hi,int depth)
            {
                if(failed) return;
                const Ratio a=between(chart.begin,chart.end,lo);
                const Ratio b=between(chart.begin,chart.end,hi);
                Polynomial q=polynomial(owner,a,b);
                try
                {
                    auto sub=isolateEndpointRelative(q);
                    if(sub.identicallyZero) { failed=true; return; }
                    for(std::size_t k=0;k<sub.roots.size();++k)
                    {
                        const bool endpoint=sub.roots[k]==0 || sub.roots[k]==1;
                        if (endpoint && ((sub.roots[k]==0 && lo>0) || (sub.roots[k]==1 && hi<1)))
                            require(MappedHermiteRootRefinement::endpointDerivativeSign(q,sub.roots[k]==1)!=0,
                                "dyadic subdivision exposed a multiple root; physical refinement is required");
                        const auto mapped=[&](double t) {
                            return ExactExpansion::add(ExactExpansion::fromDouble(lo),
                                ExactExpansion::multiply(ExactExpansion::fromDouble(hi-lo),ExactExpansion::fromDouble(t)));
                        };
                        merged.roots.push_back(lo+(hi-lo)*sub.roots[k]);
                        merged.enclosures.push_back({std::max(lo,ExactExpansion::enclosure(mapped(sub.enclosures[k].lo)).lo),
                            std::min(hi,ExactExpansion::enclosure(mapped(sub.enclosures[k].hi)).hi)});
                        exactSplitRoots.push_back(endpoint);
                    }
                }
                catch(const std::exception&)
                {
                    if(depth>=maxDepth) { failed=true; return; }
                    const double mid=0.5*(lo+hi);
                    recurse(lo,mid,depth+1);
                    recurse(mid,hi,depth+1);
                }
            };
            recurse(0.0,1.0,0);
            if(failed)
                throw std::runtime_error(std::string("dyadic root-isolation budget exhausted: ")+outer.what());
            std::vector<std::size_t> order(merged.roots.size());
            std::iota(order.begin(),order.end(),std::size_t{0});
            std::sort(order.begin(),order.end(),
                [&](std::size_t a,std::size_t b){return merged.roots[a]<merged.roots[b];});
            PolynomialSegmentRoots::Result sorted;
            bool previousExact=false;
            for(std::size_t idx:order)
            {
                const double r=merged.roots[idx];
                auto e=merged.enclosures[idx];
                if(!sorted.enclosures.empty() && e.lo<=sorted.enclosures.back().hi)
                {
                    // Only a shared EXACT dyadic endpoint proves root identity.
                    // Overlapping numerical brackets can contain distinct roots.
                    require(previousExact && exactSplitRoots[idx] && r==sorted.roots.back(),
                        "distinct dyadic root enclosures overlap; no root was merged");
                    continue;
                }
                sorted.roots.push_back(r);
                sorted.enclosures.push_back(e);
                previousExact=exactSplitRoots[idx];
            }
            sorted.identicallyZero=false;
            return sorted;
        }
    }

    void partition(const Ratio& begin, const Ratio& end)
    {
        const auto cells = incidentCells(between(begin,end,0.5));
        require(!cells.empty(), "positive-length segment has no active closed-cell coverage");
        const auto owner = cells.back();
        const int strictSign = layout_.pieceSign(owner,begin,end);
        if (strictSign != 0)
        {
            phases.push_back({begin,end,strictSign < 0});
            return;
        }
        Chart chart{begin,end};
        auto p = polynomial(owner,chart.begin,chart.end);
        if constexpr(std::is_same_v<Layout,NativeQuarticCellLayout> || std::is_same_v<Layout,NativeC2MacroCellLayout>)if(permitZero_){
            const auto zero=[](const Polynomial& q){return std::all_of(q.begin(),q.end(),
                [](const Expansion& v){return Exact::sign(v)==0;});};
            if(zero(p)){
                for(const auto& cell:cells)require(zero(polynomial(cell,chart.begin,chart.end)),
                    "adjacent native zero traces are not the same retained scalar");
                zeroTraces.push_back({chart,cells});return;
            }
        }
        const auto absolute = [](Expansion value) {
            return Exact::sign(value) < 0 ? Exact::subtract(Expansion{},value) : value;
        };
        if (Exact::sign(Exact::subtract(absolute(Exact::evaluate(p,1.0)),absolute(p[0]))) < 0)
        {
            std::swap(chart.begin,chart.end);
            p = polynomial(owner,chart.begin,chart.end);
        }
        const std::size_t chartIndex = charts.size();
        charts.push_back(chart);
        PolynomialSegmentRoots::Result isolated;
        bool repeated=false;
        try
        {
            isolated = isolatePieceWithDyadicFallback(owner, chart, p,repeated);
            hasMultipleRootCharts=hasMultipleRootCharts || (repeated && !isolated.roots.empty());
        }
        catch (const std::exception& error)
        {
            std::ostringstream message;
            const long double pb=approximate(begin),pe=approximate(end);
            message << std::setprecision(17) << "[LevelSetMappedSegment] root isolation cell=("
                << owner << "); segmentA=("<<a_.x<<","<<a_.y<<") segmentB=("<<b_.x<<","<<b_.y
                <<") piece0=("<<(a_.x+(b_.x-a_.x)*(double)pb)<<","<<(a_.y+(b_.y-a_.y)*(double)pb)
                <<") piece1=("<<(a_.x+(b_.x-a_.x)*(double)pe)<<","<<(a_.y+(b_.y-a_.y)*(double)pe)
                <<"); coefficients=";
            for (const auto& coefficient : p)
            {
                message << '[';
                for (double value : coefficient.components()) message << value << ',';
                message << ']';
            }
            message << ": " << error.what();
            throw std::runtime_error(message.str());
        }
        require(!isolated.identicallyZero, "positive-length coincident zero trace is unsupported");
        require(isolated.roots.size() == isolated.enclosures.size(), "root representative/enclosure count mismatch");

        // A line on a shared cell side has the same value trace, potentially
        // with another positive clearing factor. Prove that identity before
        // sharing the owner's root chart among its closed-cell incidences.
        int pivot = 0;
        while (pivot < static_cast<int>(p.size()) && Exact::sign(p[pivot]) == 0) ++pivot;
        require(pivot < static_cast<int>(p.size()), "nonzero restriction has no coefficient");
        for (const auto& cell : cells)
        {
            if (cell == owner) continue;
            const auto q = polynomial(cell,chart.begin,chart.end);
            require(Exact::sign(q[pivot]) == Exact::sign(p[pivot]), "adjacent mapped traces disagree in sign");
            try
            {
                require(layout_.sameTrace(p,q,pivot),
                    "adjacent mapped cell traces are not the same retained scalar");
            }
            catch (const std::exception& error)
            {
                throw std::runtime_error(std::string("[LevelSetMappedSegment] shared trace identity: ") + error.what());
            }
        }
        // The first nonzero ascending coefficient fixes the one-sided sign
        // immediately after s=0, including an exact endpoint zero. All
        // Ordinary roots below are simple. An exact square-free inventory
        // retains every distinct zero; its phases come from the ORIGINAL P.
        int phaseSign = Exact::sign(p[pivot]);
        Ratio previous = chart.begin;
        double previousIsolated = 0.0;
        const auto append = [&](Ratio lo, Ratio hi, int sign) {
            if (compare(lo,hi) == 0) return;
            if (compare(lo,hi) > 0) std::swap(lo,hi);
            require(compare(lo,hi) < 0, "root interval collapsed in its retained chart");
            phases.push_back({std::move(lo),std::move(hi),sign < 0});
        };
        for (std::size_t k = 0; k < isolated.roots.size(); ++k)
        {
            const double root = isolated.roots[k];
            const auto bracket = isolated.enclosures[k];
            require(root >= 0.0 && root <= 1.0 && bracket.lo >= 0.0 && bracket.hi <= 1.0
                && bracket.lo <= root && root <= bracket.hi, "invalid isolated root enclosure");
            require((root != 0.0 || Exact::sign(p[0]) == 0)
                && (root != 1.0 || Exact::sign(Exact::evaluate(p,1.0)) == 0),
                "interior mapped root collapsed onto a nonzero endpoint");
            for (const auto& cell : cells)
                retainRoot(cell,chartIndex,root,bracket.lo,bracket.hi,false);
            if (root == 0.0 || root == 1.0) continue;
            require(root > previousIsolated, "isolated root representatives are not strictly ordered");
            const Ratio parameter = between(chart.begin,chart.end,root);
            // Phase moments, phase fluxes and free-arc cuts retain exactly
            // the SAME refined point; no consumer keeps the coarse double.
            append(previous,parameter,phaseSign);
            previous = parameter;
            previousIsolated = root;
            if(!repeated)phaseSign = -phaseSign;
            else {
                const double next=k+1<isolated.roots.size()?isolated.enclosures[k+1].lo:1.;
                const double sample=.5*bracket.hi+.5*next;
                require(sample>bracket.hi && sample<next,"no representable phase witness between distinct roots");
                phaseSign=Exact::sign(Exact::evaluate(p,sample));
                require(phaseSign!=0,"square-free phase witness is an original-polynomial zero");
            }
        }
        append(previous,chart.end,phaseSign);
        const int lastSign = Exact::sign(Exact::evaluate(p,1.0));
        require(lastSign == 0 || lastSign == phaseSign, "endpoint phase disagrees with the complete root inventory");
    }
};
using Query = BasicQuery<UniformCellLayout>;
using QuarticQuery = BasicQuery<QuarticCellLayout>;
using NativeQuery = BasicQuery<NativeCellLayout>;
using NativeQuarticQuery = BasicQuery<NativeQuarticCellLayout>;
using NativeC2MacroQuery = BasicQuery<NativeC2MacroCellLayout>;
}
}
}
