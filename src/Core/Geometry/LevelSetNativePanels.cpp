#include "Core/Geometry/LevelSetInterfacePanels.h"

#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Geometry/Detail/LevelSetInterfaceQuadrature.h"
#include <array>
#include <algorithm>
#include <exception>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
namespace Q = Quantum::Geometry::Detail::InterfaceQuadrature;
std::array<double, 2> measures(const std::vector<LevelSetInterfacePanels::ArcSample> &samples)
{
    long double length = 0, area = 0;
    for (const auto &sample : samples)
    {
        length += sample.weight;
        area += sample.physicalWeight;
    }
    return {{Q::checked(length, "native material arclength overflow"),
             Q::checked(area, "native material surface measure overflow")}};
}

// Positive composite GL8 on the actual retained graph. Compare a rule with
// its two halves for BOTH ds and 2*pi*r*ds, then accumulate accepted leaves
// in graph order. Prefixes are local to the same accepted leaves; independent
// full-panel prefixes can exceed a poorly resolved full-panel GL8 total.
// This changes quadrature only, never the scalar, roots, graph or liquid C.
void integratePanel(const LevelSetInterfacePanels::NativeGraphChain& chain,int index,
    std::vector<LevelSetInterfacePanels::ArcSample> initial,LevelSetInterfacePanels::Panel& panel)
{
    using Samples=std::vector<LevelSetInterfacePanels::ArcSample>;
    constexpr double relativeTolerance=1e-10;
    constexpr std::size_t maximumEvaluations=65536;
    std::size_t evaluations=initial.size();
    std::array<long double,2> accumulated{{0,0}},estimatedError{{0,0}};
    const auto rule=[&](double a,double b) {
        Q::require(evaluations+8<=maximumEvaluations,"native surface quadrature work budget exhausted");
        evaluations+=8;return chain.quadrature(index,a,b);
    };
    std::function<void(double,double,Samples,int)> visit;
    visit=[&](double a,double b,Samples coarse,int depth) {
        const double middle=a+(b-a)*.5;
        Q::require(middle>a && middle<b,"native surface quadrature interval is unrepresentable");
        auto left=rule(a,middle),right=rule(middle,b);
        const auto whole=measures(coarse),lo=measures(left),hi=measures(right);
        std::array<double,2> error;
        bool resolved=true;
        for(int m=0;m<2;++m) {
            error[m]=std::abs((static_cast<long double>(lo[m])+hi[m])-whole[m]);
            resolved&=whole[m]>0 && error[m]<=relativeTolerance*whole[m];
        }
        const bool accurate=resolved;
        std::array<double,2> lastPrevious{{0,0}},lastPartial{{0,0}};
        double lastFraction=0;
        std::vector<LevelSetInterfacePanels::QuadratureNode> nodes;
        if(resolved) {
            std::array<double,2> previous{{0,0}};
            for(std::size_t k=0;k<coarse.size();++k) {
                // The coordinate is recovered from the saved chart, so the
                // charge lookup sees precisely the same retained parameter.
                const auto& sample=coarse[k];
                const double fraction=chain.graphFraction(index,sample.frame.cellChartPoint);
                Q::require(fraction>a && fraction<b,"native surface quadrature node lost its interval");
                const auto partial=measures(rule(a,fraction));
                lastFraction=fraction;lastPrevious=previous;lastPartial=partial;
                for(int m=0;m<2;++m)
                    resolved&=partial[m]>previous[m] && partial[m]<whole[m];
                if(!resolved)break;
                nodes.push_back({sample.frame,sample.weight,partial[0],sample.physicalWeight,partial[1],fraction});
                previous=partial;
            }
        }
        if(!resolved) {
            if(depth>=20) {
                std::ostringstream detail;detail<<std::setprecision(17)
                    <<"native surface quadrature failed "<<(accurate?"cumulative ordering":"accuracy")
                    <<"; interval=["<<a<<','<<b<<"], depth="<<depth<<", evaluations="<<evaluations
                    <<", ds/area="<<whole[0]<<'/'<<whole[1]
                    <<", error="<<error[0]<<'/'<<error[1]
                    <<", budget="<<relativeTolerance*whole[0]<<'/'<<relativeTolerance*whole[1];
                if(accurate)detail<<", fraction="<<lastFraction<<", previous="<<lastPrevious[0]<<'/'<<lastPrevious[1]
                    <<", partial="<<lastPartial[0]<<'/'<<lastPartial[1];
                Q::require(false,detail.str().c_str());
            }
            visit(a,middle,std::move(left),depth+1);
            visit(middle,b,std::move(right),depth+1);
            return;
        }
        for(auto& node:nodes) {
            node.localArclength=Q::checked(accumulated[0]+node.localArclength,"native cumulative arclength overflow");
            node.localPhysicalMeasure=Q::checked(accumulated[1]+node.localPhysicalMeasure,"native cumulative area overflow");
            if(!panel.quadratureNodes.empty()) {
                const auto& previous=panel.quadratureNodes.back();
                Q::require(node.nativeGraphFraction>previous.nativeGraphFraction
                    && node.localArclength>previous.localArclength
                    && node.localPhysicalMeasure>previous.localPhysicalMeasure,
                    "native cumulative surface measure is unrepresentable");
            }
            panel.quadratureNodes.push_back(std::move(node));
        }
        for(int m=0;m<2;++m){accumulated[m]+=whole[m];estimatedError[m]+=error[m];}
        ++panel.measureQuadratureIntervals;
    };
    visit(0,1,std::move(initial),0);
    panel.measure=Q::checked(accumulated[0],"native panel arclength overflow");
    panel.physicalMeasure=Q::checked(accumulated[1],"native panel area overflow");
    Q::require(!panel.quadratureNodes.empty()
        && panel.measure>panel.quadratureNodes.back().localArclength
        && panel.physicalMeasure>panel.quadratureNodes.back().localPhysicalMeasure,
        "native cumulative endpoint measure is unrepresentable");
    for(int m=0;m<2;++m)
        panel.measureQuadratureRelativeError[m]=Q::checked(estimatedError[m]/accumulated[m],
            "native surface quadrature error estimate overflow");
}
} // namespace

LevelSetInterfacePanels::LevelSetInterfacePanels(NativeGraphChain chain, NativeMetadata)
    : cellToPanel_(chain.mapping_), axisymmetric_(chain.axisymmetric()), closed_(chain.closed()),
      beginEndpointKind_(chain.beginEndpointKind()), endEndpointKind_(chain.endEndpointKind()),
      boundaryContacts_(chain.boundaryContacts()), source_(chain.source()), nativeChain_(std::move(chain))
{
    Q::require(source_ && source_->usesNativeLeaves() && nativeChain_->graphs().size() >= 3,
               "native material panels require a retained complete interface chain");
}

const std::vector<int>& LevelSetInterfacePanels::NativeGraphChain::cellToGraph() const
{
    if (std::find(mapping_.begin(), mapping_.end(), -2) != mapping_.end())
        throw std::logic_error("single-graph consumer cannot use several arcs in one cell");
    return mapping_;
}

bool LevelSetInterfacePanels::NativeGraphChain::zeroSetPartialResolved(
    int index,int axis,double first,double last) const
{
    Q::require(axis>=0 && axis<2 && std::isfinite(first) && std::isfinite(last)
        && first>=0 && first<=1 && last>=0 && last<=1 && first!=last,
        "native partial certificate requires a nonempty retained subarc");
    const auto& graph=graphs_.at(index);
    const auto a=frameAt(index,first),b=frameAt(index,last);
    Q::require(a.hasCellChartPoint && b.hasCellChartPoint,"native subarc lost its local coordinates");
    const auto ca=graph.coordinateIsX?a.cellChartPoint.x:a.cellChartPoint.y;
    const auto cb=graph.coordinateIsX?b.cellChartPoint.x:b.cellChartPoint.y;
    auto lower=graph.localLower,upper=graph.localUpper;
    double& lo=graph.coordinateIsX?lower.x:lower.y;
    double& hi=graph.coordinateIsX?upper.x:upper.y;
    // Exact comparisons of endpoint-relative coordinates avoid rounding
    // 1-tiny to 1. Bracket each end independently: the smallest common
    // dyadic ancestor can still include a pole outside the requested arc.
    const double originalLo=lo,originalHi=hi;
    for(bool lowerEnd:{true,false}) {
        double left=originalLo,right=originalHi;
        for(int depth=0;depth<52;++depth) {
            const double middle=left+(right-left)/2;
            if(!(middle>left && middle<right))break;
            const Quantum::Geometry::CellChartCoordinate cut{middle,false};
            const int aSide=Quantum::Geometry::compareCellChartCoordinate(ca,cut);
            const int bSide=Quantum::Geometry::compareCellChartCoordinate(cb,cut);
            if(lowerEnd) {
                if(aSide>0 && bSide>0)left=middle;else right=middle;
            } else {
                if(aSide<0 && bSide<0)right=middle;else left=middle;
            }
        }
        if(lowerEnd)lo=left;else hi=right;
    }
    return source_->nativeZeroSetPartialResolved(graph.nativeCellIndex,axis,lower,upper);
}

const std::vector<int>& LevelSetInterfacePanels::cellToPanel() const
{
    if (std::find(cellToPanel_.begin(), cellToPanel_.end(), -2) != cellToPanel_.end())
        throw std::logic_error("single-panel consumer cannot use several arcs in one cell");
    return cellToPanel_;
}

LevelSetInterfacePanels::LevelSetInterfacePanels(NativeGraphChain chain)
    : LevelSetInterfacePanels(std::move(chain), NativeMetadata{})
{
    const auto &native = *nativeChain_;
    const auto count = native.graphs().size();
    panels_.resize(count);
    std::vector<std::array<LevelSetHermiteGeometry::Frame, 10>> fixedFrames(count);
    const auto preparePanel = [&](int index) {
        const auto &graph = native.graphs()[index];
        try
        {
        Panel panel;
        panel.cellIndex = graph.nativeCellIndex;
        panel.coordinateIsX = graph.coordinateIsX;
        panel.beginFrame = native.frameAt(index, 0);
        panel.endFrame = native.frameAt(index, 1);
        panel.frame = native.frameAt(index, .5);
        panel.begin = panel.beginFrame.point;
        panel.end = panel.endFrame.point;
        auto samples = native.quadrature(index);
        Q::require(samples.size() == Q::gaussX.size(), "native full-arc GL8 inventory changed");
        std::array<LevelSetHermiteGeometry::Frame, 10> fixed;
        fixed[0] = panel.beginFrame;
        fixed[1] = panel.endFrame;
        for (std::size_t k = 0; k < samples.size(); ++k)
            fixed[k + 2] = samples[k].frame;
        fixedFrames[index] = std::move(fixed);
        integratePanel(native,index,std::move(samples),panel);
        panels_[index] = std::move(panel);
        }
        catch (const std::exception &dbgError)
        {
            std::ostringstream dbg;
            dbg << std::setprecision(17)
                << "[LevelSetNativePanels P0.1 diagnostic] graphIndex=" << index
                << " of " << count << ", nativeCellIndex=" << graph.nativeCellIndex
                << ", coordinateIsX=" << (graph.coordinateIsX ? 1 : 0)
                << ", chart lower=(" << graph.localLower.x << ',' << graph.localLower.y
                << "), upper=(" << graph.localUpper.x << ',' << graph.localUpper.y
                << "), begin(side=" << graph.begin.side << ",ord=" << graph.begin.rootOrdinal
                << "), end(side=" << graph.end.side << ",ord=" << graph.end.rootOrdinal << ')';
            try
            {
                const auto b = source_->nativeCellBounds(graph.nativeCellIndex);
                dbg << " | bounds=(" << b[0].x << ',' << b[0].y << ")-("
                    << b[1].x << ',' << b[1].y << ')';
                if(source_->reconstruction()==LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic) {
                    const auto dofs=source_->nativeQuarticLocalDofs(graph.nativeCellIndex);
                    dbg<<" | Q4_local_DOFs=[";
                    for(double v:dofs)dbg<<v<<',';
                    dbg<<']';
                } else {
                    // Displayed physical jets, not an exact serialized source.
                    // Keep cubic failures useful too: the Q4-only query used
                    // to throw before even the physical box was recorded.
                    dbg<<" | rounded_Q3_corner_jets_SW_SE_NE_NW=[";
                    for(const auto p:{b[0],Vector2D{b[1].x,b[0].y},b[1],Vector2D{b[0].x,b[1].y}}) {
                        const auto jet=source_->scalarJetInNativeCell(graph.nativeCellIndex,p);
                        dbg<<jet.value<<','<<jet.dx<<','<<jet.dy<<','<<jet.dxy<<';';
                    }
                    dbg<<']';
                }
            }
            catch (...) {}
            dbg << ": " << dbgError.what();
            throw std::runtime_error(dbg.str());
        }
    };
    // One immutable scalar/graph inventory, independent panel slots. Keep
    // each panel's adaptive integration and summation order independent. The
    // first panel also checks the caller's floating environment before work
    // is dispatched; worker failures are joined in original chain order.
    preparePanel(0);
    std::vector<std::exception_ptr> failures(count);
#pragma omp parallel for schedule(dynamic, 1) if(count > 8)
    for (int index = 1; index < static_cast<int>(count); ++index)
    {
        try { preparePanel(index); }
        catch (...) { failures[index] = std::current_exception(); }
    }
    for (const auto& failure : failures)
        if (failure) std::rethrow_exception(failure);
    // Publish only the complete immutable inventory. A failed constructor
    // never leaves a partially prepared chain visible to a consumer.
    nativeChain_->fixedProjectionFrames_ = std::move(fixedFrames);
    nativeChain_->prepareProjectionCoordinates();
}

const LevelSetInterfacePanels::NativeGraphChain &LevelSetInterfacePanels::nativeChain() const
{
    Q::require(nativeChain_.has_value(), "native material query requires explicit native construction");
    return *nativeChain_;
}

LevelSetInterfacePanels LevelSetInterfacePanels::rebindNative(NativeGraphChain target,
    const std::vector<int>& cellMap) const
{
    namespace Exact = Quantum::Geometry::Detail::ExactExpansion;
    Exact::validateEnvironment();
    using Reconstruction = LevelSetHermiteGeometry::Reconstruction;
    Q::require(nativeChain_ && source_ && target.source() &&
        source_->reconstruction() == Reconstruction::EnrichedQuartic &&
        target.source()->reconstruction() == Reconstruction::EnrichedQuartic &&
        cellMap.size() == source_->nativeCellCount(), "native panel rebind requires complete Q4 source identities");
    Q::require(target.graphs().size() == panels_.size() && target.closed() == closed_ &&
        target.axisymmetric() == axisymmetric_ && target.beginEndpointKind() == beginEndpointKind_ &&
        target.endEndpointKind() == endEndpointKind_ && target.boundaryContacts().size() == boundaryContacts_.size(),
        "native panel rebind changed material topology or endpoint kinds");
    const auto samePoint = [](Vector2D a, Vector2D b) { return a.x == b.x && a.y == b.y; };
    Q::require(samePoint(source_->minCoord(), target.source()->minCoord()) &&
        samePoint(source_->maxCoord(), target.source()->maxCoord()), "native panel rebind changed physical domain");
    for (std::size_t k = 0; k < boundaryContacts_.size(); ++k)
        Q::require(samePoint(boundaryContacts_[k], target.boundaryContacts()[k]),
                   "native panel rebind changed a declared contact");
    Q::require(target.symmetryPlanes()==nativeChain_->symmetryPlanes() &&
        target.endpointSymmetryPlanes_==nativeChain_->endpointSymmetryPlanes_,
        "native panel rebind changed a symmetry section identity");
    const auto count = panels_.size();
    std::vector<int> targetBySource(count, -1), sourceByTarget(count, -1);
    bool identicalRepresentation = true;
    for (std::size_t i = 0; i < count; ++i)
    {
        const auto& old = panels_[i];
        const int cell = cellMap.at(old.cellIndex);
        Q::require(cell >= 0 && static_cast<std::size_t>(cell) < target.cellToGraph().size(),
                   "retained panel has no target leaf");
        const int j = target.cellToGraph()[cell];
        Q::require(j >= 0 && static_cast<std::size_t>(j) < count && sourceByTarget[j] < 0,
                   "retained panel map is not a bijection");
        const auto& a = nativeChain_->graphs()[i]; const auto& b = target.graphs()[j];
        Q::require(a.coordinateIsX == b.coordinateIsX && a.begin.side == b.begin.side &&
            samePoint(a.localLower,b.localLower) && samePoint(a.localUpper,b.localUpper) &&
            a.controlledIntegration==b.controlledIntegration &&
            a.end.side == b.end.side && a.begin.rootOrdinal == b.begin.rootOrdinal &&
            a.end.rootOrdinal == b.end.rootOrdinal, "retained panel changed orientation or original edge-root identity");
        Q::require(source_->identicalNativeQuarticCell(old.cellIndex, *target.source(), cell),
                   "retained panel polynomial or physical bounds changed");
        // Equal rational polynomials can have different expansion histories.
        // Reusing floating root representatives additionally needs identical
        // restriction inputs and normalized metric jets, not equality alone.
        identicalRepresentation &= source_->identicalNativeQuarticRepresentation(
            old.cellIndex, *target.source(), cell);
        targetBySource[i] = j; sourceByTarget[j] = static_cast<int>(i);
    }
    for (std::size_t i = 0; i < count; ++i)
        Q::require(closed_ ? targetBySource[(i + 1) % count] == (targetBySource[i] + 1) % count
                           : targetBySource[i] == static_cast<int>(i),
                   "retained panel map changed neighbors or open endpoint order");
    if (!identicalRepresentation) return LevelSetInterfacePanels(std::move(target));

    LevelSetInterfacePanels result(std::move(target), NativeMetadata{});
    result.panels_.resize(count);
    std::vector<std::array<LevelSetHermiteGeometry::Frame, 10>> frames(count);
    for (std::size_t j = 0; j < count; ++j)
    {
        const auto& old = panels_[sourceByTarget[j]];
        auto panel = old;
        const int cell = result.nativeChain_->graphs()[j].nativeCellIndex;
        const auto reindex = [&](LevelSetHermiteGeometry::Frame& frame) {
            Q::require(frame.nativeCellIndex == old.cellIndex && frame.hasCellChartPoint &&
                !frame.hasCellLocalPoint && frame.cellI == -1 && frame.cellJ == -1,
                "retained panel frame lost native provenance");
            frame.nativeCellIndex = cell;
        };
        panel.cellIndex = cell;
        reindex(panel.beginFrame); reindex(panel.endFrame); reindex(panel.frame);
        // Projection's fixed endpoint/GL8 search inventory is independent of
        // adaptive material integration nodes. Rebind BOTH unchanged sets.
        frames[j]=nativeChain_->fixedProjectionFrames_.at(sourceByTarget[j]);
        for(auto& frame:frames[j])reindex(frame);
        for(auto& node:panel.quadratureNodes)reindex(node.frame);
        result.panels_[j] = std::move(panel);
    }
    result.nativeChain_->fixedProjectionFrames_ = std::move(frames);
    result.nativeChain_->prepareProjectionCoordinates();
    result.nativeMeasuresReused_ = true;
    return result;
}
