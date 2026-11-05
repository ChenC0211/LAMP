#include "Modules/Electric/SurfaceChargeField.h"
#include "Core/Geometry/CellChartPoint.h"

#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
void require(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[SurfaceChargeMeshTransfer] ") + reason);
}
bool sameChart(const LevelSetHermiteGeometry::Frame& a, const LevelSetHermiteGeometry::Frame& b)
{
    return a.hasCellChartPoint && b.hasCellChartPoint
        && Geometry::cellChartEqual(a.cellChartPoint, b.cellChartPoint);
}
void validatePanel(const LevelSetInterfacePanels::Panel& old, const LevelSetInterfacePanels::Panel& next)
{
    require(old.coordinateIsX == next.coordinateIsX && sameChart(old.beginFrame,next.beginFrame)
        && sameChart(old.endFrame,next.endFrame), "retained material panel changed its orientation or endpoint chart");
    require(old.measure == next.measure && old.physicalMeasure == next.physicalMeasure
        && old.quadratureNodes.size() == next.quadratureNodes.size(), "retained material panel measure changed");
    for (std::size_t q = 0; q < old.quadratureNodes.size(); ++q)
    {
        const auto& a = old.quadratureNodes[q]; const auto& b = next.quadratureNodes[q];
        require(sameChart(a.frame,b.frame) && a.nativeGraphFraction == b.nativeGraphFraction
            && a.weight == b.weight && a.physicalWeight == b.physicalWeight
            && a.localArclength == b.localArclength && a.localPhysicalMeasure == b.localPhysicalMeasure,
            "retained material panel quadrature/moments changed");
    }
}
}

SurfaceChargeMeshTransfer SurfaceChargeField::rebindMesh(const QuadtreeMesh& mesh,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> nextSnapshot,
    const std::vector<int>& cellMap, const std::vector<std::vector<double>>& history,
    const LevelSetInterfacePanels* transferredPanels) const
{
    require(source_ && nextSnapshot && panels_.usesNativeLeaves()
        && source_->usesNativeLeaves() && nextSnapshot->usesNativeLeaves()
        && source_->reconstruction() == LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
        && nextSnapshot->reconstruction() == LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic,
        "mesh rebind requires two retained native Q4 scalar sources");
    require(cellMap.size() == source_->nativeCellCount(), "cell map must use the complete source native index array");
    const auto oldLow = source_->minCoord(), oldHigh = source_->maxCoord();
    const auto nextLow = nextSnapshot->minCoord(), nextHigh = nextSnapshot->maxCoord();
    require(oldLow.x == nextLow.x && oldLow.y == nextLow.y
        && oldHigh.x == nextHigh.x && oldHigh.y == nextHigh.y,
        "mesh rebind changed the physical background domain");
    // Validate every history before discovering/reconstructing candidate panels.
    for (const auto& values : history)
        SurfaceChargeTransport::fromIntegratedCharge(charge_.panelMeasure(), values, panels_.closed());
    const auto& oldPanels = panels_.panels();
    for (const auto& panel : oldPanels)
    {
        const int next = cellMap.at(panel.cellIndex);
        require(next >= 0 && static_cast<std::size_t>(next) < nextSnapshot->nativeCellCount(),
            "charged source leaf has no retained target leaf");
        require(source_->identicalNativeQuarticCell(panel.cellIndex,*nextSnapshot,next),
            "charged leaf polynomial or physical bounds changed; a copied index is not an identity certificate");
    }

    SurfaceChargeMeshTransfer result;
    if (transferredPanels)
    {
        require(transferredPanels->usesNativeLeaves() && transferredPanels->usesSource(*nextSnapshot)
            && transferredPanels->axisymmetric() == panels_.axisymmetric(), "transferred panels use another scalar or metric");
        // Interface may have prepared this exact immutable Core geometry.
        // All original chart, measure, topology and charge/history gates below
        // still apply; only redundant discovery/rebinding is avoided.
        result.field = std::make_unique<SurfaceChargeField>(mesh, *transferredPanels, 0.0);
    }
    else
    {
        auto chain = LevelSetInterfacePanels::discoverNativeGraphs(mesh, nextSnapshot,
            panels_.axisymmetric(), panels_.nativeChain().boundaryContacts());
        result.field = std::make_unique<SurfaceChargeField>(mesh, panels_.rebindNative(std::move(chain), cellMap), 0.0);
    }
    const auto& nextPanels = result.field->panels_.panels();
    require(nextPanels.size() == oldPanels.size() && panels_.closed() == result.field->panels_.closed()
        && panels_.beginEndpointKind() == result.field->panels_.beginEndpointKind()
        && panels_.endEndpointKind() == result.field->panels_.endEndpointKind(),
        "mesh rebind changed material topology or endpoint identity");
    const std::size_t count = oldPanels.size();
    std::vector<bool> used(count, false);
    result.sourceToTargetPanel.resize(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        std::size_t j=count;
        for(int candidate:result.field->panels_.nativeChain().cellToGraphs().at(cellMap[oldPanels[i].cellIndex])) {
            const auto& next=nextPanels[candidate];
            if(oldPanels[i].coordinateIsX!=next.coordinateIsX
                || !sameChart(oldPanels[i].beginFrame,next.beginFrame)
                || !sameChart(oldPanels[i].endFrame,next.endFrame))continue;
            require(j==count,"retained charged interval has several target material panels");
            j=static_cast<std::size_t>(candidate);
        }
        require(j<count,"retained charged interval has no target material panel");
        require(j < count && !used[j], "material panel map is not a bijection");
        used[j] = true;
        result.sourceToTargetPanel[i] = j;
        validatePanel(oldPanels[i], nextPanels[j]);
    }
    for (std::size_t i = 0; i < count; ++i)
        require(panels_.closed()
                ? result.sourceToTargetPanel[(i+1)%count] == (result.sourceToTargetPanel[i]+1)%count
                : result.sourceToTargetPanel[i] == i,
            "material panel map changed chain neighbors or open endpoint order");
    const auto permute = [&](const std::vector<double>& old) {
        std::vector<double> next(count);
        for (std::size_t i = 0; i < count; ++i) next[result.sourceToTargetPanel[i]] = old[i];
        // Reuse the mature integral-state validity checks, without Q/area*area.
        const auto validated = SurfaceChargeTransport::fromIntegratedCharge(
            result.field->charge_.panelMeasure(), next, panels_.closed());
        return validated.panelCharge();
    };
    result.field->assignIntegratedCharge(permute(charge_.panelCharge()));
    result.history.reserve(history.size());
    for (const auto& values : history) result.history.push_back(permute(values));
    return result;
}
}
