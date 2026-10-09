#pragma once

#include "Modules/Interface/HermiteGalsPhaseRegion.h"
#include <memory>
#include <vector>

namespace Quantum::Interface {
struct HermiteGalsInventoryRegionConfig {
    double relativeGeometryInventoryTolerance=1e-10; // |V_geometry-M| / native capacity; NOT a conservation tolerance
};
struct HermiteGalsInventoryRegionCell {
    int cell=-1;
    long double capacity=0,geometricCapacity=0,geometricLiquid=0,geometricGas=0,inventory=0;
    long double liquidWeight=0,gasWeight=0,representedInventory=0;
};
struct HermiteGalsInventoryRegionReport {
    std::size_t adjustedCells=0;
    double relativeGeometryInventoryTolerance=0;
    long double inventory=0,representedInventory=0,geometryInventoryL1=0;
    long double maxRelativeCapacityDefect=0,maxRelativeGeometryDefect=0,maxRelativeRepresentationDefect=0;
};
struct HermiteGalsInventoryRegionMeasure {
    HermiteGalsPhaseRegionMeasure geometry; // original scalar integral, including the same leaf contributions
    long double liquid=0,gas=0;
    long double signedCorrection=0,l1Correction=0;
};
// Prepared bounded phase-content weights on the SAME immutable geometry.
// Does not update M, create a new geometry, or publish an accepted step.
class HermiteGalsInventoryRegionSource {
public:
    HermiteGalsInventoryRegionSource(const QuadtreeMesh&,const HermiteGalsPhaseRegionSource&,
        const HermiteGalsInventoryRegionConfig& = {});
    const HermiteGalsState& source() const;
    const HermiteGalsPhaseRegionSource& regions() const;
    const std::vector<HermiteGalsInventoryRegionCell>& cells() const;
    const HermiteGalsInventoryRegionReport& report() const;
    void validateMesh(const QuadtreeMesh&) const;
    HermiteGalsInventoryRegionMeasure measure(const QuadtreeMesh&,Vector2D lower,Vector2D upper) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    const Data& checked() const;
};
}
