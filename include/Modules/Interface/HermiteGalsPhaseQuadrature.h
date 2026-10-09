#pragma once

#include "Modules/Interface/HermiteGalsPublication.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace Quantum::Interface
{
struct HermiteGalsPhaseQuadratureConfig
{
    double relativeMomentTolerance = 1e-9;
    double relativeVolumeTolerance = 1e-8;
    int maximumDepth = 24;
    std::size_t maximumSliceQueriesPerCell = 8192;
    std::size_t maximumPoints = 2000000;
};
struct HermiteGalsPhaseQuadraturePoint
{
    int cell = -1; // Native cell index, never a CCH DOF.
    Vector2D point{0,0};
    double physicalWeight = 0; // dV including 2*pi*r in RZ; strictly positive.
    bool liquid = false; // Certified phi<0 on the retained scalar.
};
struct HermiteGalsPhaseQuadratureCell
{
    int cell = -1;
    long double liquid = 0, gas = 0; // Measured from unmodified quadrature weights.
    long double relativeVolumeDifference = 0, relativeMomentEstimate = 0;
};
struct HermiteGalsPhaseQuadratureReport
{
    std::size_t sliceQueries = 0, points = 0;
    long double maxRelativeVolumeDifference = 0, maxRelativeMomentEstimate = 0;
};
// Read-only positive volume cubature of ONE publication. Subdivision and roots
// query its original native Q4; no scalar fit, phase inventory or weight repair.
class HermiteGalsPhaseQuadrature
{
public:
    static HermiteGalsPhaseQuadrature build(const QuadtreeMesh&, HermiteGalsPublication,
        const HermiteGalsPhaseQuadratureConfig& = {});
    const HermiteGalsPublication& publication() const noexcept;
    HermiteGalsIdentity geometryIdentity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const std::vector<HermiteGalsPhaseQuadraturePoint>& points() const noexcept;
    const std::vector<HermiteGalsPhaseQuadratureCell>& cells() const noexcept;
    const HermiteGalsPhaseQuadratureReport& report() const noexcept;
    void validateMesh(const QuadtreeMesh&) const;
private:
    struct Data;
    explicit HermiteGalsPhaseQuadrature(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
