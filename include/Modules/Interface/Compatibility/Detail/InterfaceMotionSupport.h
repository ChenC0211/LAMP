#pragma once

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include <array>
#include <optional>
#include <set>
#include <vector>

class InterfaceTransportDomain;
class LevelSetInterfacePanels;
class NativeLeafDomain;
class UniformLeafGrid;

// Private, single-motion workspace. All geometry/layout inputs are synchronous
// borrows and must outlive it. No scalar, mesh or accepted-state copy is owned.
// The caller has already validated the source/domain/layout correspondence.
class InterfaceMotionSupport
{
public:
    struct Band
    {
        double fullWidth, zeroWidth, portWidth;
        Band(double requiredHalfWidthCells, bool quartic, double maxh);
    };
    enum class Stage { Transported, Corrected };

    InterfaceMotionSupport(const QuadtreeMesh& mesh, const InterfaceTransportDomain& transport,
        const LevelSetHermiteGeometry::Snapshot& source, const UniformLeafGrid* grid,
        const std::vector<std::array<int, 2>>& logical, double maxh, Band band);

    // Select a certified old-scalar band once, before any characteristic.
    // Later stages only grow the required cell set and certify this SAME band.
    void prepare(const LevelSetInterfacePanels& previous);
    void checkCandidate(const LevelSetInterfacePanels& panels,
        const LevelSetHermiteGeometry::Snapshot& geometry, Stage stage);

    const Band& band() const { return band_; }
    const std::set<int>& pressureCells() const { return pressureCells_; }
    const std::set<int>& supportCells() const { return supportCells_; }

private:
    void append(const LevelSetInterfacePanels& panels);
    const QuadtreeMesh& mesh_;
    const NativeLeafDomain& domain_;
    const LevelSetHermiteGeometry::Snapshot& source_;
    const UniformLeafGrid* grid_;
    const std::vector<std::array<int, 2>>& logical_;
    const double maxh_, supportHalo_;
    Band band_;
    std::set<int> pressureCells_, supportCells_;
    std::optional<CartesianLeafLocator> locator_;
};
