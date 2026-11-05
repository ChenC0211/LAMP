#pragma once

#include "Core/Geometry/ActiveHermiteSensitivity.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Geometry {

// Immutable, geometry-neutral completion of a full uniform formal patch.
// This does not create transported state, certify Interface identity/time, or
// enable the four-piece scalar. The typed Interface owner supplies those checks.
class ActiveQuarticHermiteCompletion {
public:
    using SourceWeight = ActiveQuarticHermiteSensitivity;
    using SymmetryPlane = LevelSetInterfacePanels::SymmetryPlane;

    struct Row {
        SourceWeight::Kind kind = SourceWeight::Kind::Vertex;
        int index = -1;
        // Output identity: vertex [value,dx,dy,dxy,dxx,dyy,dxxy,dxyy,dxxyy],
        // edge [value,normalDerivative,normalSecond], center [value].
        // sourceWeights always use the original formal source components only.
        int component = 0;
        double value = 0.0;
        bool symmetryConstrained = false;
        std::vector<SourceWeight> sourceWeights;
    };

    // Requires a complete, nonhanging uniform patch, at least 2 cells per axis.
    // Edge source identities retain their input patch ordinals.
    static ActiveQuarticHermiteCompletion fromUniformPatch(
        const QuadtreeMesh& mesh,
        const ActiveQuarticHermitePatch& patch,
        const std::vector<SymmetryPlane>& symmetryPlanes = {});

    const ActiveQuarticHermitePatch& sourcePatch() const;
    const std::vector<Row>& rows() const;
    // Physical datum order is x-fast [V0,Vmid,V1,D10,D11,D20,D21]^2.
    // Cells reference unique shared vertex/edge rows, without copied closures.
    const std::array<std::size_t, 49>& cellRows(int nativeCellIndex) const;
    std::array<double, 49> physicalCellDofs(int nativeCellIndex) const;
    const std::vector<SymmetryPlane>& symmetryPlanes() const;

    // Apply the same frozen linear completion to admissible formal rates.
    // Original 25 datum rows remain selectors. Declared symmetry validates
    // their zero odd rates; only added odd derivative rows are exact zero rows.
    std::vector<double> applyRates(
        const ActiveQuarticHermitePatch& formalRates) const;
    std::vector<SourceWeight> pullbackCell(
        int nativeCellIndex,
        const std::array<long double, 49>& completedWeights) const;
    std::vector<SourceWeight> pullbackRows(
        const std::vector<long double>& completedWeights) const;

    // Verify actual mesh identity, exact formal source data/layout and the
    // explicitly declared symmetry. Revision/full clock belong to Interface.
    void validateFor(
        const QuadtreeMesh& mesh,
        const ActiveQuarticHermitePatch& patch,
        const std::vector<SymmetryPlane>& symmetryPlanes = {}) const;

private:
    struct Data;
    explicit ActiveQuarticHermiteCompletion(std::shared_ptr<const Data> data);
    const Data& data() const;
    std::shared_ptr<const Data> data_;
};

} // namespace LAMP::Geometry
