#pragma once

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include <cstdint>
#include <memory>
#include <vector>

// Retained complete-uniform scalar and SW--NE PL diagnostic area fraction.
// No regularization width, pressure gradient, velocity row or dual measure.
// In RZ the fraction remains meridional area, not revolved phase volume.
class InterfaceGeometry
{
public:
    InterfaceGeometry(const QuadtreeMesh& mesh, const VertField<double>& scalar,
        bool axisymmetric = false, CartesianGridScope scope = CartesianGridScope::RootOnly,
        LevelSetHermiteGeometry::Reconstruction reconstruction = LevelSetHermiteGeometry::Reconstruction::Cubic);
    const LevelSetHermiteGeometry::Snapshot& snapshot() const { return *source_; }
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source() const { return source_; }
    const std::vector<double>& positiveFraction() const { return fraction_; }
    bool axisymmetric() const { return axisymmetric_; }
    void validateSource(const QuadtreeMesh& mesh, const VertField<double>& scalar) const;
private:
    const QuadtreeMesh* mesh_;
    bool axisymmetric_;
    CartesianGridScope scope_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
    std::vector<double> fraction_;
    std::vector<std::uint64_t> words_;
};
