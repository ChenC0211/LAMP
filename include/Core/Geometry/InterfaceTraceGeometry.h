#pragma once

#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

// Ordered geometric queries on one retained scalar. No native velocity row,
// pressure gradient, dual measure, material or physical traction lives here.
// Consumers validate each source-local jet for their reconstruction policy;
// this collection does not certify a nearest-point projection algorithm.
class InterfaceTraceGeometry
{
public:
    struct Point
    {
        LevelSetHermiteGeometry::Frame frame;
        double curvatureTangentialDerivative = 0.0;
        double endpointKktResidual = 0.0;
    };
    InterfaceTraceGeometry(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        std::vector<Point> points) : source_(std::move(source)), points_(std::move(points))
    {
        if (!source_) throw std::invalid_argument("interface trace queries require a retained scalar");
    }
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& source() const noexcept { return source_; }
    const std::vector<Point>& points() const noexcept { return points_; }
private:
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source_;
    std::vector<Point> points_;
};
