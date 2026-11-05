#pragma once

#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Mesh
{
struct NativeLeafOverlapPiece
{
    enum class Kind { Unchanged,Refined,Coarsened };
    Kind kind=Kind::Unchanged;
    int sourceCell=-1,targetCell=-1;
    Geometry::PlanarCellCut::Rectangle rectangle;
    long double capacity=0;
};
// Geometry-only coverage of selected whole nested leaves on the same
// background lattice. Both selections must cover the same physical region;
// every active leaf is checked independently. Inactive native groups and
// rectangles are empty/zero. Owns indices/measures, neither mesh nor stocks.
// Uses two actual quadtree lookups per active leaf, not all-pairs intersection.
class NativeLeafOverlap
{
public:
    static NativeLeafOverlap build(const QuadtreeMesh&,const NativeLeafDomain&,
        const QuadtreeMesh&,const NativeLeafDomain&);
    const std::vector<NativeLeafOverlapPiece>& pieces()const noexcept{return pieces_;}
    const std::vector<std::vector<std::size_t>>& bySource()const noexcept{return bySource_;}
    const std::vector<std::vector<std::size_t>>& byTarget()const noexcept{return byTarget_;}
    const std::vector<Geometry::PlanarCellCut::Rectangle>& sourceRectangles()const noexcept{return source_;}
    const std::vector<Geometry::PlanarCellCut::Rectangle>& targetRectangles()const noexcept{return target_;}
    long double maximumRelativeCoverageResidual()const noexcept{return coverage_;}
    void validateFor(const QuadtreeMesh&,const QuadtreeMesh&)const;
private:
    std::shared_ptr<const NativeMeshIdentity> sourceIdentity_,targetIdentity_;
    std::vector<Geometry::PlanarCellCut::Rectangle> source_,target_;
    std::vector<NativeLeafOverlapPiece> pieces_;
    std::vector<std::vector<std::size_t>> bySource_,byTarget_;
    long double coverage_=0;
};
}
