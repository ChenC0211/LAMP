#pragma once

#include "Modules/Interface/GeometricVofPublication.h"
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace LAMP::Interface
{
enum class VofSmoothGeometryFailure {
    SourceKind, UnsupportedMesh, UnresolvedCell, BoundaryEndpoint,
    CornerEndpoint, UnmatchedEndpoint, InvalidTangent, InfeasibleVolume,
    InvalidCurve, VolumeResidual, ForeignSource
};
class VofSmoothGeometryError : public std::runtime_error
{
public:
    VofSmoothGeometryError(VofSmoothGeometryFailure,int nativeCell,const char*);
    VofSmoothGeometryFailure reason() const noexcept { return reason_; }
    int nativeCell() const noexcept { return nativeCell_; }
private:
    VofSmoothGeometryFailure reason_;
    int nativeCell_;
};
struct VofSmoothGeometryConfig
{
    long double relativeVolumeTolerance=1e-12L;
    // Frozen, transverse boundary intersections only; no moving contact model.
    bool allowBoundaryEndpoints=false;
};
struct VofSmoothWeightedPoint
{
    Geometry::PlanarCellCut::Point point;
    // Actual physical dA (XY) or 2*pi*r*dA (RZ), never a normalized weight.
    long double weight=0;
};
struct VofSmoothNode
{
    int nativeRow=-1;
    Geometry::PlanarCellCut::Point point,tangent; // Shared unit direction.
};
struct VofSmoothArc
{
    using Point=Geometry::PlanarCellCut::Point;
    int nativeCell=-1;
    std::array<std::size_t,2> nodes{};
    // Authoritative quartic Bezier controls, liquid on the left. End controls
    // are the shared owned nodes; convex-hull/monotonicity certificates apply.
    std::array<Point,5> controlPoints{};
    // Rounded power projection x(t)=sum power[k]*t^k for algebraic consumers.
    // sample()/derivative() evaluate the authoritative controls directly.
    std::array<Point,5> power{};
    unsigned graphAxis=0;
    long double volumeCorrection=0; // Signed graph-transverse bump, length units.
    Point sample(long double t) const;
    Point derivative(long double t) const;
    Point normal(long double t) const; // Unit liquid->gas normal.
};
struct VofSmoothCell
{
    int nativeCell=-1,arc=-1;
    Geometry::PlanarCellCut::Rectangle rectangle;
    std::array<long double,2> inventory{},measured{};
    std::array<std::vector<VofSmoothWeightedPoint>,2> phaseQuadrature;
};
struct VofSmoothGeometryReport
{
    std::size_t mixedCells=0,nodes=0,boundaryEndpoints=0;
    // Passive counts for the direct XY path; RZ keeps its declared policy.
    std::size_t xyDirectCandidates=0,xyGeometryCertifications=0,
        xyPhaseQuadratureBuilds=0,xyFitPhaseMeasurements=0,xyVolumeNewtonSteps=0;
    unsigned maximumCertificateSubdivisionDepth=0;
    long double maximumOwnPhaseRelativeResidual=0;
};
struct VofSmoothStripMeasure
{
    // phase uses the original shared-FT, minor-preserving normalization.
    // rawPhase and absoluteError describe the independent positive measures.
    std::array<long double,2> phase{},rawPhase{},absoluteError{};
    long double requestedVolume=0,relativeCoverageResidual=0;
    std::size_t polynomialRoots=0;
};
// Immutable geometry derived from actual independent VOF stocks. It supplies
// read-only phase measures and electric geometry; no query changes stocks.
// The actual LocalPlic publication remains the source owner, including clock.
class VofSmoothGeometry
{
public:
    static VofSmoothGeometry prepare(const QuadtreeMesh&,
        const GeometricVofPublication&,const VofSmoothGeometryConfig& = {});
    const GeometricVofPublication& source() const;
    GeometricVofIdentity sourceIdentity() const;
    GeometricVofTimeLayer sourceTime() const;
    const std::vector<VofSmoothCell>& cells() const;
    const std::vector<VofSmoothNode>& nodes() const;
    const std::vector<VofSmoothArc>& arcs() const;
    const VofSmoothGeometryReport& report() const;
    // Classify the actual retained graph/rectangle closure, never the source
    // PLIC plane. A point on an indistinguishable curve trace rejects.
    unsigned phaseAt(int nativeCell,Geometry::PlanarCellCut::Point) const;
    // Actual whole-face inward strip on a uniform XY cell: normalAxis=0/1,
    // upperFace selects E/N, and volume is in [0, actual cell capacity].
    // Face anchor, requested FT and tangent width remain separate through
    // exact root/phase decisions, including sub-ulp world widths.
    // No inventory or geometry changes; ambiguous/unrepresentable measures reject.
    VofSmoothStripMeasure measureNativeStrip(const QuadtreeMesh&,
        const GeometricVofPublication& actualSource,int nativeCell,unsigned normalAxis,
        bool upperFace,long double actualSweptVolume) const;
    bool sameInstance(const VofSmoothGeometry&) const noexcept;
    void validateFor(const QuadtreeMesh&,const GeometricVofPublication&) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit VofSmoothGeometry(std::shared_ptr<const Data>);
};
}
