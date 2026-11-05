#pragma once

#include "Modules/Interface/HermiteGalsPublication.h"
#include "Solvers/Poisson/ElectricLocalOperator.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Electric
{
// P1: complete UniformLeaves, one connected run of regular noncoincident graphs in each cut
// native cell and one connected element per phase. No small-cut agglomeration
// or AMR is supplied. A failed geometric/conditioning prerequisite rejects.
struct ElectricPhaseGeometryConfig
{
    double relativeMatrixTolerance = 1e-9;
    int maximumRefinements = 4; // coarse plus at most this many tighter rules.
    double minimumRelativePhaseArea = 1e-6; // ordinary meridional dA / cell dA.
    std::size_t maximumPoints = 2000000;
    std::size_t maximumSliceQueriesPerCell = 8192;
    int maximumDepth = 24;
};
struct ElectricPhaseSkeletonFace
{
    enum class Kind { Cartesian, Interface };
    Kind kind = Kind::Cartesian;
    int boundarySide = -1; // W,S,E,N; -1 is a shared interior/interface face.
    bool liquid = false; // Only meaningful for Cartesian phase-subfaces.
    // Original publication's positive-coordinate native half-face partition.
    // Each retained nonempty range is a shared physical subface; it is never
    // combined by comparing rounded roots or endpoint display coordinates.
    int nativeRow = -1, half = -1, range = -1;
    int component = -1, graph = -1; // Original same-source interface graph.
    double length = 0; // Meridional ds, also actual h_F.
    // Cartesian normal is positive coordinate; interface normal is phi<0 ->
    // phi>0. Element copies change only the outward normal's sign.
    // coordinate is NORMALIZED ARCLENGTH on this full shared subface.
    std::vector<ElectricDiscretization::FaceQuadraturePoint> nodes;
};
struct ElectricPhaseElement
{
    int nativeCell = -1;
    bool liquid = false;
    ElectricDiscretization::LocalGeometry geometry;
    std::vector<int> faces; // Shared face IDs in geometry.faces order.
};
// Passive wall-clock seconds accumulated over every coarse, rejected and
// accepted quadrature refinement. Stage fields are disjoint; total includes
// all build work and remaining is its unclassified remainder. Timings do not
// participate in quadrature resolution or numerical acceptance.
struct ElectricPhaseGeometryTimings
{
    double total_s = 0;
    // Positive volume cubature, phase-node distribution, basis and P5 moments.
    double phaseQuadrature_s = 0;
    double cartesianFaces_s = 0;
    double interfaceFaces_s = 0;
    // Unit-permittivity local matrices and the constant reconstruction checks.
    double localOperators_s = 0;
    double matrixComparison_s = 0;
    double remaining_s = 0;
};
struct ElectricPhaseGeometryReport
{
    std::size_t elements = 0, cartesianSubfaces = 0, interfaceFaces = 0;
    std::size_t volumePoints = 0, surfacePoints = 0;
    int acceptedRefinement = 0;
    double minimumRelativePhaseArea = 1;
    double maxRelativeReconstructionDifference = 0;
    double maxRelativeGradientGramDifference = 0;
    double maxRelativeStiffnessDifference = 0;
    double maxRelativeCondensedDifference = 0;
    double maxRelativeOrdinaryP5MomentDifference = 0;
    double maxConstantReconstructionDefect = 0;
    double maxRelativeConstantEnergyDefect = 0;
    ElectricPhaseGeometryTimings timings;
};
class ElectricPhaseGeometry
{
public:
    static ElectricPhaseGeometry build(const QuadtreeMesh&,Interface::HermiteGalsPublication,
        const ElectricPhaseGeometryConfig& = {});
    const Interface::HermiteGalsPublication& publication() const;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    bool axisymmetric() const;
    const std::vector<ElectricPhaseElement>& elements() const;
    const std::vector<ElectricPhaseSkeletonFace>& faces() const;
    const ElectricPhaseGeometryReport& report() const;
    void validateMesh(const QuadtreeMesh&) const;
    // Identity and every clock field are compared; equal sampled scalar values
    // cannot substitute for the immutable publication owning this geometry.
    void validateSource(const QuadtreeMesh&,const Interface::HermiteGalsPublication&) const;
private:
    struct Data;
    explicit ElectricPhaseGeometry(std::shared_ptr<const Data>);
    const Data& checked() const;
    std::shared_ptr<const Data> data_;
};
}
