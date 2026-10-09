#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cstddef>
#include <memory>
#include <vector>

class QuadtreeMesh;
class HydroFittedDomain;
template<typename T> class CellField;
template<typename T> class VertField;
namespace HydroFittedDetail { class BoundaryMomentum; }

// Immutable geometry-dependent routing. A native coarse/fine stencil is used
// only with its original active halo. Boundary-reaching components use shared
// physical flux segments, joined to old rows only across a regular side.
class HydroFittedMomentumLayout
{
public:
    static std::shared_ptr<const HydroFittedMomentumLayout> build(
        const QuadtreeMesh& mesh, const HydroFittedDomain& domain);
    const std::vector<unsigned char>& nativeInteriorRows() const noexcept { return nativeInterior_; }
    const HydroFittedDetail::BoundaryMomentum* boundary() const noexcept { return boundary_.get(); }
    void validate(const QuadtreeMesh& mesh, const HydroFittedDomain& domain) const;
private:
    std::vector<unsigned char> nativeInterior_;
    std::shared_ptr<const HydroFittedDetail::BoundaryMomentum> boundary_;
    std::shared_ptr<const HydroFittedDomain> domain_;
};

struct HydroFittedMomentumSystem
{
    Eigen::SparseMatrix<double> matrix;
    Eigen::VectorXd rhs;
    std::vector<double> densityByVelocityRow;
    std::size_t regularNativeRows = 0, boundaryPieceRows = 0;
    std::size_t adaptiveNativeRows = 0;
};

// Native indices throughout. Produced by the existing Hydro RZ assembly for
// selected interior rows only; no independent momentum discretization here.
struct HydroFittedNativeInterior
{
    std::shared_ptr<const HydroFittedMomentumLayout> layout;
    Eigen::SparseMatrix<double,Eigen::RowMajor> matrix;
    Eigen::VectorXd rhs;
    std::vector<double> density, viscosity;
};

// Fixed face-fitted RZ CN/AB2 momentum assembly. Adaptive interior rows must be
// supplied from the old native assembler where its halo is active. Boundary
// coarse/fine rows integrate paired physical fluxes in the same CN/AB2 system.
// Their complete equations use the old dual-volume row normalization, retaining
// the common pressure gradient and projection mass. The vectors
// use full native indexing, but solid-only entries are never read. Wall/axis
// native constraints must be exactly zero. Body fields are force per volume;
// all material, geometry, pressure, velocity and time inputs use one consistent
// unit system. No historical pressure gradient is included in this result.
//
// PressureReservoir explicitly supplies incoming material and the normal
// first-inner-native / tangential zero-shear donor model. PressureOutlet still
// rejects incoming flow in both histories. The optional dimensional velocity
// tolerance is supplied per compact free row by the driver (e.g. uncancelled
// pressure-impulse roundoff); an empty vector means exact zero tolerance.
// It never clips velocity or drops an accepted momentum flux.
//
// The caller supplies vertex material values formed from ACTIVE fluid only.
// The returned system has compact free-velocity rows and columns. Assembly
// either returns a complete finite system or throws without publishing state.
HydroFittedMomentumSystem assembleHydroFittedMomentum(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain, double dt,
    const CellField<double>& cellRho, const CellField<double>& cellMu,
    const VertField<double>& vertexRho, const VertField<double>& vertexMu,
    const Eigen::VectorXd& currentNative, const Eigen::VectorXd& previousNative,
    const CellField<double>& bodyX, const CellField<double>& bodyY,
    const std::vector<double>& outletBackflowToleranceByVelocityRow = {},
    const HydroFittedNativeInterior* nativeInterior = nullptr);
