#pragma once
#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include <Eigen/SuperLUSupport>
#include <mutex>

namespace Quantum::Hydrodynamics {
// Private assembly of the same constrained streamfunction on actual leaves.
// No independently evolving velocity or pressure operator is owned here.
struct HydroCompatibleMomentSpline::NativeLayout {
    using Matrix=HydroCompatibleMomentSpace::Matrix;
    struct Cell {int index;Vector2D lower,size;};
    struct Face {InterfaceNormalFluxView::FaceChart chart;int cell;bool horizontal,boundary,wall,impermeable;};
    CartesianLeafLocator locator;
    std::shared_ptr<const NativeMeshIdentity> identity;
    std::vector<Cell> cells;
    std::vector<int> slotByCell;
    std::vector<Face> faces;
    // Immutable exterior support derived from the physical boundary owner.
    // Order: west, east, south, north; tangential intervals are sorted.
    std::array<std::vector<std::array<double,2>>,4> outletSpans;
    Vector2D lower,upper;
    bool rz,generalTrace,permeable;
    int dofs=0;
    std::vector<std::array<int,9>> cellInteriorCoordinates;
    Matrix polynomial,system,data,allConstraints,allData,momentMap,kinematicConstraints;
    std::vector<int> coordinateRows;
    Eigen::VectorXd coordinateValues;
    Eigen::SuperLU<Matrix> solver;
    mutable std::mutex mutex;
    struct MaterialCache;
    struct CoordinateAssembly;
    struct ProjectionRequest {
        const NativeLayout& source;
        const Eigen::MatrixXd& coefficients;
        Eigen::VectorXd load;
        double sourceNormSquared=0;
    };
    mutable std::shared_ptr<MaterialCache> materialCache;
    mutable std::shared_ptr<CoordinateAssembly> assemblyCache;
    NativeLayout(const QuadtreeMesh&,const NativeLeafDomain&,const HydroCompatibleMomentSpace&,const HydroCompatibleMomentBoundary*,bool prepareMomentLift=true);
    Eigen::VectorXd solve(const Eigen::VectorXd&) const;
    Eigen::MatrixXd coefficients(const Eigen::VectorXd&) const;
    Eigen::VectorXd fieldCoordinates(const Eigen::MatrixXd&) const;
    Eigen::MatrixXd restrictCoefficients(const NativeLayout&,const Eigen::MatrixXd&) const;
    Eigen::Matrix<double,6,81> cellJets(int,Vector2D) const;
    Eigen::Matrix<double,3,9> cartesianAxis(int,int,double) const;
    NativeMaterialForms materialForms(const HydroCompatibleMomentMaterial&,const Eigen::MatrixXd&,double,
        const HydroCompatibleMomentComplement* =nullptr,ProjectionRequest* =nullptr) const;
    int locate(Vector2D) const;
    void essentialJet(Eigen::Matrix<double,6,1>&,Vector2D) const;
    FirstJet firstJet(const Eigen::MatrixXd&,Vector2D) const;
    FirstJet transportJet(const Eigen::MatrixXd&,Vector2D) const;
    InterfaceNormalFluxView normalFlux(const Eigen::MatrixXd&,bool velocityProducts=false) const;
    NativeFaceProducts faceProducts(const Eigen::MatrixXd&,int) const;
    NativeFaceTrace faceTrace(int) const;
    Matrix normalBoundaryConstraints(const HydroCompatibleMomentBoundary&) const;
    Matrix localCoordinateMap(const HydroCompatibleMomentComplement*,std::vector<int>& offsets) const;
    std::shared_ptr<CoordinateAssembly> coordinateAssembly(const HydroCompatibleMomentComplement*) const;
    NativeVolumeTransport volumeTransport(const HydroCompatibleMomentMassCurrent&,const HydroCompatibleMomentComplement*) const;
    NativeMomentumTransport momentumTransport(const HydroCompatibleMomentMassCurrent&,const HydroCompatibleMomentBoundary&,
        const HydroCompatibleMomentComplement*,const Eigen::VectorXd&) const;
    Eigen::VectorXd transpose(const Eigen::MatrixXd&,const std::vector<PointForce>&) const;
    Eigen::VectorXd nativeTranspose(const std::vector<PointForce>&) const;
    Matrix pointWorkGram(const std::vector<PointForce>&) const;
};
}
