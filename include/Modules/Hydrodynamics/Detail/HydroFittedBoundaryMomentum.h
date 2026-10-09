#pragma once

#include "Modules/Hydrodynamics/Detail/HydroFittedMomentumDuals.h"
#include "Modules/Hydrodynamics/Detail/HydroNativeVelocityTrace.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <array>
#include <memory>
#include <utility>
#include <vector>

class QuadtreeMesh;
class HydroFittedDomain;
struct HydroFittedMomentumSystem;
struct HydroMomentumAssembly;
template<class T> class CellField;
template<class T> class VertField;

namespace HydroFittedDetail
{
// A linear trace of existing native velocities. Derivatives and interpolation
// reuse the native first-order least-squares reconstruction; no ghost fields.
class BoundaryMomentum
{
public:
    static std::shared_ptr<const BoundaryMomentum> build(const QuadtreeMesh& mesh,
        const HydroFittedDomain& domain, std::vector<unsigned char>& nativeInterior,
        std::vector<unsigned char> boundaryRows);
    bool contains(int native) const noexcept { return selected_[native] != 0; }
    // Adds only the selected rows; the caller retains old regular/interior
    // rows. One physical common flux is emitted with opposite signs before
    // each row's old-volume/physical-volume normalization.
    void append(const QuadtreeMesh& mesh, const HydroFittedDomain& domain, double dt,
        const CellField<double>& rho, const CellField<double>& mu,
        const VertField<double>& vertexRho, const VertField<double>& vertexMu,
        const Eigen::VectorXd& current, const Eigen::VectorXd& previous,
        const CellField<double>& bodyX, const CellField<double>& bodyY,
        HydroFittedMomentumSystem& system, HydroMomentumAssembly& assembly) const;

private:
    enum class Kind { Normal, Shear, PressurePair, General, Exterior };
    struct Connection
    {
        int negative = -1, positive = -1, normal = 0;
        Kind kind = Kind::General;
        std::vector<int> parts;
        int materialCell = -1, materialVertex = -1;
        double area = 0.0;
        // du_component/dn + du_normal/dcomponent, or 2*du_normal/dn.
        MomentumTrace strain;
        MomentumTrace currentStrain; // Legacy explicit shear's actual distance
        MomentumTrace advector;
        MomentumTrace incoming;
    };
    MomentumDuals duals_;
    std::vector<unsigned char> selected_;
    std::vector<Connection> connections_;
    std::vector<MomentumTrace> inward_;
    std::vector<MomentumTrace> radialSlope_;
    std::vector<std::array<double, 2>> ownerVolumes_;
    void prepare(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
        std::vector<unsigned char>& nativeInterior);
};
}
