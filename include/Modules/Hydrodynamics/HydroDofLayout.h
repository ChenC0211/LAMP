#pragma once

#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace LAMP::Hydrodynamics
{
class HydroTwoPhaseState;
class HydroIntegralDivergence;
namespace Detail { class HydroStateMeshValidationBatch; }
// Native HYDRO indexing, before boundary elimination or a pressure gauge.
// This immutable value contract owns no mesh/face geometry, clock or solver.
// Build belongs to lamp_hydrodynamics; geometric validation stays in Core.
class HydroDofLayout
{
public:
    enum class NormalAxis : unsigned char { First, Second }; // x/r, y/z

    static HydroDofLayout build(const QuadtreeMesh& mesh, const NativeLeafDomain& domain);

    // Copies share one immutable layout identity. A fresh build has a new
    // identity even for an equivalent mesh; no implicit field remapping.
    HydroDofLayout(const HydroDofLayout&) = default;
    HydroDofLayout& operator=(const HydroDofLayout&) = default;
    bool sameIdentity(const HydroDofLayout& other) const noexcept { return data_ == other.data_; }
    void validateMesh(const QuadtreeMesh& mesh) const { data_->domain.validateMeshIdentity(mesh); }
    const NativeLeafDomain& domain() const noexcept { return data_->domain; }
    int velocityCount() const noexcept { return static_cast<int>(data_->nativeFaceByVelocity.size()); }
    int pressureCount() const noexcept { return static_cast<int>(data_->domain.cellIndexByRow().size()); }
    std::size_t primaryScalarCount() const noexcept
    { return static_cast<std::size_t>(velocityCount()) + static_cast<std::size_t>(pressureCount()); }

    // Valid but unselected native entities map to -1; invalid indices reject.
    int velocityDof(int nativeFace) const
    { return index(data_->velocityByNativeFace, nativeFace, "native face index out of range"); }
    int nativeFace(int velocityDof) const
    { return index(data_->nativeFaceByVelocity, velocityDof, "velocity DOF out of range"); }
    int pressureDof(int nativeCell) const
    { return index(data_->domain.cellRowByIndex(), nativeCell, "native cell index out of range"); }
    int nativeCell(int pressureDof) const
    { return index(data_->domain.cellIndexByRow(), pressureDof, "pressure DOF out of range"); }
    NormalAxis normalAxis(int velocityDof) const
    { return index(data_->normalAxis, velocityDof, "velocity DOF out of range"); }

private:
    struct Data
    {
        NativeLeafDomain domain;
        std::vector<int> velocityByNativeFace, nativeFaceByVelocity;
        std::vector<NormalAxis> normalAxis;
    };
    explicit HydroDofLayout(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
    template<class T> static T index(const std::vector<T>& values, int row, const char* reason)
    {
        if (row < 0 || static_cast<std::size_t>(row) >= values.size())
            throw std::invalid_argument(std::string("[HydroDofLayout] ") + reason);
        return values[static_cast<std::size_t>(row)];
    }
    std::shared_ptr<const Data> data_;
};

// Exactly Nu+Np primary scalars, copied from caller data and exposed read-only.
// No gradients, C1 coefficients, face moments, ghost primary state or inventory.
// This is a layout-bound data value, not an accepted/candidate time-step state.
class HydroNativeValues
{
public:
    HydroNativeValues(HydroDofLayout layout, const std::vector<double>& normalVelocity,
        const std::vector<double>& pressure)
        : layout_(std::move(layout)), normalVelocity_(normalVelocity), pressure_(pressure)
    {
        if (normalVelocity_.size() != static_cast<std::size_t>(layout_.velocityCount()) ||
            pressure_.size() != static_cast<std::size_t>(layout_.pressureCount()))
            throw std::invalid_argument("[HydroNativeValues] velocity/pressure dimensions do not match layout");
        for (const auto* values : {&normalVelocity_, &pressure_})
            for (double value : *values) if (!std::isfinite(value))
                throw std::invalid_argument("[HydroNativeValues] nonfinite native value");
    }
    const HydroDofLayout& layout() const noexcept { return layout_; }
    const std::vector<double>& normalVelocity() const noexcept { return normalVelocity_; }
    const std::vector<double>& pressure() const noexcept { return pressure_; }
    void validateFor(const HydroDofLayout& expected, const QuadtreeMesh& mesh) const
    {
        if (!layout_.sameIdentity(expected))
            throw std::invalid_argument("[HydroNativeValues] different layout identity; explicit transfer required");
        expected.validateMesh(mesh);
    }
private:
    HydroDofLayout layout_;
    std::vector<double> normalVelocity_, pressure_;
};
namespace Detail
{
// Private to one State validation invocation. The ADDRESS alone is not proof:
// only a completed original full LIVE scan installs the actual immutable
// NativeMeshIdentity binding. Never retain across any user callback.
class HydroStateMeshValidationBatch
{
    friend class LAMP::Hydrodynamics::HydroTwoPhaseState;
    friend class LAMP::Hydrodynamics::HydroIntegralDivergence;
    HydroStateMeshValidationBatch()=default;
    HydroStateMeshValidationBatch(const HydroStateMeshValidationBatch&)=delete;
    HydroStateMeshValidationBatch& operator=(const HydroStateMeshValidationBatch&)=delete;
    const QuadtreeMesh* mesh_=nullptr;
    std::shared_ptr<const NativeMeshIdentity> identity_;
    void clear() noexcept {mesh_=nullptr;identity_.reset();}
    void rememberAfterValidation(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) noexcept
    {mesh_=&mesh;identity_=domain.meshIdentity();}
    bool sameBinding(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) const noexcept
    {return mesh_==&mesh&&identity_&&domain.meshIdentity()==identity_;}
    void validateLayout(const QuadtreeMesh& mesh,const HydroDofLayout& layout)
    {
        if(sameBinding(mesh,layout.domain()))return;
        layout.validateMesh(mesh);
        rememberAfterValidation(mesh,layout.domain());
    }
    void validateValues(const QuadtreeMesh& mesh,const HydroNativeValues& values,const HydroDofLayout& expected)
    {
        if(!values.layout().sameIdentity(expected)){values.validateFor(expected,mesh);return;}
        validateLayout(mesh,expected);
    }
};
}
}
