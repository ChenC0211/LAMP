#pragma once

#include <memory>
#include <cstddef>

class Hydrodynamics;
class HydroFittedDomain;
class NativeMeshIdentity;
class InterfaceState;
namespace Quantum::Geometry::Parametric {template<std::size_t> struct PolynomialCurve;}
enum class HydroDomainKind;
namespace LevelSetHermiteGeometry { class Snapshot; }

enum class HydroForceEvolution { Frozen, TimeDependent };

// Immutable provenance published together with a completed force assembly.
// Shares existing mesh identity / retained interface; copies neither mesh nor
// field values and is not a second geometry-version registry. Native payloads
// remain caller-owned and must not be modified while this load is borrowed.
class HydroForceSource
{
public:
    static std::shared_ptr<const HydroForceSource> fitted(const HydroFittedDomain& domain,
        double sourceTime, HydroForceEvolution evolution,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry = {});

    double sourceTime() const noexcept { return sourceTime_; }
    HydroForceEvolution evolution() const noexcept { return evolution_; }
    const std::shared_ptr<const LevelSetHermiteGeometry::Snapshot>& interfaceGeometry() const noexcept
    { return interfaceGeometry_; }
    const std::shared_ptr<const Quantum::Geometry::Parametric::PolynomialCurve<5>>& parametricGeometry() const noexcept
    { return parametricGeometry_; }

private:
    friend class Hydrodynamics;
    friend class HydroFittedForce;
    HydroForceSource(HydroDomainKind kind, std::shared_ptr<const NativeMeshIdentity> meshIdentity,
        double sourceTime, HydroForceEvolution evolution,
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry,
        std::shared_ptr<const Quantum::Geometry::Parametric::PolynomialCurve<5>> parametricGeometry = {});
    void validate(HydroDomainKind kind,
        const std::shared_ptr<const NativeMeshIdentity>& expectedMesh,
        double stepBegin, double dt, const LevelSetHermiteGeometry::Snapshot* expectedInterface,
        const InterfaceState* expectedParametric = nullptr) const;

    HydroDomainKind kind_;
    std::shared_ptr<const NativeMeshIdentity> meshIdentity_;
    double sourceTime_;
    HydroForceEvolution evolution_;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> interfaceGeometry_;
    std::shared_ptr<const Quantum::Geometry::Parametric::PolynomialCurve<5>> parametricGeometry_;
};
