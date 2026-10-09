#include "Modules/Hydrodynamics/HydroFittedForce.h"
#include "Modules/Hydrodynamics/HydroFittedCapillary.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <cmath>
#include <stdexcept>
#include <utility>

struct HydroFittedForce::Data
{
    Data(const QuadtreeMesh& mesh, std::shared_ptr<const HydroForceSource> source)
        : x(mesh,0.0), y(mesh,0.0), source(std::move(source)) {}
    VEdgeField<double> x;
    HEdgeField<double> y;
    std::shared_ptr<const HydroForceSource> source;
};

HydroFittedForce HydroFittedForce::combine(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const HydroFittedCapillaryState& capillary,
    const VEdgeField<double>& externalX, const HEdgeField<double>& externalY, double scale)
{
    domain.validateMeshIdentity(mesh);
    if (!capillary.source || !std::isfinite(scale) || scale <= 0.0)
        throw std::invalid_argument("fitted force combination requires capillary provenance and positive finite scale");
    capillary.source->validate(HydroDomainKind::FittedRZ, domain.geometry().meshIdentity(),
        capillary.source->sourceTime(), 0.0, capillary.source->interfaceGeometry().get());
    auto storage = std::make_shared<Data>(mesh,capillary.source);
    const auto add = [scale](const auto& base, const auto& external, auto& total)
    {
        if (!base.value || !external.value || base.size() != total.size() || external.size() != total.size())
            throw std::invalid_argument("fitted force component storage does not match the receiving mesh");
        for (int i=0; i<total.size(); ++i)
        {
            if (base.label(i) != total.label(i) || external.label(i) != total.label(i))
                throw std::invalid_argument("fitted force component has stale native labels");
            // Preserve addition for already normalized input and the original
            // separate conversion then addition for external traction units.
            const double converted = scale == 1.0 ? external.value[i] : external.value[i] * scale;
            const double sum = base.value[i] + converted;
            if (!std::isfinite(converted) || !std::isfinite(sum)
                || (external.value[i] != 0.0 && converted == 0.0))
                throw std::invalid_argument("unrepresentable fitted force normalization or sum");
            total.value[i] = sum;
        }
    };
    add(capillary.forceX, externalX, storage->x);
    add(capillary.forceY, externalY, storage->y);
    for (const auto& face : domain.faces())
    {
        if (face.velocityRow >= 0) continue;
        const double value = face.horizontal ? storage->y.value[face.nativeRow-mesh.vEdgeIndexNum()]
            : storage->x.value[face.nativeRow];
        if (value != 0.0)
            throw std::invalid_argument("fitted force reached a constrained or inactive velocity row");
    }
    HydroFittedForce result;
    result.data_ = std::move(storage); // Publish only after every check succeeds.
    return result;
}

const HydroFittedForce::Data& HydroFittedForce::data() const
{
    if (!data_) throw std::logic_error("fitted force has not been assembled");
    return *data_;
}
std::shared_ptr<const HydroForceSource> HydroFittedForce::source() const noexcept
{ return data_ ? data_->source : nullptr; }
HydroReadOnlyField<VEdgeField<double>> HydroFittedForce::x() const
{ return HydroReadOnlyField<VEdgeField<double>>(data().x); }
HydroReadOnlyField<HEdgeField<double>> HydroFittedForce::y() const
{ return HydroReadOnlyField<HEdgeField<double>>(data().y); }
HydroForceDensity HydroFittedForce::density() const
{
    const auto& stored = data();
    HydroForceDensity result;
    result.faceX = &stored.x; result.faceY = &stored.y;
    result.source = stored.source; result.evaluationTime = stored.source->sourceTime();
    return result;
}
