#include "Electrospray/Coupling/ElectrosprayForces.h"
#include "Modules/Electricity/FittedElectricCaseConfig.h"
#include "Core/Materials/MaterialTable.h"
#include "Electrospray/Setup/ElectrosprayElectrostatics.h"
#include "Modules/Hydrodynamics/HydroFittedCapillary.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Units/UnitsContext.h"

#include <cmath>
#include <chrono>
#include <stdexcept>
#include <utility>

ElectrosprayDielectricForces::ElectrosprayDielectricForces(const QuadtreeMesh& mesh)
    : electric(mesh), polarization(mesh) {}

ElectrosprayConductingForces::ElectrosprayConductingForces(const QuadtreeMesh& mesh)
    : capillary(mesh), maxwell(mesh) {}

ElectrosprayDielectricForces assembleElectrosprayDielectricForces(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const ElectrosprayElectricState& initial,
    const HydroFittedCapillaryState& capillary, double referenceDensitySI, const FittedElectricCaseConfig& electricConfig,
    const UnitsContext& units)
{
    constexpr double eps0 = 8.8541878128e-12;
    initial.dielectric.validateMesh(mesh);
    if (units.mode() != UnitsMode::kSI || electricConfig.interfaceDiscretization != "cellHarmonic"
        || electricConfig.initialSurfaceChargeSI != 0.0 || !initial.dielectric.report().success || !initial.dielectric.report().usedSuperLU)
        throw std::invalid_argument("spray polarization requires the actual SuperLU q=0 cellHarmonic SI electric solution");
    const double length = units.scales().L0_m, voltage = units.scales().Phi0_V;
    const double velocity = length / units.scales().T0_s;
    const double pressure = referenceDensitySI * velocity * velocity;
    const double fieldScale = voltage / length;
    const double electricToHydro = eps0 * fieldScale * fieldScale / pressure;
    ElectrosprayDielectricForces result(mesh);
    if (!capillary.source) throw std::invalid_argument("spray dielectric loads require their capillary source provenance");
    if (initial.dielectric.time()!=capillary.source->sourceTime()
        || (initial.dielectric.potential().source()
            && initial.dielectric.potential().source()!=capillary.source->interfaceGeometry()))
        throw std::invalid_argument("dielectric field and capillary load must share their actual scalar and time");
    result.electric = buildFittedDielectricPolarization(mesh, domain.geometry(), initial.dielectric.graph(), initial.dielectric.potential(),
        initial.dielectric.epsilonByCell(), initial.dielectric.epsilonBySegment(), initial.dielectric.boundaries(), electricConfig.polarizationConsistency);
    result.polarization = buildFittedScalarGradientLoad(mesh, domain, result.electric.geometry,
        result.electric.permittivityByCell, result.electric.coefficientByNative, electricToHydro);
    result.total = HydroFittedForce::combine(mesh,domain,capillary,
        result.polarization.forceX,result.polarization.forceY);
    return result;
}

ElectrosprayConductingForces buildElectrosprayConductingForces(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, const SurfaceChargeSource& charge,
    const FittedSharpDielectricResult& electric, const PhaseMaterial& liquid, const PhaseMaterial& gas,
    double surfaceTensionSI, const FittedElectricCaseConfig& electricConfig,
    double capillaryHalfWidthCells, const UnitsContext& units, double sourceTime)
{
    constexpr double eps0 = 8.8541878128e-12;
    if (!snapshot || units.mode() != UnitsMode::kSI || electricConfig.model != "leakyDielectric"
        || electricConfig.interfaceDiscretization != "sharpHermiteP3" || !electric.report.success
        || !electric.report.usedSuperLU || !electric.report.usedSharpInterface
        || !electric.report.usedSurfaceCharge)
        throw std::invalid_argument("full spray Maxwell loads require the actual charged SuperLU P3 solution and SI scales");
    const double length = units.scales().L0_m;
    const double velocity = length / units.scales().T0_s;
    const double pressure = liquid.densitySI * velocity * velocity;
    const double fieldScale = units.scales().Phi0_V / length;
    const double electricToHydro = eps0 * fieldScale * fieldScale / pressure;
    const double sigma = surfaceTensionSI / (pressure * length);
    for (double value : {length, velocity, pressure, fieldScale, electricToHydro})
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument("invalid charged spray force reference scale");
    if (!std::isfinite(sigma) || sigma < 0.0
        || (surfaceTensionSI > 0.0 && sigma == 0.0)
        || !std::isfinite(capillaryHalfWidthCells) || capillaryHalfWidthCells <= 0.0)
        throw std::invalid_argument("invalid charged spray capillary coefficient or layer width");

    ElectrosprayConductingForces result(mesh);
    const auto source = HydroForceSource::fitted(domain, sourceTime, HydroForceEvolution::TimeDependent, snapshot);
    result.electricToHydro = electricToHydro;
    const bool native = snapshot->usesNativeLeaves();
    const auto spacing = native ? panels.nativeChain().minimumActiveSpacing() : snapshot->cellSize();
    const double halfWidth = capillaryHalfWidthCells * std::max(spacing.x, spacing.y);
    const auto geometryStarted = std::chrono::steady_clock::now();
    result.geometry = native
        ? Quantum::Hydrodynamics::buildNativeFittedInterfaceLoadGeometry(mesh,domain,snapshot,panels,halfWidth,
            Quantum::Hydrodynamics::FittedTransverseGradientClosure::ActiveHermiteWallGradient)
        : Quantum::Hydrodynamics::buildFittedInterfaceLoadGeometry(mesh,domain,snapshot,panels,halfWidth,
            Quantum::Hydrodynamics::FittedTransverseGradientClosure::ActiveHermiteWallGradient);
    result.geometrySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - geometryStarted).count();
    auto& capillary = result.capillary;
    assembleFittedCapillaryForce(mesh, domain, result.geometry, sigma, capillary);
    capillary.source = source;

    DielectricConfig config;
    config.model = electricConfig.model; config.interfaceDiscretization = electricConfig.interfaceDiscretization;
    config.permittivityNegative = liquid.relativePermittivity;
    config.permittivityPositive = gas.relativePermittivity;
    // This reconstruction consumes no linear/flux tolerance: its local trace
    // and traction identity gates are fixed arithmetic checks in Electricity.
    config.axisymmetric = true;
    const auto maxwellStarted = std::chrono::steady_clock::now();
    result.electric = native
        ? reconstructNativeFittedSharpDielectricTractions(mesh,result.geometry.traceGeometry(),panels,
            config,electric.potential,charge)
        : reconstructFittedSharpDielectricTractions(mesh,result.geometry.traceGeometry(),panels,
            CartesianGridScope::UniformLeaves,config,electric.potential,charge);
    result.maxwell = Quantum::Hydrodynamics::buildFittedInterfaceTractionLoad(mesh, domain, result.geometry,
        result.electric.geometry, result.electric.tractionComponents());
    result.maxwellSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - maxwellStarted).count();
    result.total = HydroFittedForce::combine(mesh,domain,capillary,
        result.maxwell.forceX,result.maxwell.forceY,electricToHydro);
    return result;
}
