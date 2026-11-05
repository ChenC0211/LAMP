#include "Pic/Setup/PicRunSetup.h"
#include "Core/Units/UnitsContext.h"
#include "Modules/Plasma/IO/DictionaryNormalization.h"
#include "Modules/Plasma/PlasmaRunInputs.h"
#include "Modules/Plasma/PlasmaCollisionInputs.h"
#include <utility>

PicUnitsSetup preparePicUnits(const PicSharedInput& input)
{
    auto units = input.units ? input.units : std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());
    std::string normWarn;
    const UnitsScales effectiveScales = resolvePlasmaReferenceScales(*units, input.normalization.config, normWarn);

    if (units->mode() == UnitsMode::kSI)
        units = std::make_shared<UnitsContext>(UnitsContext::MakeSI(effectiveScales));
    else
        units = std::make_shared<UnitsContext>(UnitsContext::MakeCodeUnits());

    return {std::move(units),std::move(normWarn)};
}

PicRunSetup preparePicRun(const PicSharedInput& input,const PicUnitsSetup& units,const PicGeometryClosure& geometry)
{
    PicRunSetup setup;
    setup.solvers = input.solvers;
    setup.units = units.units;
    setup.writeSIOutputs = input.writeSIOutputs;
    auto normLoad = input.normalization;
    // Track where normalization settings came from (for user-facing reports).
    bool normExplicit = normLoad.hasDict;
    const std::string normSourceConfig = input.normalizationSource;
    std::string normSource = normSourceConfig;
    std::string normDimSource = normLoad.hasDimension ? normSourceConfig : "";
    std::string normSlabSource = normLoad.hasSlabThickness ? normSourceConfig : "";

    // Geometry.Domain may supply dimensional closure metadata not set in Normalization.
    bool geometryProvidedClosure = false;
    if (!normLoad.hasDimension || !normLoad.hasSlabThickness)
    {
        const auto& geom = geometry;
        if (!normLoad.hasDimension && geom.hasDimension)
        {
            normLoad.config.dimension = geom.dimension;
            normDimSource = "geometry.dict/Domain";
            geometryProvidedClosure = true;
        }
        if (!normLoad.hasSlabThickness && geom.hasSlabThickness)
        {
            normLoad.config.slabThickness_m = geom.slabThickness_m;
            normSlabSource = "geometry.dict/Domain";
            geometryProvidedClosure = true;
        }
    }
    if (geometryProvidedClosure)
    {
        normExplicit = true;
        if (!normLoad.hasDict)
            normSource = "geometry.dict/Domain";
    }

    if (normDimSource.empty())
        normDimSource = "default";
    if (normLoad.config.dimension == 2)
    {
        if (normSlabSource.empty())
            normSlabSource = "default";
    }
    else if (normSlabSource.empty())
    {
        normSlabSource = "n/a";
    }

    setup.normalization = std::make_shared<NormalizationContext>(
        NormalizationContext::MakeFromUnits(setup.units, normLoad.config));

    setup.normalizationInput = std::move(normLoad);
    setup.normalizationExplicit = normExplicit;
    setup.normalizationSource = std::move(normSource);
    setup.dimensionSource = std::move(normDimSource);
    setup.slabSource = std::move(normSlabSource);
    return setup;
}

PlasmaRunInputs bindPicRunInputs(const PicRunSetup& setup)
    {
        PlasmaRunInputs inputs;
        inputs.units = setup.units;
        inputs.normalization = setup.normalization;
        if (setup.plasmaCollisions)
        {
            inputs.gas = setup.plasmaCollisions->gas;
            inputs.crossSections = setup.plasmaCollisions->crossSections;
            inputs.reactions = setup.plasmaCollisions->reactions;
            inputs.collisions = setup.plasmaCollisions->collisions;
        }
        inputs.writeSIOutputs = setup.writeSIOutputs;
        return inputs;
    }

