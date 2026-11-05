#pragma once

#include "Core/Units/NormalizationContext.h"
#include "Modules/Plasma/PicCaseInput.h"
#include <memory>
#include <optional>
#include <string>

class UnitsContext;
struct PlasmaCollisionInputs;
namespace Solvers { struct SolversConfig; }

struct PicAppConfig
{
    PicIntegrator integrator = PicIntegrator::Explicit;
};
struct PicNormalizationInput
{
    NormalizationConfig config;
    bool hasDict = false;
    bool hasDimension = false;
    bool hasSlabThickness = false;
    bool hasAny = false;
};
// Parsed shared input only; no runtime module, main clock or physical state.
struct PicSharedInput
{
    std::shared_ptr<const Solvers::SolversConfig> solvers;
    std::shared_ptr<const UnitsContext> units;
    PicNormalizationInput normalization;
    std::string normalizationSource, primaryPicPath;
    bool writeSIOutputs = false;
};
struct PicGeometryClosure
{
    bool hasDimension = false, hasSlabThickness = false;
    int dimension = 2;
    double slabThickness_m = 1.0;
};

PicAppConfig readPicAppConfig(const std::string& runDictPath="./system/run.dict");
// Original whitelist, cache extraction and default loaders. IO creates the
// original directories before reading shared dictionaries, as before.
std::optional<PicSharedInput> readPicSharedInput();
// Called after reference scales are resolved, at the original geometry stage.
PicGeometryClosure readPicGeometryClosure(const PicSharedInput&,const UnitsContext&);
// Called after unit/normalization reports, retaining the original failure prefix.
std::shared_ptr<const PlasmaCollisionInputs> readPicRunCollisions(const PicSharedInput&);
std::optional<PicCaseInput> readPicIntegratorInput(PicIntegrator,
    const UnitsContext*,const Solvers::SolversConfig*);
