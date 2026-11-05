#pragma once

#include "Pic/IO/PicInput.h"
#include <memory>
#include <string>

class NormalizationContext;
struct PlasmaRunInputs;

// Runtime contexts built only from parsed typed inputs.
struct PicRunSetup
{
    std::shared_ptr<const Solvers::SolversConfig> solvers;
    std::shared_ptr<const UnitsContext> units;
    std::shared_ptr<const NormalizationContext> normalization;
    std::shared_ptr<const PlasmaCollisionInputs> plasmaCollisions;
    bool writeSIOutputs = false;
    // Passive report provenance; no file or console policy lives in Setup.
    PicNormalizationInput normalizationInput;
    bool normalizationExplicit = false;
    std::string normalizationSource, dimensionSource, slabSource;
};

struct PicUnitsSetup
{
    std::shared_ptr<const UnitsContext> units;
    std::string warning;
};
PicUnitsSetup preparePicUnits(const PicSharedInput&);
PicRunSetup preparePicRun(const PicSharedInput&,const PicUnitsSetup&,const PicGeometryClosure&);
PlasmaRunInputs bindPicRunInputs(const PicRunSetup&);
