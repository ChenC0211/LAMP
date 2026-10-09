#pragma once

#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicOverrides.h"
#include "Modules/Plasma/Boundary/MovingEBConfig.h"
#include "Modules/Plasma/Boundary/LegacyGridMotion.h"
#include "Solvers/Poisson/PoissonFVSolverCellCentered.h"

enum class PicIntegrator { Explicit, Implicit };

// Prepared once by the caller. The legacy geometry/mesh records are retained
// because their existing construction APIs consume them directly. Running PIC
// does not reread system dictionaries or infer input presence from its cwd.
struct PicCaseInput
{
    // Diagnostic provenance only; the physical driver never opens this path.
    std::string source = "prepared PIC input";
    PicConfig pic;
    ImplicitPicOverrides implicitOverrides;
    DictionaryGeometry geometry;
    DictionaryMesh mesh;
    MovingEBConfig movingBoundary;
    LegacyGridMotionConfig gridMotion;
    PoissonFVSolverCellCentered::Options poisson;
    bool geometryFilePresent = false;
    bool haveGeometry = false;
};
