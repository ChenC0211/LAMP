#pragma once

#include <string>
#include <optional>
#include <map>
#include "Modules/Interface/Compatibility/IO/DictionaryInterface.h"

class DictionaryGeometry;
class DictionaryMesh;
struct HydrodynamicsConfig;
struct InterfaceNodalConfig;

// Strict modern complete-domain input, shared by Hydro and EHD. Read geometry
// before mesh, preserving the supported subset and original rejection order.
// Paths are supplied by the caller; this preflight creates no physical state
// or output. On failure the supplied input records may be partially filled,
// as in the original readers; discard them rather than retrying in place.
// A consumer supplying adaptation receives the shared Interface controls for
// closed conservative Hydro; other consumers retain the fixed-root contract.
namespace HydroGeometryInput
{
// Analytic initial shape only:
// r(z)=R*(1+delta*cos(halfWaves*pi*(z-zmin)/Lz))/sqrt(1+delta^2/2).
// The two axial domain ends are physical reflection sections, not walls or
// periodic pairs. No alternative geometry owner persists after initialization.
struct AxisymmetricColumn {double radius=0,relativeAmplitude=0;int halfWaves=1;};
struct ParametricShape {
    Vector2D center,radii;
    int segments=32;
    bool suppliedCap=false,necked=false;
    double neckRatio=1;
    // Formula strings are initial data, never an alternative evolving field.
    bool analytic=false;
    std::string firstExpression,secondExpression;
    std::map<std::string,double> constants;
    bool implicit=false;
    std::string implicitExpression;
    Vector2D seed;
    double traceTolerance=0; // Omitted: 1e-9 times the shorter domain side.
};
struct ParametricLoop : ParametricShape {std::vector<ParametricShape> components;};
int load(const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface, DictionaryGeometry& geometry,
    DictionaryMesh& mesh, const std::string& geometryPath, const std::string& meshPath,
    std::optional<AxisymmetricColumn>* column=nullptr,
    std::optional<InterfaceAdaptationControls>* adaptation=nullptr,
    std::optional<ParametricLoop>* parametric=nullptr);
}
