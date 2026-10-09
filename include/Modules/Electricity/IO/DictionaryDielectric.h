#pragma once
#include "Modules/Electricity/DielectricElectrostatics.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include <optional>
namespace DictionaryInput { struct Node; }

// Legacy EHD expands its original single tolerance at this input boundary.
DielectricConfig readDielectricConfig(const std::string& path = "./system/electricity.dict");
// Formal inputs share the physical parser, with independent Conservation
// controls/defaults. The caller owns the Node's containing document.
DielectricConfig readDielectricConfig(const DictionaryInput::Node& block);
DielectricConfig readModularDielectricConfig(const std::string& path = "./system/electricity.dict");

// Fixed full-domain Electricity input, in code units. No run/output settings,
// Hydro properties or mesh/field state. Geometry and mesh retain their shared
// readers and existing capabilities; sharp half-width is an audit only.
struct DielectricDomainInput
{
    DielectricConfig electric;
    DictionaryGeometry geometry;
    DictionaryMesh mesh;
    double diffuseThickness = 0.0;
    std::optional<double> ignoredSharpHalfWidth;
};

// dt is used only for the original explicit charge-relaxation preflight.
// Paths may be outside the current directory; no cwd changes or output files.
DielectricDomainInput readDielectricDomainInput(double dt,
    const std::string& electricPath = "system/electricity.dict",
    const std::string& geometryPath = "./system/geometry.dict",
    const std::string& meshPath = "system/mesh.dict");
