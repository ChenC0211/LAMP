#pragma once

#include "Core/IO/TecplotWriter.h"

class QuadtreeMesh;
template<typename T> class VertField;

namespace TecplotMeshOutput
{
// Copies the existing native-indexed nodal values and actual mesh connectivity.
// No interpolation, unit conversion, geometry repair or physics calculation.
// The caller sets time, strand and physical metadata before serialization.
TecplotWriter::Zone nodalZone(const QuadtreeMesh& mesh,
    const std::vector<std::string>& names, const std::vector<VertField<double>>& fields);
}
