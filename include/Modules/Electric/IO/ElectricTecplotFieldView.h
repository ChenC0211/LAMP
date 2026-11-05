#pragma once

#include "Core/IO/TecplotWriter.h"
#include <array>
#include <vector>

class QuadtreeMesh;
namespace LAMP::Electric
{
class ElectricPhaseResult;
// Read-only display mesh from the retained electric phase geometry. Curved graph
// boundaries become sampled chords; no phase averaging, fitting, new solve or
// geometry publication is performed. One FE triangle zone per present phase.
// Separate element-local nodes preserve both one-sided reconstructed P3 fields.
// Current display scope requires one simple boundary loop per phase element;
// a hole or unresolved display connection is reported rather than filled.
std::vector<TecplotWriter::Zone> buildElectricTecplotFieldView(const QuadtreeMesh&,
    const ElectricPhaseResult&,const std::array<int,2>& phaseIds);
}
