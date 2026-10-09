#pragma once

#include "Core/Base/Vector2D.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <optional>

enum class HydroBoundarySide { West, South, East, North };
enum class HydroBoundaryType { NormalVelocity, PressureOutlet, FreeSlip };

namespace HydroBoundaryDetail
{
// tau_nt=0 gives the missing wall tangential trace from the first actual
// interior native value. Keep both historical evaluation orders: boundary
// half-duals first form the derivative, while interior rows divide last.
inline double zeroShearTrace(double interior, double normalSpacing, double outwardSign,
    double normalDerivative)
{
    return interior - 0.5 * normalSpacing * outwardSign * normalDerivative;
}
inline double zeroShearTrace(double interior, double normalSpacing, double outwardSign,
    double normalDifference, double tangentSpacing)
{
    return interior - 0.5 * normalSpacing * outwardSign * normalDifference / tangentSpacing;
}
}

// Static patches on a rectangular domain, independent of case/material names.
// start/end are physical tangential coordinates, aligned with cell faces.
// NormalVelocity.value is OUTWARD normal velocity (negative means supply).
// PressureOutlet.value specifies full traction (-p I + 2 mu D)n = -value*n;
// it is not a simultaneous prescribed velocity or a cell pressure gauge.
// FreeSlip is impermeable, u.n=0, with natural zero tangential traction.
struct HydroBoundaryPatch
{
    std::string name;
    HydroBoundarySide side = HydroBoundarySide::South;
    HydroBoundaryType type = HydroBoundaryType::NormalVelocity;
    double start = 0.0;
    double end = 0.0;
    double value = 0.0;
    // Incoming material: -1 negative phase, +1 positive phase, 0 unspecified. Case adapters
    // bind the existing Interface portPhases declaration to this owner.
    int incomingPhase = 0;
    // Exterior momentum reservoir. NormalVelocity defaults to value*n;
    // a pressure outlet needs explicit data before it may admit backflow.
    std::optional<Vector2D> incomingVelocity;
};

inline void validateHydroBoundaryPatches(const std::vector<HydroBoundaryPatch>& patches,
    Vector2D lower, Vector2D upper, Vector2D spacing, bool axisymmetric, bool allowFreeSlip = false)
{
    if (patches.empty()) return; // Preserve the existing closed-domain contract.
    const auto require = [](bool ok, const char* reason) {
        if (!ok) throw std::invalid_argument(std::string("[HydroBoundary] ") + reason);
    };
    require(std::isfinite(lower.x) && std::isfinite(lower.y) && std::isfinite(upper.x)
        && std::isfinite(upper.y) && std::isfinite(spacing.x) && std::isfinite(spacing.y)
        && upper.x > lower.x && upper.y > lower.y && spacing.x > 0.0 && spacing.y > 0.0,
        "invalid domain or spacing");
    bool hasOutlet = false, hasPort = false;
    for (std::size_t i = 0; i < patches.size(); ++i)
    {
        const auto& patch = patches[i];
        require(patch.incomingPhase>=-1&&patch.incomingPhase<=1,"invalid incoming pure phase");
        require(!patch.incomingVelocity||(std::isfinite(patch.incomingVelocity->x)&&std::isfinite(patch.incomingVelocity->y)),"nonfinite incoming velocity");
        require(!patch.name.empty() && patch.name.find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") == std::string::npos,
            "patch names must be nonempty filename-safe identifiers");
        require(patch.side == HydroBoundarySide::West || patch.side == HydroBoundarySide::South
            || patch.side == HydroBoundarySide::East || patch.side == HydroBoundarySide::North,
            "unknown patch side");
        require(patch.type == HydroBoundaryType::NormalVelocity || patch.type == HydroBoundaryType::PressureOutlet
            || (allowFreeSlip && patch.type == HydroBoundaryType::FreeSlip),
            "unknown patch type");
        if (patch.type == HydroBoundaryType::FreeSlip)
            require(patch.value == 0 && patch.incomingPhase == 0 && !patch.incomingVelocity,
                "freeSlip requires zero normal velocity and no exterior phase/velocity data");
        require(!(axisymmetric && lower.x == 0.0 && patch.side == HydroBoundarySide::West),
            "the radial symmetry axis has its own regularity conditions and cannot be a boundary patch");
        const bool horizontal = patch.side == HydroBoundarySide::South || patch.side == HydroBoundarySide::North;
        const double lo = horizontal ? lower.x : lower.y;
        const double hi = horizontal ? upper.x : upper.y;
        const double h = horizontal ? spacing.x : spacing.y;
        require(std::isfinite(patch.start) && std::isfinite(patch.end) && std::isfinite(patch.value)
            && patch.start >= lo && patch.end <= hi && patch.end > patch.start,
            "patch interval/value is nonfinite, empty or outside its side");
        for (double coordinate : {patch.start, patch.end})
        {
            const double position = (coordinate - lo) / h;
            require(std::isfinite(position) && std::abs(position - std::round(position))
                <= 128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(position)),
                "patch endpoints must align with the fixed uniform mesh");
        }
        for (std::size_t j = 0; j < i; ++j)
        {
            require(patches[j].name != patch.name, "duplicate patch name");
            require(patches[j].side != patch.side || patch.start >= patches[j].end || patch.end <= patches[j].start,
                "boundary patches overlap");
        }
        hasOutlet = hasOutlet || patch.type == HydroBoundaryType::PressureOutlet;
        hasPort = hasPort || patch.type != HydroBoundaryType::FreeSlip;
    }
    require(!hasPort || hasOutlet, "permeable boundary patches require at least one pressure outlet");
}
