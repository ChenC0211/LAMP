#pragma once
#include <stdexcept>
#include <string>

namespace Quantum::Geometry
{
// A valid scalar can still be under-resolved by the retained graph chart.
// Callers may refine/retry; they must not publish a partial zero set.
class InterfaceGeometryResolution : public std::runtime_error
{
public:
    enum class Reason { CriticalCell, RootIncidence, GraphPairing, SmallComponent, MixedResolution };
    InterfaceGeometryResolution(Reason reason, int cell, const std::string& detail)
        : std::runtime_error("[InterfaceGeometry] refinement required: " + detail),
          reason_(reason), cell_(cell) {}
    Reason reason() const noexcept { return reason_; }
    int cellIndex() const noexcept { return cell_; }
private:
    Reason reason_;
    int cell_;
};
}
