#pragma once

#include <string>
#include <vector>
#include <optional>
#include "Core/Base/Vector2D.h"

// Input for the ordinary nodal interface. Independent of flow material,
// momentum boundaries, main-clock and output settings. Legacy readers may
// bind named inlet profiles to the owning flow's boundary patches; numerical
// Interface state receives the resolved constraints and transport options.
struct InterfaceNodalConfig
{
    struct Inflow
    {
        std::string patch;
        // phi = offset + slope*(s-origin), in code length units.
        double offset = 0.0, slope = 0.0, origin = 0.0;
        double valueAt(double coordinate) const { return offset + slope*(coordinate-origin); }
    };
    struct ScalarBoundary
    {
        // Empty preserves the original path without physical wall contacts.
        // prescribedScalar is not a contact-angle or contact-line force law.
        std::string model;
        std::vector<Vector2D> fixedContacts;
        std::vector<Inflow> inflow;
    };

    bool advect = false;
    std::string transport = "nodal";
    struct PortPhase { std::string patch; int sign=0; };
    std::vector<PortPhase> portPhases; // Whole-step pure phase: liquid=-1, gas=+1.
    int reinitializeEvery = 0;
    ScalarBoundary boundary;
    struct NeckBreakup {
        int maximumLevel=0;
        double radiusCells=0;
    };
    // Explicit finest-grid RZ neck model: upper section detachment and the
    // uncharged source-bound four-cell interior split. Charged interior
    // material/Maxwell continuation still requires its own coupling path.
    // Consumers must enforce this finest level and a shrinking physical neck.
    std::optional<NeckBreakup> neckBreakup;
};
