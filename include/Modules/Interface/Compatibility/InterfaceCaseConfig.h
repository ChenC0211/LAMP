#pragma once

#include <vector>
#include <string>

// Fitted free-surface selection and contact constraints. The App resolves
// dictionary IDs to geometry before calling Interface numerical services.
struct FittedInterfaceCaseConfig
{
    std::string fluidInterfaceLevelSet;
    int movingBoundaryId = -1;
    std::string contactLineModel;
    std::vector<int> pinnedContactPointIds;
};

// Independent prescribed-input case data, not a main clock or output policy.
// Geometry IDs are resolved once by the App; numerical services take the
// resulting domain, coordinates and frozen fields, not this input object.
struct InterfaceCaseConfig
{
    enum class Transport { PinnedToAxis, Conservative };
    Transport transport = Transport::PinnedToAxis;
    int excludedPhase = -1, pinnedPoint = -1; // 1-based geometry dictionary IDs
    double initialHeight = 0, peakSpeed = 0, requiredHalfWidthCells = 0; // code units
    std::vector<int> stationaryWalls, symmetry, singlePhasePorts;
};
