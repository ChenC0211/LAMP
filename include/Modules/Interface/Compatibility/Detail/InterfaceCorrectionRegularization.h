#pragma once
#include "Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"

namespace InterfaceVolumeCorrection
{
// Module policy, using only same-source Core differential basis queries.
QuarticRegularization geometricRegularization(const QuadtreeMesh& mesh,
    const InterfaceState& state,double strength);
}
