#pragma once
#include "Modules/Interface/Compatibility/InterfaceState.h"
#include <limits>

// Read-only Fourier diagnostic on the published RZ surface. A valid report
// requires one radial graph spanning the two declared axial sections; every
// native chart must use z as its monotone coordinate. GL8 moments are numerical
// quadrature, not exact integrals or a topology-event closure.
struct InterfaceColumnMode
{
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source;
    bool singleRadiusGraph=false;
    double axialLength=0;
    double meanRadius=std::numeric_limits<double>::quiet_NaN();
    double cosineAmplitude=std::numeric_limits<double>::quiet_NaN();
    double relativeAmplitude=std::numeric_limits<double>::quiet_NaN();
    static InterfaceColumnMode evaluate(const InterfaceState&,double zmin,double zmax,double referenceRadius,int halfWaves=1);
};
