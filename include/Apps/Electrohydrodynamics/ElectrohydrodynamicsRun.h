#pragma once

#include "Electrohydrodynamics/IO/ElectrohydrodynamicsInput.h"
#include "RunConfig.h"

// Owns the joint electric/fluid transaction and its accepted output views.
int runElectrohydrodynamics(const ElectrohydrodynamicsInput&,const RunConfig&);
