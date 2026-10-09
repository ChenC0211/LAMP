#pragma once
#include "Modules/Plasma/Boundary/MovingEBConfig.h"
class DictionaryGeometry;

// Interpret the retained legacy section. Absent input returns disabled defaults.
MovingEBConfig readMovingEBConfig(const DictionaryGeometry& geometry);
