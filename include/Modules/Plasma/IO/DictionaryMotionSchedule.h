#pragma once
#include "Modules/Plasma/Boundary/LegacyGridMotion.h"
class DictionaryGeometry;
LegacyGridMotionConfig readMotionSchedule(const DictionaryGeometry& geometry);
