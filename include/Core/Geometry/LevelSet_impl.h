#pragma once

#include "Core/Geometry/LevelSet.h"

using namespace std;

template<typename Mesh, typename>
LevelSet::LevelSet(Mesh& mesh) : VertField<double>(mesh, 0.0) {}
