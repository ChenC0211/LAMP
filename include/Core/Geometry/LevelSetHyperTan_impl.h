#pragma once

using namespace std;

template<typename Mesh>
LevelSetHyperTan::LevelSetHyperTan(Mesh& mesh) : VertField<double>(mesh) {}