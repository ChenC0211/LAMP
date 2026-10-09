#pragma once

#include <iostream>
#include <vector>
#include <string>
#include "Core/Base/Vector2D.h"
#include "Core/Field/VertField.h"
#include "Core/Base/Structures.h"

class LevelSetHyperTan : public VertField<double>
{
public:
    LevelSetHyperTan();
    template<typename Mesh>
    LevelSetHyperTan(Mesh& mesh);
    ~LevelSetHyperTan();

    inline double thickness() const { return thickness_; }
    inline void setThickness(double thickness) { thickness_ = thickness; }

private:
    double thickness_ = 0.04;

};

#include "Core/Geometry/LevelSetHyperTan_impl.h"
