#pragma once

#include <iostream>
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/Vert.h"
#include "Core/Base/Base.h"

namespace Interpolation
{
    template<typename Mesh, typename Field>
    double biQuadInterpolation
    (
        const Mesh& mesh, 
        Field& field, 
        const Cell& cell, 
        const Vector2D& coord
    );

    template<typename Mesh, typename Field>
    double biLinearInterpolation
    (
        const Mesh& mesh, 
        Field& field, 
        const Cell& cell, 
        const Vector2D& coord
    );

    template<typename Mesh, typename Field>
    double averageInterpolation
    (
        const Mesh& mesh, 
        Field& field, 
        int preI
    );

    double leastSquareInterpolation
    (
        const std::vector<Vector2D>& points,
        const std::vector<double>& values,
        const Vector2D& aimCoord,
        bool secondOrder
    );

    std::vector<double> leastSquareInterpolationCoeff
    (
        const std::vector<Vector2D>& points,
        const Vector2D& aimCoord,
        bool secondOrder = false
    );

    double diracInterpolation
    (
        const Vector2D& LagrangianCoord,
        const Vector2D& EulerianCoord,
        const Vector2D& h
    );
};

#include "Core/Base/Interpolation_impl.h"
