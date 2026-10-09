#pragma once

#include <iostream>
#include <array>
#include "Core/Base/Enumeration.h"
#include "Core/Mesh/Vert.h"
#include "Core/Mesh/Edge.h"
#include "Core/Mesh/Cell.h"
#include "Core/Lagrangian/InterfaceNode.h"
#include "Core/Base/Utility.h"
#include <Eigen/Dense>

namespace FiniteDifference
{
    template<typename Field, typename Mesh>
    double dx(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dy(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dxx(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dyy(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dxy(const Field& dxField, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dyx(const Field& dyField, const Mesh& mesh, const Vert& vert);

    template<typename Field, typename Mesh>
    double dx_backward(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dy_backward(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dx_forward(const Field& field, const Mesh& mesh, const Vert& vert);
    template<typename Field, typename Mesh>
    double dy_forward(const Field& field, const Mesh& mesh, const Vert& vert);

    template<typename Field1, typename Mesh, typename Field2>
    double dx_backward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dxxField);
    template<typename Field1, typename Mesh, typename Field2>
    double dy_backward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dyyField);
    template<typename Field1, typename Mesh, typename Field2>
    double dx_forward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dxxField);
    template<typename Field1, typename Mesh, typename Field2>
    double dy_forward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dyyField);

    template<typename LevelSet, typename Mesh>
    std::array<double, 6> IIMDiscrete(const LevelSet& LS, const Mesh& mesh, std::vector<InterfaceNode>& nodes, Vert& vert, InterfaceNode& node, const int relatedVertLabel[6], const int inPhaseNum, const std::vector<double>& b);
    template<typename LevelSet>
    void IIMDiscrete_RHS(const LevelSet& LS, std::vector<InterfaceNode>& nodes, Vert& vert, InterfaceNode& node, DiscCoeff& discCoeff, std::array<double, 6>& d, const double& constant = 0.0);
};

namespace FiniteDifferenceForLevelSet
{
    template<typename Field1, typename Mesh, typename Field2>
    double dx_backward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dxxField);
    template<typename Field1, typename Mesh, typename Field2>
    double dy_backward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dyyField);
    template<typename Field1, typename Mesh, typename Field2>
    double dx_forward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dxxField);
    template<typename Field1, typename Mesh, typename Field2>
    double dy_forward(const Field1& field, const Mesh& mesh, const Vert& vert, const Field2& dyyField);
}

#include "Core/Discrete/FiniteDifference_impl.h"
