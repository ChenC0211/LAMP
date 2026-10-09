#pragma once

#include "Core/Base/Base.h"
#include "Core/Discrete/FiniteDifference.h"
using namespace std;

namespace FiniteDifference
{

template<typename Field, typename Mesh>
static void getOperatorForX(const Field& field, const Mesh& mesh, const Vert& vert, double* s, double* U) 
{
    U[0] = field.value[vert.index()];

    if (vert.type() == TVertType::west)
    {
        Cell& supCell0 = mesh.cells(vert.vertToCell[0]);
        Vert& supVert0 = mesh.verts(supCell0.cellToVert[0]);
        Vert& supVert3 = mesh.verts(supCell0.cellToVert[3]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert3.index()] - (field.value[neiVert1.index()] + field.value[neiVert3.index()])) * 0.5 + U[0];
        U[2] = field.value[neiVert2.index()];

        s[0] = supVert0.coord().x - vert.coord().x;
        s[1] = neiVert2.coord().x - vert.coord().x;
    }
    else if (vert.type() == TVertType::east)
    {
        Cell& supCell2 = mesh.cells(vert.vertToCell[2]);
        Vert& supVert1 = mesh.verts(supCell2.cellToVert[1]);
        Vert& supVert2 = mesh.verts(supCell2.cellToVert[2]);
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = field.value[neiVert0.index()];
        U[2] = (field.value[supVert1.index()] + field.value[supVert2.index()] - (field.value[neiVert1.index()] + field.value[neiVert3.index()])) * 0.5 + U[0];

        s[0] = neiVert0.coord().x - vert.coord().x;
        s[1] = supVert1.coord().x - vert.coord().x;
    }
    else if (vert.boundary() == Boundary::west || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::northWest)
    {
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);
        Vert& neiVert2To2 = mesh.verts(neiVert2.neiVertices[2]);

        U[1] = field.value[neiVert2To2.index()];
        U[2] = field.value[neiVert2.index()];

        s[0] = neiVert2To2.coord().x - vert.coord().x;
        s[1] = neiVert2.coord().x - vert.coord().x;
    }
    else if (vert.boundary() == Boundary::east || vert.boundary() == Boundary::southEast || vert.boundary() == Boundary::northEast)
    {
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert0To0 = mesh.verts(neiVert0.neiVertices[0]);
        
        U[1] = field.value[neiVert0.index()];
        U[2] = field.value[neiVert0To0.index()];

        s[0] = neiVert0.coord().x - vert.coord().x;
        s[1] = neiVert0To0.coord().x - vert.coord().x;
    }
    else
    {
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);

        U[1] = field.value[neiVert0.index()];
        U[2] = field.value[neiVert2.index()];

        s[0] = neiVert0.coord().x - vert.coord().x;
        s[1] = neiVert2.coord().x - vert.coord().x;
    }
}

template<typename Field, typename Mesh>
static void getOperatorForY(const Field& field, const Mesh& mesh, const Vert& vert, double* s, double* U)
{
    U[0] = field.value[vert.index()];

    if (vert.type() == TVertType::south)
    {
        Cell& supCell1 = mesh.cells(vert.vertToCell[1]);
        Vert& supVert0 = mesh.verts(supCell1.cellToVert[0]);
        Vert& supVert1 = mesh.verts(supCell1.cellToVert[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert1.index()] - (field.value[neiVert0.index()] + field.value[neiVert2.index()])) * 0.5 + U[0];
        U[2] = field.value[neiVert3.index()];

        s[0] = supVert0.coord().y - vert.coord().y;
        s[1] = neiVert3.coord().y - vert.coord().y;
    }
    else if (vert.type() == TVertType::north)
    {
        Cell& supCell3 = mesh.cells(vert.vertToCell[3]);
        Vert& supVert2 = mesh.verts(supCell3.cellToVert[2]);
        Vert& supVert3 = mesh.verts(supCell3.cellToVert[3]);
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]); 

        U[1] = field.value[neiVert1.index()];
        U[2] = (field.value[supVert2.index()] + field.value[supVert3.index()] - (field.value[neiVert0.index()] + field.value[neiVert2.index()])) * 0.5 + U[0];

        s[0] = neiVert1.coord().y - vert.coord().y;
        s[1] = supVert2.coord().y - vert.coord().y;
    }
    else if (vert.boundary() == Boundary::south || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::southEast)
    {
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);
        Vert& neiVert3To3 = mesh.verts(neiVert3.neiVertices[3]);

        U[1] = field.value[neiVert3To3.index()];
        U[2] = field.value[neiVert3.index()];

        s[0] = neiVert3To3.coord().y - vert.coord().y;
        s[1] = neiVert3.coord().y - vert.coord().y;
    }
    else if (vert.boundary() == Boundary::north || vert.boundary() == Boundary::northWest || vert.boundary() == Boundary::northEast)
    {
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert1To1 = mesh.verts(neiVert1.neiVertices[1]);

        U[1] = field.value[neiVert1.index()];
        U[2] = field.value[neiVert1To1.index()];

        s[0] = neiVert1.coord().y - vert.coord().y;
        s[1] = neiVert1To1.coord().y - vert.coord().y;
    }
    else
    {
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = field.value[neiVert1.index()];
        U[2] = field.value[neiVert3.index()];

        s[0] = neiVert1.coord().y - vert.coord().y;
        s[1] = neiVert3.coord().y - vert.coord().y;
    }
}

}

namespace FiniteDifference
{

template<typename Field, typename Mesh>
double dx(const Field& field, const Mesh& mesh, const Vert& vert)
{
    double s[2] = { 0.0, 0.0 };
    double U[3] = { 0.0, 0.0, 0.0 };

    getOperatorForX(field, mesh, vert, s, U);

    double result = -(s[1] + s[0]) / (s[0] * s[1]) * U[0]
        + s[1] / (s[0] * (s[1] - s[0])) * U[1]
        + s[0] / (s[1] * (s[0] - s[1])) * U[2];

    return result;
}

template<typename Field, typename Mesh>
double dy(const Field& field, const Mesh& mesh, const Vert& vert)
{
    double s[2] = { 0.0, 0.0 };
    double U[3] = { 0.0, 0.0, 0.0 };

    getOperatorForY(field, mesh, vert, s, U);

    double result = -(s[1] + s[0]) / (s[0] * s[1]) * U[0]
        + s[1] / (s[0] * (s[1] - s[0])) * U[1]
        + s[0] / (s[1] * (s[0] - s[1])) * U[2];

    return result;
}

template<typename Field, typename Mesh>
double dxx(const Field& field, const Mesh& mesh, const Vert& vert)
{
    double s[2] = { 0.0, 0.0 };
    double U[3] = { 0.0, 0.0, 0.0 };

    getOperatorForX(field, mesh, vert, s, U);

    double result = 2.0 / (s[0] * s[1]) * U[0]
        - 2.0 / (s[0] * (s[1] - s[0])) * U[1]
        - 2.0 / (s[1] * (s[0] - s[1])) * U[2];

    return result;
}

template<typename Field, typename Mesh>
double dyy(const Field& field, const Mesh& mesh, const Vert& vert)
{
    double s[2] = { 0.0, 0.0 };
    double U[3] = { 0.0, 0.0, 0.0 };

    getOperatorForY(field, mesh, vert, s, U);

    double result = 2.0 / (s[0] * s[1]) * U[0]
        - 2.0 / (s[0] * (s[1] - s[0])) * U[1]
        - 2.0 / (s[1] * (s[0] - s[1])) * U[2];

    return result;
}

template<typename Field, typename Mesh>
double dxy(const Field& dxField, const Mesh& mesh, const Vert& vert)
{
    double s[2] = { 0.0, 0.0 };
    double U[3] = { 0.0, 0.0, 0.0 };

    getOperatorForY(dxField, mesh, vert, s, U);

    double result = -(s[1] + s[0]) / (s[0] * s[1]) * U[0]
        + s[1] / (s[0] * (s[1] - s[0])) * U[1]
        + s[0] / (s[1] * (s[0] - s[1])) * U[2];

    return result;
}

template<typename Field, typename Mesh>
double dyx(const Field& dyField, const Mesh& mesh, const Vert& vert)
{
    double s[2] = { 0.0, 0.0 };
    double U[3] = { 0.0, 0.0, 0.0 };

    getOperatorForX(dyField, mesh, vert, s, U);

    double result = -(s[1] + s[0]) / (s[0] * s[1]) * U[0]
        + s[1] / (s[0] * (s[1] - s[0])) * U[1]
        + s[0] / (s[1] * (s[0] - s[1])) * U[2];

    return result;
}

template<typename Field, typename Mesh>
double dx_backward(const Field &field, const Mesh &mesh, const Vert &vert)
{
    double U[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];

    if (vert.type() == TVertType::west)
    {
        Cell& supCell0 = mesh.cells(vert.vertToCell[0]);
        Vert& supVert0 = mesh.verts(supCell0.cellToVert[0]);
        Vert& supVert3 = mesh.verts(supCell0.cellToVert[3]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert3.index()]) * 0.5;
        s = supVert0.coord().x - vert.coord().x;
    }
    else if (vert.boundary() == Boundary::west || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::northWest)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);

        U[1] = field.value[neiVert0.index()];
        s = neiVert0.coord().x - vert.coord().x;
    }

    double result = (U[1] - U[0]) / s;

    return result;
}

template<typename Field, typename Mesh>
double dy_backward(const Field &field, const Mesh &mesh, const Vert &vert)
{
    double U[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];

    if (vert.type() == TVertType::south)
    {
        Cell& supCell1 = mesh.cells(vert.vertToCell[1]);
        Vert& supVert0 = mesh.verts(supCell1.cellToVert[0]);
        Vert& supVert1 = mesh.verts(supCell1.cellToVert[1]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert1.index()]) * 0.5;
        s = supVert0.coord().y - vert.coord().y;
    }
    else if (vert.boundary() == Boundary::south || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::southEast)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);

        U[1] = field.value[neiVert1.index()];
        s = neiVert1.coord().y - vert.coord().y;
    }

    double result = (U[1] - U[0]) / s;

    return result;
}

template<typename Field, typename Mesh>
double dx_forward(const Field &field, const Mesh &mesh, const Vert &vert)
{
    double U[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];

    if (vert.type() == TVertType::east)
    {
        Cell& supCell2 = mesh.cells(vert.vertToCell[2]);
        Vert& supVert1 = mesh.verts(supCell2.cellToVert[1]);
        Vert& supVert2 = mesh.verts(supCell2.cellToVert[2]);

        U[1] = (field.value[supVert1.index()] + field.value[supVert2.index()]) * 0.5;
        s = supVert1.coord().x - vert.coord().x;
    }
    else if (vert.boundary() == Boundary::east || vert.boundary() == Boundary::southEast || vert.boundary() == Boundary::northEast)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);

        U[1] = field.value[neiVert2.index()];
        s = neiVert2.coord().x - vert.coord().x;
    }

    double result = (U[1] - U[0]) / s;

    return result;
}

template<typename Field, typename Mesh>
double dy_forward(const Field &field, const Mesh &mesh, const Vert &vert)
{
    double U[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];

    if (vert.type() == TVertType::north)
    {
        Cell& supCell3 = mesh.cells(vert.vertToCell[3]);
        Vert& supVert2 = mesh.verts(supCell3.cellToVert[2]);
        Vert& supVert3 = mesh.verts(supCell3.cellToVert[3]);

        U[1] = (field.value[supVert2.index()] + field.value[supVert3.index()]) * 0.5;
        s = supVert2.coord().y - vert.coord().y;
    }
    else if (vert.boundary() == Boundary::north || vert.boundary() == Boundary::northWest || vert.boundary() == Boundary::northEast)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = field.value[neiVert3.index()];
        s = neiVert3.coord().y - vert.coord().y;
    }

    double result = (U[1] - U[0]) / s;

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dx_backward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dxxField)
{
    double U[2] = { 0.0 };
    double dxx[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dxx[0] = dxxField.value[vert.index()];

    if (vert.type() == TVertType::west)
    {
        Cell& supCell0 = mesh.cells(vert.vertToCell[0]);
        Vert& supVert0 = mesh.verts(supCell0.cellToVert[0]);
        Vert& supVert3 = mesh.verts(supCell0.cellToVert[3]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert3.index()]) * 0.5;
        s = supVert0.coord().x - vert.coord().x;
    }
    else if (vert.boundary() == Boundary::west || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::northWest)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);

        U[1] = field.value[neiVert0.index()];
        dxx[1] = dxxField.value[neiVert0.index()];
        s = neiVert0.coord().x - vert.coord().x;
    }

    double result = (U[1] - U[0]) / s + 0.5 * s * minmod(dxx[0], dxx[1]);

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dy_backward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dyyField)
{
    double U[2] = { 0.0 };
    double dyy[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dyy[0] = dyyField.value[vert.index()];

    if (vert.type() == TVertType::south)
    {
        Cell& supCell1 = mesh.cells(vert.vertToCell[1]);
        Vert& supVert0 = mesh.verts(supCell1.cellToVert[0]);
        Vert& supVert1 = mesh.verts(supCell1.cellToVert[1]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert1.index()]) * 0.5;
        s = supVert0.coord().y - vert.coord().y;
    }
    else if (vert.boundary() == Boundary::south || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::southEast)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);

        U[1] = field.value[neiVert1.index()];
        dyy[1] = dyyField.value[neiVert1.index()];
        s = neiVert1.coord().y - vert.coord().y;
    }

    double result = (U[1] - U[0]) / s + 0.5 * s * minmod(dyy[0], dyy[1]);

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dx_forward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dxxField)
{
    double U[2] = { 0.0 };
    double dxx[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dxx[0] = dxxField.value[vert.index()];

    if (vert.type() == TVertType::east)
    {
        Cell& supCell2 = mesh.cells(vert.vertToCell[2]);
        Vert& supVert1 = mesh.verts(supCell2.cellToVert[1]);
        Vert& supVert2 = mesh.verts(supCell2.cellToVert[2]);

        U[1] = (field.value[supVert1.index()] + field.value[supVert2.index()]) * 0.5;
        s = supVert1.coord().x - vert.coord().x;
    }
    else if (vert.boundary() == Boundary::east || vert.boundary() == Boundary::southEast || vert.boundary() == Boundary::northEast)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);

        U[1] = field.value[neiVert2.index()];
        dxx[1] = dxxField.value[neiVert2.index()];
        s = neiVert2.coord().x - vert.coord().x;
    }

    double result = (U[1] - U[0]) / s - 0.5 * s * minmod(dxx[0], dxx[1]);

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dy_forward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dyyField)
{
    double U[2] = { 0.0 };
    double dyy[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dyy[0] = dyyField.value[vert.index()];

    if (vert.type() == TVertType::north)
    {
        Cell& supCell3 = mesh.cells(vert.vertToCell[3]);
        Vert& supVert2 = mesh.verts(supCell3.cellToVert[2]);
        Vert& supVert3 = mesh.verts(supCell3.cellToVert[3]);

        U[1] = (field.value[supVert2.index()] + field.value[supVert3.index()]) * 0.5;
        s = supVert2.coord().y - vert.coord().y;
    }
    else if (vert.boundary() == Boundary::north || vert.boundary() == Boundary::northWest || vert.boundary() == Boundary::northEast)
    {
        return 0.0;
    }
    else
    {
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = field.value[neiVert3.index()];
        dyy[1] = dyyField.value[neiVert3.index()];
        s = neiVert3.coord().y - vert.coord().y;
    }

    double result = (U[1] - U[0]) / s - 0.5 * s * minmod(dyy[0], dyy[1]);

    return result;
}

template<typename LevelSet, typename Mesh>
std::array<double, 6> IIMDiscrete(const LevelSet &LS, const Mesh &mesh, 
    std::vector<InterfaceNode>& nodes, Vert &vert, InterfaceNode& node, 
    const int relatedVertLabel[6], DiscCoeff& discCoeff, const int inPhaseNum, 
    const std::vector<double>& b)
{
    const double LS0 = LS.value[vert.index()];
    const Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));
    
    Eigen::VectorXd _b(6);
    for (int i = 0; i < 6; i++)
    {
        _b[i] = b[i];
    }

    double k_pos = LS0 < 0.0 ? node.iimCoeff.k_vac : node.iimCoeff.k_liq;
    double k_neg = LS0 < 0.0 ? node.iimCoeff.k_liq : node.iimCoeff.k_vac;
    double jump_k = sgn(LS0) * (node.iimCoeff.k_liq - node.iimCoeff.k_vac);
    double jump_f = sgn(LS0) * (node.iimCoeff.f_liq - node.iimCoeff.f_vac);
    double f_neg = LS0 < 0.0 ? node.iimCoeff.f_liq : node.iimCoeff.f_vac;
    double beta_n = LS0 < 0.0 ? node.iimCoeff.beta_liq : node.iimCoeff.beta_vac;
    double beta_p = LS0 < 0.0 ? node.iimCoeff.beta_vac : node.iimCoeff.beta_liq;
    double beta_p_inv = 1.0 / beta_p;
    double beta_n_p_inv = beta_n * beta_p_inv;

    std::array<double, 6> d = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    Eigen::MatrixXd A(6, 6);
    double _xi[6] = { 0.0 };
    double _eta[6] = { 0.0 };
    // if (vert.neiVertices[0] == -1)
    // {
    //     if (vert.relatedVertsNum != 5)
    //     {
    //         cout << "[FiniteDifference] ERROR: unexpected number of related vertices." << endl;
    //         exit(0);
    //     }

    //     for (int col = 0; col < 5; col++)
    //     {
    //         Vert& relatedVert = mesh.verts(relatedVertLabel[col]);

    //         _xi[col] = (relatedVert.coord().x - node.coord().x) * node.normal().x + (relatedVert.coord().y - node.coord().y) * node.normal().y;
    //         _eta[col] = -(relatedVert.coord().x - node.coord().x) * node.normal().y + (relatedVert.coord().y - node.coord().y) * node.normal().x;

    //         if (col < inPhaseNum)
    //         {
    //             A(0, col) = 1.0;
    //             A(1, col) = _xi[col];
    //             A(2, col) = _eta[col];
    //             A(3, col) = 0.5 * _xi[col] * _xi[col];
    //             A(4, col) = 0.5 * _eta[col] * _eta[col];
    //             A(5, col) = _xi[col] * _eta[col];
    //         }
    //         else
    //         {
    //             if (fabs(node.coord().x) < 0.1 * cell.size().x)// 对称轴附近特殊处理
    //             {
    //                 A(0, col) = 1.0 - jump_k * 0.5 * _xi[col] * _xi[col];
    //                 A(1, col) = beta_n_p_inv * _xi[col] + (beta_n_p_inv - 1.0) * node.curvatureInPlane() * _xi[col] * _xi[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * 0.5 * _eta[col] * _eta[col];
    //                 A(2, col) = _eta[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * _xi[col] * _eta[col];
    //                 A(3, col) = 0.5 * _xi[col] * _xi[col];
    //                 A(4, col) = 0.5 * _eta[col] * _eta[col];
    //                 A(5, col) = beta_n_p_inv * _xi[col] * _eta[col];

    //             }
    //             else
    //             {
    //                 A(0, col) = 1.0 - jump_k * 0.5 * _xi[col] * _xi[col];
    //                 A(1, col) = beta_n_p_inv * _xi[col] + (beta_n_p_inv - 1.0) * (node.curvatureInPlane() - node.normal().x / node.coord().x) * 0.5 * _xi[col] * _xi[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * 0.5 * _eta[col] * _eta[col];
    //                 A(2, col) = _eta[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * _xi[col] * _eta[col];
    //                 A(3, col) = 0.5 * _xi[col] * _xi[col];
    //                 A(4, col) = 0.5 * _eta[col] * _eta[col];
    //                 A(5, col) = beta_n_p_inv * _xi[col] * _eta[col];
    //             }
    //         }
    //     }
    //     Vert& relatedVert = mesh.verts(vert.neiVertices[2]);
    //     _xi[5] = (-relatedVert.coord().x - node.coord().x) * node.normal().x + (relatedVert.coord().y - node.coord().y) * node.normal().y;
    //     _eta[5] = -(-relatedVert.coord().x - node.coord().x) * node.normal().y + (relatedVert.coord().y - node.coord().y) * node.normal().x;

    //     bool neiVert2InPhase = false;
    //     if (LS.value[relatedVert.index()] * LS.value[vert.index()] > 0.0)
    //     {
    //         A(0, 5) = 1.0;
    //         A(1, 5) = _xi[5];
    //         A(2, 5) = _eta[5];
    //         A(3, 5) = 0.5 * _xi[5] * _xi[5];
    //         A(4, 5) = 0.5 * _eta[5] * _eta[5];
    //         A(5, 5) = _xi[5] * _eta[5];
    //         neiVert2InPhase = true;
    //     }
    //     else
    //     {
    //         if (fabs(node.coord().x) < 0.1 * cell.size().x)// 对称轴附近特殊处理
    //         {
    //             A(0, 5) = 1.0 - jump_k * 0.5 * _xi[5] * _xi[5];
    //             A(1, 5) = beta_n_p_inv * _xi[5] + (beta_n_p_inv - 1.0) * node.curvatureInPlane() * _xi[5] * _xi[5] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * 0.5 * _eta[5] * _eta[5];
    //             A(2, 5) = _eta[5] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * _xi[5] * _eta[5];
    //             A(3, 5) = 0.5 * _xi[5] * _xi[5];
    //             A(4, 5) = 0.5 * _eta[5] * _eta[5];
    //             A(5, 5) = beta_n_p_inv * _xi[5] * _eta[5];
    //         }
    //         else
    //         {
    //             A(0, 5) = 1.0 - jump_k * 0.5 * _xi[5] * _xi[5];
    //             A(1, 5) = beta_n_p_inv * _xi[5] + (beta_n_p_inv - 1.0) * (node.curvatureInPlane() - node.normal().x / node.coord().x) * 0.5 * _xi[5] * _xi[5] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * 0.5 * _eta[5] * _eta[5];
    //             A(2, 5) = _eta[5] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * _xi[5] * _eta[5];
    //             A(3, 5) = 0.5 * _xi[5] * _xi[5];
    //             A(4, 5) = 0.5 * _eta[5] * _eta[5];
    //             A(5, 5) = beta_n_p_inv * _xi[5] * _eta[5];
    //         }
    //     }

    //     Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
    //     Eigen::VectorXd x = lu.solve(_b);

    //     discCoeff.discNum = 5;
    //     for (int i = 0; i < 6; i++)
    //     {
    //         discCoeff.discOpr[i] = 0.0;
    //     }
    //     for (int i = 0; i < 5; i++)
    //     {
    //         discCoeff.discOpr[i] = x[i];
    //         if (relatedVertLabel[i] == vert.neiVertices[2])
    //         {
    //             discCoeff.discOpr[i] += x[5];
    //         }
    //     }

    //     // 右端项
    //     double b2(0.0), b4(0.0), b6(0.0), b8(0.0), b10(0.0), b12(0.0);
    //     int n = neiVert2InPhase ? 5 : 6;
    //     for (int i = inPhaseNum; i < n; i++)
    //     {
    //         b2 += x[i];
    //         b4 += _xi[i] * x[i];
    //         b6 += _eta[i] * x[i];
    //         b8 += 0.5 * _xi[i] * _xi[i] * x[i];
    //         b10 += 0.5 * _eta[i] * _eta[i] * x[i];
    //         b12 += _xi[i] * _eta[i] * x[i];
    //     }

    //     if (fabs(node.coord().x) < 0.1 * cell.size().x)
    //     {
    //         d0 = b2 - k_pos * b8;
    //         d1 = b6 + node.curvatureInPlane() * b12;
    //         d2 = b10 - 2.0 * b8;
    //         d3 = beta_p_inv * (b4 + 2.0 * node.curvatureInPlane() * b8 - node.curvatureInPlane() * b10);
    //         d4 = b12 * beta_p_inv;
    //         d5 = b8 * jump_f;
    //     }
    //     else
    //     {
    //         d0 = b2 - k_pos * b8;
    //         d1 = b6 + node.normal().y / node.coord().x * b8 + node.curvatureInPlane() * b12;
    //         d2 = b10 - b8;
    //         d3 = beta_p_inv * (b4 + (node.curvatureInPlane() - node.normal().x / node.coord().x) * b8 - node.curvatureInPlane() * b10);
    //         d4 = b12 * beta_p_inv;
    //         d5 = b8 * jump_f;
    //     }

    //     double eta1 = -(nodes[node.relatedNodes[1]].coord().x - node.coord().x) * node.normal().y + (nodes[node.relatedNodes[1]].coord().y - node.coord().y) * node.normal().x;
    //     double eta2 = -(nodes[node.relatedNodes[2]].coord().x - node.coord().x) * node.normal().y + (nodes[node.relatedNodes[2]].coord().y - node.coord().y) * node.normal().x;

    //     discCoeff.rhs
    //         = (d0 + (2.0 * d2 - (eta1 + eta2) * d1) / (eta1 * eta2)) * omega[0]
    //         + ((eta2 * d1 - 2.0 * d2) / (eta1 * (eta2 - eta1))) * omega[1]
    //         + (2.0 * d2 - eta1 * d1) / (eta2 * (eta2 - eta1)) * omega[2]
    //         + (d3 - (eta1 + eta2) / (eta1 * eta2) * d4) * lambda[0]
    //         + (eta2 / (eta1 * (eta2 - eta1)) * d4) * lambda[1]
    //         + (-eta1 / (eta2 * (eta2 - eta1)) * d4) * lambda[2]
    //         + d5 - constant;
    // }
    // else
    {
        double A_inv = 1.0 / (1.0 + node.normal().x * node.normal().x);
        double B_A_inv = (1.0 + node.normal().y * node.normal().y) * A_inv;
        double C_A_inv = (2.0 * node.normal().x * node.normal().y) * A_inv;
        if (vert.relatedVertsNum != 6)
        {
            cout << "[FiniteDifference] ERROR: unexpected number of related vertices." << endl;
            exit(0);
        }

        for (int col = 0; col < 6; col++)
        {
            Vert& relatedVert = mesh.verts(relatedVertLabel[col]);

            _xi[col] = (relatedVert.coord().x - node.coord().x) * node.normal().x + (relatedVert.coord().y - node.coord().y) * node.normal().y;
            _eta[col] = -(relatedVert.coord().x - node.coord().x) * node.normal().y + (relatedVert.coord().y - node.coord().y) * node.normal().x;

            if (col < inPhaseNum)
            {
                A(0, col) = 1.0;
                A(1, col) = _xi[col];
                A(2, col) = _eta[col];
                A(3, col) = 0.5 * _xi[col] * _xi[col];
                A(4, col) = 0.5 * _eta[col] * _eta[col];
                A(5, col) = _xi[col] * _eta[col];
            }
            else
            {
                if (fabs(node.coord().x) < 0.1 * cell.size().x)// 对称轴附近特殊处理
                {
                    A(0, col) = 1.0 - A_inv * jump_k * 0.5 * _xi[col] * _xi[col];
                    A(1, col) = beta_n_p_inv * _xi[col] + B_A_inv * (beta_n_p_inv - 1.0) * node.curvatureInPlane() * 0.5 * _xi[col] * _xi[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * 0.5 * _eta[col] * _eta[col];
                    A(2, col) = _eta[col] + C_A_inv * node.curvatureInPlane() * (1.0 - beta_n_p_inv) * 0.5 * _xi[col] * _xi[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * _xi[col] * _eta[col];
                    A(3, col) = 0.5 * _xi[col] * _xi[col];
                    A(4, col) = 0.5 * _eta[col] * _eta[col];
                    A(5, col) = C_A_inv * (beta_n_p_inv - 1.0) * 0.5 * _xi[col] * _xi[col] + beta_n_p_inv * _xi[col] * _eta[col];
                }
                else
                {
                    A(0, col) = 1.0 - jump_k * 0.5 * _xi[col] * _xi[col];
                    A(1, col) = beta_n_p_inv * _xi[col] + (beta_n_p_inv - 1.0) * (node.curvatureInPlane() - node.normal().x / node.coord().x) * 0.5 * _xi[col] * _xi[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * 0.5 * _eta[col] * _eta[col];
                    A(2, col) = _eta[col] + (1.0 - beta_n_p_inv) * node.curvatureInPlane() * _xi[col] * _eta[col];
                    A(3, col) = 0.5 * _xi[col] * _xi[col];
                    A(4, col) = 0.5 * _eta[col] * _eta[col];
                    A(5, col) = beta_n_p_inv * _xi[col] * _eta[col];
                }
            }
        }

        Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
        Eigen::VectorXd x = lu.solve(_b);

        discCoeff.discNum = 6;
        for (int i = 0; i < 6; i++)
        {
            discCoeff.discOpr[i] = x[i];
        }

        // 右端项
        double b2(0.0), b4(0.0), b6(0.0), b8(0.0), b10(0.0), b12(0.0);
        for (int i = inPhaseNum; i < 6; i++)
        {
            b2 += x[i];
            b4 += _xi[i] * x[i];
            b6 += _eta[i] * x[i];
            b8 += 0.5 * _xi[i] * _xi[i] * x[i];
            b10 += 0.5 * _eta[i] * _eta[i] * x[i];
            b12 += _xi[i] * _eta[i] * x[i];
        }

        if (fabs(node.coord().x) < 0.1 * cell.size().x)
        {
            d[0] = b2 - A_inv * k_pos * b8;
            d[1] = b6 + (C_A_inv * b8 + b12) * node.curvatureInPlane();
            d[2] = b10 - B_A_inv * b8;
            d[3] = beta_p_inv * (b4 + (B_A_inv * b8 - b10) * node.curvatureInPlane());
            d[4] = (C_A_inv * b8 + b12) * beta_p_inv;
            d[5] = A_inv * b8 * jump_f;
        }
        else
        {
            d[0] = b2 - k_pos * b8;
            d[1] = b6 + node.normal().y / node.coord().x * b8 + node.curvatureInPlane() * b12;
            d[2] = b10 - b8;
            d[3] = 1.0 * beta_p_inv * (b4 + (node.curvatureInPlane() - node.normal().x / node.coord().x) * b8 - node.curvatureInPlane() * b10);
            d[4] = b12 * beta_p_inv;
            d[5] = b8 * jump_f;
        }
    }

    return d;
}

template<typename LevelSet>
void IIMDiscrete_RHS(const LevelSet& LS, vector<InterfaceNode>& nodes, Vert& vert, InterfaceNode& node, DiscCoeff& discCoeff, std::array<double, 6>& d, const double& constant)
{
    const double LS0 = LS.value[vert.index()];

    double omega[3] = { 0.0, 0.0, 0.0 };
    double lambda[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++)
    {
        InterfaceNode& relatedNode = nodes[node.relatedNodes[i]];
        omega[i] = sgn(-LS0) * relatedNode.iimCoeff.omega;
        lambda[i] = sgn(-LS0) * relatedNode.iimCoeff.lambda;
    }

    double eta1 = -(nodes[node.relatedNodes[1]].coord().x - node.coord().x) * node.normal().y + (nodes[node.relatedNodes[1]].coord().y - node.coord().y) * node.normal().x;
    double eta2 = -(nodes[node.relatedNodes[2]].coord().x - node.coord().x) * node.normal().y + (nodes[node.relatedNodes[2]].coord().y - node.coord().y) * node.normal().x;

    discCoeff.rhs
        = (d[0] + (2.0 * d[2] - (eta1 + eta2) * d[1]) / (eta1 * eta2)) * omega[0]
        + ((eta2 * d[1] - 2.0 * d[2]) / (eta1 * (eta2 - eta1))) * omega[1]
        + (2.0 * d[2] - eta1 * d[1]) / (eta2 * (eta2 - eta1)) * omega[2]
        + (d[3] - (eta1 + eta2) / (eta1 * eta2) * d[4]) * lambda[0]
        + (eta2 / (eta1 * (eta2 - eta1)) * d[4]) * lambda[1]
        + (-eta1 / (eta2 * (eta2 - eta1)) * d[4]) * lambda[2]
        + d[5] - constant;
}

}

namespace FiniteDifferenceForLevelSet
{

template<typename Field1, typename Mesh, typename Field2>
double dx_backward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dxxField)
{
    double U[2] = { 0.0 };
    double dxx[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dxx[0] = dxxField.value[vert.index()];

    if (vert.type() == TVertType::west)
    {
        Cell& supCell0 = mesh.cells(vert.vertToCell[0]);
        Vert& supVert0 = mesh.verts(supCell0.cellToVert[0]);
        Vert& supVert3 = mesh.verts(supCell0.cellToVert[3]);
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert3.index()] - (field.value[neiVert1.index()] + field.value[neiVert3.index()])) * 0.5;
        s = supVert0.coord().x - vert.coord().x;

        return U[1] / s;
    }
    else if (vert.boundary() == Boundary::west || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::northWest)
    {
        return 0.0;
    }
    else
    {
        Edge& hEdge = mesh.hEdges(vert.vertToHEdge[0]);
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        dxx[1] = dxxField.value[neiVert0.index()];
        s = -hEdge.length();

        if (hEdge.cut())
        {
            U[1] = 0.0;
            s += hEdge.cutLength();
        }
        else
        { 
            U[1] = field.value[neiVert0.index()];    
        }
    }

    double result = (U[1] - U[0]) / s - 0.5 * s * minmod(dxx[0], dxx[1]);

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dy_backward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dyyField)
{
    double U[2] = { 0.0 };
    double dyy[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dyy[0] = dyyField.value[vert.index()];

    if (vert.type() == TVertType::south)
    {
        Cell& supCell1 = mesh.cells(vert.vertToCell[1]);
        Vert& supVert0 = mesh.verts(supCell1.cellToVert[0]);
        Vert& supVert1 = mesh.verts(supCell1.cellToVert[1]);
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);

        U[1] = (field.value[supVert0.index()] + field.value[supVert1.index()] - (field.value[neiVert0.index()] + field.value[neiVert2.index()])) * 0.5;
        s = supVert0.coord().y - vert.coord().y;

        return U[1] / s;
    }
    else if (vert.boundary() == Boundary::south || vert.boundary() == Boundary::southWest || vert.boundary() == Boundary::southEast)
    {
        return 0.0;
    }
    else
    {
        Edge& vEdge = mesh.vEdges(vert.vertToVEdge[0]);
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        dyy[1] = dyyField.value[neiVert1.index()];
        s = -vEdge.length();

        if (vEdge.cut())
        {
            U[1] = 0.0;
            s += vEdge.cutLength();
        }
        else
        {
            U[1] = field.value[neiVert1.index()]; 
        }
    }

    double result = (U[1] - U[0]) / s - 0.5 * s * minmod(dyy[0], dyy[1]);

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dx_forward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dxxField)
{
    double U[2] = { 0.0 };
    double dxx[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dxx[0] = dxxField.value[vert.index()];

    if (vert.type() == TVertType::east)
    {
        Cell& supCell2 = mesh.cells(vert.vertToCell[2]);
        Vert& supVert1 = mesh.verts(supCell2.cellToVert[1]);
        Vert& supVert2 = mesh.verts(supCell2.cellToVert[2]);
        Vert& neiVert1 = mesh.verts(vert.neiVertices[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);

        U[1] = (field.value[supVert1.index()] + field.value[supVert2.index()] - (field.value[neiVert1.index()] + field.value[neiVert3.index()])) * 0.5;
        s = supVert1.coord().x - vert.coord().x;

        return U[1] / s;
    }
    else if (vert.boundary() == Boundary::east || vert.boundary() == Boundary::southEast || vert.boundary() == Boundary::northEast)
    {
        return 0.0;
    }
    else
    {
        Edge& hEdge = mesh.hEdges(vert.vertToHEdge[1]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);
        dxx[1] = dxxField.value[neiVert2.index()];

        if (hEdge.cut())
        {
            U[1] = 0.0;
            s = hEdge.cutLength();
        }
        else
        {
            U[1] = field.value[neiVert2.index()];
            s = hEdge.length();
        }
    }

    double result = (U[1] - U[0]) / s - 0.5 * s * minmod(dxx[0], dxx[1]);

    return result;
}

template<typename Field1, typename Mesh, typename Field2>
double dy_forward(const Field1 &field, const Mesh &mesh, const Vert &vert, const Field2& dyyField)
{
    double U[2] = { 0.0 };
    double dyy[2] = { 0.0 };
    double s = 0.0;

    U[0] = field.value[vert.index()];
    dyy[0] = dyyField.value[vert.index()];

    if (vert.type() == TVertType::north)
    {
        Cell& supCell3 = mesh.cells(vert.vertToCell[3]);
        Vert& supVert2 = mesh.verts(supCell3.cellToVert[2]);
        Vert& supVert3 = mesh.verts(supCell3.cellToVert[3]);
        Vert& neiVert0 = mesh.verts(vert.neiVertices[0]);
        Vert& neiVert2 = mesh.verts(vert.neiVertices[2]);

        U[1] = (field.value[supVert2.index()] + field.value[supVert3.index()] - (field.value[neiVert0.index()] + field.value[neiVert2.index()])) * 0.5;
        s = supVert2.coord().y - vert.coord().y;

        return U[1] / s;
    }
    else if (vert.boundary() == Boundary::north || vert.boundary() == Boundary::northWest || vert.boundary() == Boundary::northEast)
    {
        return 0.0;
    }
    else
    {
        Edge& vEdge = mesh.vEdges(vert.vertToVEdge[1]);
        Vert& neiVert3 = mesh.verts(vert.neiVertices[3]);
        dyy[1] = dyyField.value[neiVert3.index()];

        if (vEdge.cut())
        {
            U[1] = 0.0;
            s = vEdge.cutLength();
        }
        else
        {
            U[1] = field.value[neiVert3.index()];
            s = vEdge.length();
        }
    }

    double result = (U[1] - U[0]) / s - 0.5 * s * minmod(dyy[0], dyy[1]);

    return result;
}

}