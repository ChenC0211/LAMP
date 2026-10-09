#pragma once

#include <cstdint>
#include <cstring>
#include "Core/Base/Base.h"
#include "Modules/Electricity/Legacy/Electricity.h"
#include "Modules/Electricity/LegacyNodalMaxwell.h"
#include "Core/Base/Enumeration.h"
#include "Modules/Hydrodynamics/LegacyNodalState.h"

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> Trip;

template<typename Mesh>
Electricity::Electricity(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    ensureLevelSetCount(static_cast<int>(iGeometry.levelSets.size()));
    phi.initialize(mesh, 0.0);
    Ex.initialize(mesh, 0.0);
    Ey.initialize(mesh, 0.0);
    spaceRhoe.initialize(mesh, 0.0);
}

template<typename Mesh>
void Electricity::initialize(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    cout << endl << "[Electricity] Initializing electric field..." << endl;

    const int levelSetCount = static_cast<int>(iGeometry.levelSets.size());
    ensureLevelSetCount(levelSetCount);
    for (int i = 0; i < levelSetCount; i++)
    {
        const int nodeCount = static_cast<int>(sharpInterface.geometry().nodeIndexToLabel[i].size());
        rhoe[i].initialize(nodeCount, 0.0);
        iPhi[i].initialize(nodeCount, 0.0);
        EnVac[i].initialize(nodeCount, 0.0);
        EnLiq[i].initialize(nodeCount, 0.0);
        Et[i].initialize(nodeCount, 0.0);
        accCharge[i].initialize(nodeCount, 0.0);
        expansivity[i].initialize(nodeCount, 0.0);
        condCurrent[i].initialize(nodeCount, 0.0);
        convCurrent[i].initialize(nodeCount, 0.0);
        evapCurrent[i].initialize(nodeCount, 0.0);
        dtRhoe[i].initialize(nodeCount, 0.0);
        preRhoe[i].initialize(nodeCount, 0.0);
        prePreRhoe[i].initialize(nodeCount, 0.0);
        tempRhoe[i].initialize(nodeCount, 0.0);
        Vt[i].initialize(nodeCount, 0.0);
        DtVt[i].initialize(nodeCount, 0.0);
        preAccCharge[i].initialize(nodeCount, 0.0);
        preInterfaceRhoe[i].initialize(nodeCount, 0.0);
        prePreInterfaceRhoe[i].initialize(nodeCount, 0.0);
        preNormal[i].initialize(nodeCount, Vector2D(0.0, 0.0));
        prePreNormal[i].initialize(nodeCount, Vector2D(0.0, 0.0));
        nStressMaxwell[i].initialize(nodeCount, 0.0);
        tStressMaxwell[i].initialize(nodeCount, 0.0);
    }
    spaceRhoe.initialize(mesh, 0.0);

    buildEleMatrixRZ(mesh, sharpInterface, iGeometry);
    buildEleRHSRZ(mesh, sharpInterface, iGeometry);

    SparseMatrix<double> phiMatrix(mesh.vertIndexNum(), mesh.vertIndexNum());
    phiMatrix.setFromTriplets(phiTrip_.begin(), phiTrip_.end());

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(phiMatrix);

    VectorXd x_phi = solver.solve(b_phi_);

    phi.initialize(mesh, 0.0);
    Ex.initialize(mesh, 0.0);
    Ey.initialize(mesh, 0.0);
    res.initialize(mesh, 0.0);
// #pragma omp parallel for num_threads(thread_num)
    double _maxRes = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(max: _maxRes)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        phi.value[vertI] = x_phi[vertI];

        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);
        double _r = vert.coord().x;
        double _z = vert.coord().y;

        double _r2 = _r * _r;
        double _u1 = cos(_r2) * cos(_z);   // u^-
        double _u2 = cos(_r2) * sin(_z);   // u^+
        if (iGeometry.levelSets[1].value[vertI] < 0.0)
        {
            res.value[vertI] = _u1 - phi.value[vertI];
            _maxRes = max(_maxRes, fabs(res.value[vertI]));
        }
        else
        {
            res.value[vertI] = _u2 - phi.value[vertI];
            _maxRes = max(_maxRes, fabs(res.value[vertI]));
        }
    }
    cout << "[Electricity] Potential solve completed." << endl;
    cout << "[Electricity] Max residual: " << _maxRes << endl;

    computeInterfaceCoeff(mesh, sharpInterface, iGeometry);
    computeInterfaceCoeff_RHS(mesh, sharpInterface, iGeometry);
    computeElectricField(mesh, sharpInterface, iGeometry);

    cout << "[Electricity] Electric field initialized." << endl;
}

template<typename Mesh>
void Electricity::solve(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry, vector<ParticleType>& particleTypes, VertField<double>& Vx, VertField<double>& Vy, VertField<double>& preVx, VertField<double>& preVy, const double& dt, const double& preDt, const int _loop)
{
    cout << "[Electricity] Solving electric field..." << endl;
    ensureLevelSetCount(static_cast<int>(iGeometry.levelSets.size()));
    for (int i = 0; i < iGeometry.levelSets.size(); i++)
    {
        preAccCharge[i] = accCharge[i];

        prePreInterfaceRhoe[i] = preInterfaceRhoe[i];
        preInterfaceRhoe[i] = rhoe[i];

        accCharge[i].initialize(sharpInterface.geometry().nodeIndexToLabel[i].size(), 0.0);
        
        rhoe[i].initialize(sharpInterface.geometry().nodeIndexToLabel[i].size(), 0.0);
        nextRhoe[i].initialize(sharpInterface.geometry().nodeIndexToLabel[i].size(), 0.0);
    }

    vector<ParticleType> tempParticleTypes;
    for (int i = 0; i < particleTypes.size(); i++)
    {
        tempParticleTypes.push_back(ParticleType(particleTypes[i].name(), particleTypes[i].mass(), particleTypes[i].charge()));
    }

    computeInterfaceVar(mesh, sharpInterface, iGeometry, 1, Vx, Vy, preVx, preVy, dt, preDt, _loop);

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
    {
        const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[1][i];
        InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, nodeLabel);

        rhoe[1].value[i] = preRhoe[1].value[i];
    }

    // 不参与迭代的基本粒子背景空间电荷密度
    computeSpaceChargeDensity(mesh, particleTypes, iGeometry);

    buildEleMatrixRZ(mesh, sharpInterface, iGeometry);

    SparseMatrix<double> phiMatrix(mesh.vertIndexNum(), mesh.vertIndexNum());
    phiMatrix.setFromTriplets(phiTrip_.begin(), phiTrip_.end());

    SparseLU<SparseMatrix<double>> solver;
    solver.compute(phiMatrix);

    computeInterfaceCoeff(mesh, sharpInterface, iGeometry);

    tempRhoe[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    iPhi[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    EnVac[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    EnLiq[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    Et[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    phi.initialize(mesh, 0.0);
    Ex.initialize(mesh, 0.0);
    Ey.initialize(mesh, 0.0);
    
    cout << "[Electricity] Computing surface charge residual..." << endl;
    double rhoeResidual = 0.0;
    int _iter = 0;
    do 
    {
        _iter++;

#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
        {
            const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[1][i];
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, nodeLabel);

            tempRhoe[1].value[i] = rhoe[1].value[i];
        }

        buildEleRHSRZ(mesh, sharpInterface, iGeometry);

        VectorXd x_phi = solver.solve(b_phi_);
#pragma omp parallel for num_threads(thread_num)
        for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
        {
            phi.value[vertI] = x_phi[vertI];
        }

        computeInterfaceCoeff_RHS(mesh, sharpInterface, iGeometry);
        computeElectricField(mesh, sharpInterface, iGeometry);
        computeInterfaceRhoe(sharpInterface, 1, dt, preDt, _loop);

#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
        {
            const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[1][i];
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, nodeLabel);

            rhoe[1].value[i] = (1.0 - rhoeRelaxCoeff_) * tempRhoe[1].value[i] + rhoeRelaxCoeff_ * nextRhoe[1].value[i];
        }

        rhoeResidual = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(max:rhoeResidual)
        for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
        {
            double tempResidual = fabs(rhoe[1].value[i] - tempRhoe[1].value[i]);
            if (tempResidual > rhoeResidual)
            {
                rhoeResidual = tempResidual;
            }
        }

        if (_iter % 50 == 0)
        {
            cout << "[Electricity] Iter=" << _iter << ", surface charge residual=" << rhoeResidual << endl;
        }
    } while (rhoeResidual > 1e-6 && _iter < 500);
    cout << "[Electricity] Iter=" << _iter << ", surface charge residual=" << rhoeResidual << endl;

    condCurrent[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    convCurrent[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    evapCurrent[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    dtRhoe[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    nStressMaxwell[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
    tStressMaxwell[1].initialize(sharpInterface.geometry().nodeIndexToLabel[1].size(), 0.0);
#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
    {
        const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[1][i];
        InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, nodeLabel);
        Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));

        condCurrent[1].value[i] = nonDim_epsilon / number_lambda * EnLiq[1].value[i];

        double _tempDtRhoe =  sharpInterface.interfaceInterpolation(node.coord(), rhoe[1], 1, true, false, Time::current, node.normal());
        if (node.coord().x < 0.1 * cell.size().x)
        {
            convCurrent[1].value[i] = rhoe[1].value[i] * DtVt[1].value[i] + Vt[1].value[i] * _tempDtRhoe;
        }
        else
        {
            convCurrent[1].value[i] = rhoe[1].value[i] * DtVt[1].value[i] + Vt[1].value[i] * _tempDtRhoe + Vt[1].value[i] * rhoe[1].value[i] * (-node.normal().y) / node.coord().x;
        }

        if (EnVac[1].value[i] > 0.0)
        {
            evapCurrent[1].value[i] = rhoe[1].value[i] * chi * exp(CHI_cation * (xi_cation * sqrt(fabs(EnVac[1].value[i])) - 1.0));
        }
        else
        {
            evapCurrent[1].value[i] = rhoe[1].value[i] * chi * exp(CHI_anion * (xi_anion * sqrt(fabs(EnVac[1].value[i])) - 1.0));
        }

        dtRhoe[1].value[i] = condCurrent[1].value[i] - convCurrent[1].value[i] - evapCurrent[1].value[i];

        nStressMaxwell[1].value[i] = number_Bond * 0.5 * ((EnVac[1].value[i] * EnVac[1].value[i] - Et[1].value[i] * Et[1].value[i]) - nonDim_epsilon * (EnLiq[1].value[i] * EnLiq[1].value[i]  - Et[1].value[i] * Et[1].value[i]));
        tStressMaxwell[1].value[i] = number_Bond * Et[1].value[i] * rhoe[1].value[i];
    }

    cout << "[Electricity] Surface current solved." << endl;
    
    LegacyNodalMaxwell::add(sharpInterface.geometry().nodeIndexToLabel[1], sharpInterface.geometry().allNodes(1),
        nStressMaxwell[1], tStressMaxwell[1], sharpInterface.Fxs[1], sharpInterface.Fys[1]);

    // 发射粒子及粒子运动
    for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
    {
        const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[1][i];
        InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, nodeLabel);

        for (int j = 0; j < tempParticleTypes.size(); j++)
        {
            if (tempParticleTypes[j].charge() * EnVac[1].value[i] > 0.0)
            {
                sharpInterface.setEvaporationVelocity(1, nodeLabel, emitParticle(mesh, node, tempParticleTypes[j], 1, dt));
            }
        }
    }

    cout << "[Electricity] Particle emission completed." << endl;

#pragma omp parallel for num_threads(thread_num)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);
        if (vert.onCutEdge[1])
        {
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, vert.nearestNodeLabel[1]);

            vector<double> _bnn(6, 0.0);
            _bnn[3] = 1.0;
            array<double, 6> _dnn = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.nnCoeff, vert.inPhaseNum, _bnn);
            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.nnCoeff, _dnn);

            vector<double> _btt(6, 0.0);
            _btt[4] = 1.0;
            array<double, 6> _dtt = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.ttCoeff, vert.inPhaseNum, _btt);
            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.ttCoeff, _dtt);
        }
    }

    for(int i = 0; i < tempParticleTypes.size(); i++)
    {
        moveParticle(mesh, sharpInterface, iGeometry, 1, tempParticleTypes[i], dt, true);
    }

    for(int i = 0; i < particleTypes.size(); i++)
    {
        moveParticle(mesh, sharpInterface, iGeometry, 1, particleTypes[i], dt, false);
    }

    cout << "[Electricity] Particle push completed." << endl;

    // 粒子统计
    for (int i = 0; i < particleTypes.size(); i++)
    {
        particleTypes[i].number += tempParticleTypes[i].number;
        particleTypes[i].particle.insert(particleTypes[i].particle.end(), make_move_iterator(tempParticleTypes[i].particle.begin()), make_move_iterator(tempParticleTypes[i].particle.end()));

        cout << "[Electricity] Emitted species: " << tempParticleTypes[i].name() << ", count: " << tempParticleTypes[i].number << endl;
        cout << "[Electricity] Total species: " << particleTypes[i].name() << ", total count: " << particleTypes[i].number << endl;
    }

    cout << "[Electricity] Electric field solve completed." << endl;
}

template<typename Mesh>
void Electricity::computeInterfaceVar(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, VertField<double>& Vx, VertField<double>& Vy, VertField<double>& preVx, VertField<double>& preVy, const double& dt, const double& preDt, const int _loop)
{
    preRhoe[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), 0.0);
    prePreRhoe[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), 0.0);
    expansivity[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), 0.0);
    Vt[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), 0.0);
    DtVt[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), 0.0);
    preNormal[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), Vector2D(0.0, 0.0));
    prePreNormal[levelSetI].initialize(sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(), Vector2D(0.0, 0.0));
    if (_loop == 1)
    {
        return;
    }

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(); i++)
    {
        const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[levelSetI][i];
        InterfaceNode& node = sharpInterface.geometry().interfaceNode(levelSetI, nodeLabel);

        // pre
        Vector2D preCoord = iGeometry.getPreCoord(mesh, node.coord(), Vx, Vy, preVx, preVy, dt);
        const Cell& preCell = mesh.cells(mesh.findCellLabel(preCoord));

        double _LSPre = Interpolation::biQuadInterpolation(mesh, iGeometry.levelSetsPre[levelSetI], preCell, preCoord);
        double _LSPre_dx = Interpolation::biLinearInterpolation(mesh, iGeometry.levelSetsPre_dx[levelSetI], preCell, preCoord);
        double _LSPre_dy = Interpolation::biLinearInterpolation(mesh, iGeometry.levelSetsPre_dy[levelSetI], preCell, preCoord);

        double mod_inv = mod(_LSPre_dx, _LSPre_dy);
        preNormal[levelSetI].value[i] = Vector2D(_LSPre_dx * mod_inv, _LSPre_dy * mod_inv);
        Vector2D _preNodeCoord = preCoord - preNormal[levelSetI].value[i] * _LSPre;
        if (_preNodeCoord.x < mesh.minCoord.x)
        {
            _preNodeCoord.x = 2.0 * mesh.minCoord.x - _preNodeCoord.x;
        }
        node.setPreCoord(_preNodeCoord);

        // prePre
        Vector2D prePreCoord = iGeometry.getPrePreCoord(mesh, node.coord(), Vx, Vy, preVx, preVy, dt, preDt);
        const Cell& prePreCell = mesh.cells(mesh.findCellLabel(prePreCoord));

        double _LSPrePre = Interpolation::biQuadInterpolation(mesh, iGeometry.levelSetsPrePre[levelSetI], prePreCell, prePreCoord);
        double _LSPrePre_dx = Interpolation::biLinearInterpolation(mesh, iGeometry.levelSetsPrePre_dx[levelSetI], prePreCell, prePreCoord);
        double _LSPrePre_dy = Interpolation::biLinearInterpolation(mesh, iGeometry.levelSetsPrePre_dy[levelSetI], prePreCell, prePreCoord);

        mod_inv = mod(_LSPrePre_dx, _LSPrePre_dy);
        prePreNormal[levelSetI].value[i] = Vector2D(_LSPrePre_dx * mod_inv, _LSPrePre_dy * mod_inv);
        Vector2D prePreNodeCoord = prePreCoord - prePreNormal[levelSetI].value[i] * _LSPrePre;
        if (prePreNodeCoord.x < mesh.minCoord.x)
        {
            prePreNodeCoord.x = 2.0 * mesh.minCoord.x - prePreNodeCoord.x;
        }
        node.setPrePreCoord(prePreNodeCoord);

        // 表面电荷密度
        preRhoe[levelSetI].value[i] = sharpInterface.interfaceInterpolation(node.preCoord(), preInterfaceRhoe[levelSetI], levelSetI, false, false, Time::previous, preNormal[levelSetI].value[i]);
        prePreRhoe[levelSetI].value[i] = sharpInterface.interfaceInterpolation(node.prePreCoord(), prePreInterfaceRhoe[levelSetI], levelSetI, false, false, Time::prePrevious, prePreNormal[levelSetI].value[i]);

        // 电荷积累
        accCharge[levelSetI].value[i] = sharpInterface.interfaceInterpolation(node.preCoord(), preAccCharge[levelSetI], levelSetI, false, false, Time::previous, preNormal[levelSetI].value[i]);

        // 膨胀率
        const Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));
        Vector2D extCoord1(node.coord().x - cell.size().x * node.normal().x, node.coord().y - cell.size().y * node.normal().y);
        Vector2D extCoord2(node.coord().x - 2.0 * cell.size().x * node.normal().x, node.coord().y - 2.0 * cell.size().y * node.normal().y);
        if (extCoord1.x < mesh.minCoord.x)
        {
            extCoord1.x = 2.0 * mesh.minCoord.x - extCoord1.x;
        }
        if (extCoord2.x < mesh.minCoord.x)
        {
            extCoord2.x = 2.0 * mesh.minCoord.x - extCoord2.x;
        }

        Vector2D exVel1 = Vector2D(Interpolation::biQuadInterpolation(mesh, Vx, cell, extCoord1), Interpolation::biQuadInterpolation(mesh, Vy, cell, extCoord1));
        Vector2D exVel2 = Vector2D(Interpolation::biQuadInterpolation(mesh, Vx, cell, extCoord2), Interpolation::biQuadInterpolation(mesh, Vy, cell, extCoord2));

        double extVn0 = sharpInterface.Uxs[levelSetI].value[node.index()] * node.normal().x + sharpInterface.Uys[levelSetI].value[node.index()] * node.normal().y;
        double extVn1 = exVel1.dotProduct(node.normal());
        double extVn2 = exVel2.dotProduct(node.normal());

        expansivity[levelSetI].value[i] = (3.0 * extVn0 - 4.0 * extVn1 + extVn2) / (2.0 * cell.size().mod());

        // 切向速度
        Vt[levelSetI].value[i] = -sharpInterface.Uxs[levelSetI].value[node.index()] * node.normal().y + sharpInterface.Uys[levelSetI].value[node.index()] * node.normal().x;
        DtVt[levelSetI].value[i] = sharpInterface.interfaceInterpolation(node.coord(), sharpInterface.Uxs[levelSetI], levelSetI, true, true, Time::current, node.normal());
    }
    cout << "[Electricity] Interface variables computed." << endl;
}

template<typename Mesh>
void Electricity::buildEleMatrixRZ(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    phiTrip_.clear();

    vector<vector<Eigen::Triplet<double>>> _phiCoeff(mesh.vertIndexNum());
    // for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); levelSetI++)
    // {
        for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
        {
            const int nodeLabel = sharpInterface.geometry().nodeIndexToLabel[1][i];
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, nodeLabel);

            double _r = node.coord().x;
            double _z = node.coord().y;

            // 解析函数
            double _r2 = _r * _r;
            double _u1 = cos(_r2) * cos(_z);   // u^-
            double _u2 = cos(_r2) * sin(_z);   // u^+

            // 一阶导数
            double _u1_dr = -2.0 * _r * sin(_r2) * cos(_z);
            double _u1_dz = -cos(_r2) * sin(_z);
            double _u2_dr = -2.0 * _r * sin(_r2) * sin(_z);
            double _u2_dz = cos(_r2) * cos(_z);

            // 二阶导数
            double _u1_drr = -2.0 * sin(_r2) * cos(_z) - 4.0 * _r * _r * cos(_r2) * cos(_z);
            double _u1_dzz = -cos(_r2) * cos(_z);
            double _u2_drr = -2.0 * sin(_r2) * sin(_z) - 4.0 * _r * _r * cos(_r2) * sin(_z);
            double _u2_dzz = -cos(_r2) * sin(_z);

            double _lapl_u1 = 0.0;
            double _lapl_u2 = 0.0;
            if (node.coord().x < 0.1 * mesh.cells(mesh.findCellLabel(node.coord())).size().x)
            {
                _lapl_u1 = 2.0 * _u1_drr + _u1_dzz;
                _lapl_u2 = 2.0 * _u2_drr + _u2_dzz;
            }
            else
            {
                _lapl_u1 = _u1_drr + _u1_dzz + _u1_dr / _r;
                _lapl_u2 = _u2_drr + _u2_dzz + _u2_dr / _r;
            }

            double _u1_dn = _u1_dr * node.normal().x + _u1_dz * node.normal().y;
            double _u2_dn = _u2_dr * node.normal().x + _u2_dz * node.normal().y;

            node.iimCoeff.beta_liq = beta_liq_;
            node.iimCoeff.beta_vac = beta_vac_;
            node.iimCoeff.k_liq = k_;
            node.iimCoeff.k_vac = k_;
            node.iimCoeff.f_liq = _lapl_u1 + node.iimCoeff.k_liq * _u1;
            node.iimCoeff.f_vac = _lapl_u2 + node.iimCoeff.k_vac * _u2;

            node.iimCoeff.omega = _u2 - _u1;
            node.iimCoeff.lambda = _u2_dn - _u1_dn;
        }
    // }

    dlap_.initialize(mesh.vertIndexNum(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
#pragma omp parallel for num_threads(thread_num)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);

        vert.relatedVertsNum = 0;
        for (int i = 0; i < 6; i++)
        {
            vert.relatedVerts[i] = -1;
        }
        vert.inPhaseNum = 0;

        if (vert.onCutEdge[1])
        {
            // 寻找相关节点
            vector<int> _relatedVerts = findRelatedVertLabel(mesh, vert, 1);
            vert.relatedVertsNum = static_cast<int>(_relatedVerts.size());
            for (int i = 0; i < vert.relatedVertsNum; i++)
            {
                vert.relatedVerts[i] = _relatedVerts[i];
            }
            vert.inPhaseNum = sortRelatedVertLabel(mesh, iGeometry, vert, 1, vert.relatedVerts);
        }
    }

#pragma omp parallel for num_threads(thread_num)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);

        int row = vertI;

        // if (
        //     vert.boundary() == Boundary::north || vert.boundary() == Boundary::northEast || vert.boundary() == Boundary::northWest
        //     ||
        //     vert.boundary() == Boundary::south || vert.boundary() == Boundary::southEast || vert.boundary() == Boundary::southWest
        // )// 边界
        if (vert.boundary() != Boundary::non)
        {
            _phiCoeff[row].push_back(Trip(row, vert.index(), 1.0));
            continue;
        }
        
        if (vert.onCutEdge[1])// IIM离散
        {
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, vert.nearestNodeLabel[1]);
            const Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));

            vector<double> _b(6, 0.0);
            if (node.coord().x < 0.1 * cell.size().x)
            {
                _b[3] = 2.0;
                _b[4] = 1.0;
            }
            else
            {
                _b[1] = node.normal().x / node.coord().x;
                _b[2] = -node.normal().y / node.coord().x;
                _b[3] = 1.0;
                _b[4] = 1.0;
            }

            dlap_.value[vertI] = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.lapCoeff, vert.inPhaseNum, _b);

            for (int i = 0; i < vert.relatedVertsNum; i++)
            {
                Vert& relatedVert = mesh.verts(vert.relatedVerts[i]);
                _phiCoeff[row].push_back(Trip(row, relatedVert.index(), vert.lapCoeff.discOpr[i]));
            }
            _phiCoeff[row].push_back(Trip(row, vert.index(), k_));
        }
        else if (vert.type() == TVertType::non)
        {
            double _coeff[4] = { 0.0 };
            // 东西
            if (vert.neiVertices[0] != -1 && vert.neiVertices[2] != -1)
            {
                Vert& vert0 = mesh.verts(vert.neiVertices[0]);
                Vert& vert2 = mesh.verts(vert.neiVertices[2]);
                double s0 = vert.coord().x - vert0.coord().x;
                double s2 = vert2.coord().x - vert.coord().x;
                
                _coeff[0] = 2.0 / (s0 * (s0 + s2)) - s2 / (vert.coord().x * s0 * (s0 + s2));
                _coeff[2] = 2.0 / (s2 * (s0 + s2)) + s0 / (vert.coord().x * s2 * (s0 + s2));
                _phiCoeff[row].push_back(Trip(row, vert0.index(), _coeff[0]));
                _phiCoeff[row].push_back(Trip(row, vert2.index(), _coeff[2]));
            }
            else if (vert.neiVertices[2] == -1)
            {
                Vert& vert0 = mesh.verts(vert.neiVertices[0]);
                double s0 = vert.coord().x - vert0.coord().x;
                double s2 = s0;
                
                _coeff[0] = 2.0 / (s0 * (s0 + s2)) - s2 / (vert.coord().x * s0 * (s0 + s2));
                _coeff[2] = 2.0 / (s2 * (s0 + s2)) + s0 / (vert.coord().x * s2 * (s0 + s2));
                _phiCoeff[row].push_back(Trip(row, vert0.index(), _coeff[0]));
                _phiCoeff[row].push_back(Trip(row, vert0.index(), _coeff[2]));
            }
            else if (vert.neiVertices[0] == -1)// 轴对称边界
            {
                Vert& vert2 = mesh.verts(vert.neiVertices[2]);
                double s0 = vert2.coord().x - vert.coord().x;
                double s2 = s0;
                
                _coeff[0] = 4.0 / (s0 * (s0 + s2));
                _coeff[2] = 4.0 / (s2 * (s0 + s2));
                _phiCoeff[row].push_back(Trip(row, vert2.index(), _coeff[0]));
                _phiCoeff[row].push_back(Trip(row, vert2.index(), _coeff[2]));
            }

            // 南北
            if (vert.neiVertices[1] != -1 && vert.neiVertices[3] != -1)
            {
                Vert& vert1 = mesh.verts(vert.neiVertices[1]);
                Vert& vert3 = mesh.verts(vert.neiVertices[3]);
                double s1 = vert.coord().y - vert1.coord().y;
                double s3 = vert3.coord().y - vert.coord().y;
                
                _coeff[1] = 2.0 / (s1 * (s1 + s3));
                _coeff[3] = 2.0 / (s3 * (s1 + s3));
                _phiCoeff[row].push_back(Trip(row, vert1.index(), _coeff[1]));
                _phiCoeff[row].push_back(Trip(row, vert3.index(), _coeff[3]));
            }
            else if (vert.neiVertices[3] == -1)
            {
                Vert& vert1 = mesh.verts(vert.neiVertices[1]);
                double s1 = vert.coord().y - vert1.coord().y;
                double s3 = s1;
                
                _coeff[1] = 2.0 / (s1 * (s1 + s3));
                _coeff[3] = 2.0 / (s3 * (s1 + s3));
                _phiCoeff[row].push_back(Trip(row, vert1.index(), _coeff[1]));
                _phiCoeff[row].push_back(Trip(row, vert1.index(), _coeff[3]));
            }
            else if (vert.neiVertices[1] == -1)
            {
                Vert& vert3 = mesh.verts(vert.neiVertices[3]);
                double s1 = vert3.coord().y - vert.coord().y;
                double s3 = s1;
                
                _coeff[1] = 2.0 / (s1 * (s1 + s3));
                _coeff[3] = 2.0 / (s3 * (s1 + s3));
                _phiCoeff[row].push_back(Trip(row, vert3.index(), _coeff[1]));
                _phiCoeff[row].push_back(Trip(row, vert3.index(), _coeff[3]));
            }

            _phiCoeff[row].push_back(Trip(row, vert.index(), k_ - (_coeff[0] + _coeff[1] + _coeff[2] + _coeff[3])));
        }
        else
        {
            Vert* relatedVerts[6];
            if (vert.type() == TVertType::west)
            {
                relatedVerts[0] = &vert;
                relatedVerts[1] = &mesh.verts(vert.neiVertices[1]);
                relatedVerts[2] = &mesh.verts(vert.neiVertices[2]);
                relatedVerts[3] = &mesh.verts(vert.neiVertices[3]);
                relatedVerts[4] = &mesh.verts(vert.supVert[0]);
                relatedVerts[5] = &mesh.verts(vert.supVert[1]);
            }
            else if (vert.type() == TVertType::east)
            {
                relatedVerts[0] = &vert;
                relatedVerts[1] = &mesh.verts(vert.neiVertices[0]);
                relatedVerts[2] = &mesh.verts(vert.neiVertices[1]);
                relatedVerts[3] = &mesh.verts(vert.neiVertices[3]);
                relatedVerts[4] = &mesh.verts(vert.supVert[0]);
                relatedVerts[5] = &mesh.verts(vert.supVert[1]);
            }
            else if (vert.type() == TVertType::south)
            {
                relatedVerts[0] = &vert;
                relatedVerts[1] = &mesh.verts(vert.neiVertices[0]);
                relatedVerts[2] = &mesh.verts(vert.neiVertices[2]);
                relatedVerts[3] = &mesh.verts(vert.neiVertices[3]);
                relatedVerts[4] = &mesh.verts(vert.supVert[0]);
                relatedVerts[5] = &mesh.verts(vert.supVert[1]);
            }
            else if (vert.type() == TVertType::north)
            {
                relatedVerts[0] = &vert;
                relatedVerts[1] = &mesh.verts(vert.neiVertices[0]);
                relatedVerts[2] = &mesh.verts(vert.neiVertices[1]);
                relatedVerts[3] = &mesh.verts(vert.neiVertices[2]);
                relatedVerts[4] = &mesh.verts(vert.supVert[0]);
                relatedVerts[5] = &mesh.verts(vert.supVert[1]);
            }

            MatrixXd A(6, 6);
            VectorXd b(6);

            for (int i = 0; i < 6; i++)
            {
                double _dx = relatedVerts[i]->coord().x - vert.coord().x;
                double _dy = relatedVerts[i]->coord().y - vert.coord().y;
                A(0, i) = 1.0;
                A(1, i) = _dx;
                A(2, i) = _dy;
                A(3, i) = 0.5 * _dx * _dx;
                A(4, i) = 0.5 * _dy * _dy;
                A(5, i) = _dx * _dy;

                b(i) = 0.0;
            }
            b(1) = 1.0 / vert.coord().x;
            b(3) = 1.0;
            b(4) = 1.0;

            VectorXd coeff = A.colPivHouseholderQr().solve(b);
            for (int i = 0; i < 6; i++)
            {
                _phiCoeff[row].push_back(Trip(row, relatedVerts[i]->index(), coeff(i)));
            }
            _phiCoeff[row].push_back(Trip(row, vert.index(), k_));
        }
    }

    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        phiTrip_.insert(phiTrip_.end(), _phiCoeff[i].begin(), _phiCoeff[i].end());
    }
}

template<typename Mesh>
void Electricity::buildEleRHSRZ(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    b_phi_.resize(mesh.vertIndexNum());

#pragma omp parallel for num_threads(thread_num)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);

        int row = vertI;

        double _r = vert.coord().x;
        double _z = vert.coord().y;

        double _r2 = _r * _r;
        double _u1 = cos(_r2) * cos(_z);   // u^-
        double _u2 = cos(_r2) * sin(_z);   // u^+

        // if (
        //     vert.boundary() == Boundary::north || vert.boundary() == Boundary::northEast || vert.boundary() == Boundary::northWest
        //     ||
        //     vert.boundary() == Boundary::south || vert.boundary() == Boundary::southEast || vert.boundary() == Boundary::southWest
        // )// 边界
        if (vert.boundary() != Boundary::non)
        {
            if (iGeometry.levelSets[1].value[vert.index()] < 0.0)
            {
                b_phi_[row] = _u1;
            }
            else
            {
                b_phi_[row] = _u2;
            }

            continue;
        }
        
        if (vert.onCutEdge[1])// IIM离散
        {
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, vert.nearestNodeLabel[1]);

            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.lapCoeff, dlap_.value[vertI]);

            b_phi_[row] = vert.lapCoeff.rhs;
        }
        else
        {
            // vert.lapCoeff.rhs = 0.0;
            // b_phi_[row] = vert.lapCoeff.rhs;

            // 一阶导数
            double _u1_dr = -2.0 * _r * sin(_r2) * cos(_z);
            double _u1_dz = -cos(_r2) * sin(_z);
            double _u2_dr = -2.0 * _r * sin(_r2) * sin(_z);
            double _u2_dz = cos(_r2) * cos(_z);

            // 二阶导数
            double _u1_drr = -2.0 * sin(_r2) * cos(_z) - 4.0 * _r * _r * cos(_r2) * cos(_z);
            double _u1_dzz = -cos(_r2) * cos(_z);
            double _u2_drr = -2.0 * sin(_r2) * sin(_z) - 4.0 * _r * _r * cos(_r2) * sin(_z);
            double _u2_dzz = -cos(_r2) * sin(_z);

            double _lapl_u1 = 0.0;
            double _lapl_u2 = 0.0;
            if (vert.neiVertices[0] == -1 && vert.type() == TVertType::non)
            {
                _lapl_u1 = 2.0 * _u1_drr + _u1_dzz;
                _lapl_u2 = 2.0 * _u2_drr + _u2_dzz;
            }
            else
            {
                _lapl_u1 = _u1_drr + _u1_dzz + _u1_dr / _r;
                _lapl_u2 = _u2_drr + _u2_dzz + _u2_dr / _r;
            }

            if (iGeometry.levelSets[1].value[vert.index()] < 0.0)
            {
                b_phi_[row] = _lapl_u1 + k_ * _u1;
            }
            else
            {
                b_phi_[row] = _lapl_u2 + k_ * _u2;
            }
        }
    }    
}

template<typename Mesh>
void Electricity::computeInterfaceCoeff(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    dself_.initialize(mesh.vertIndexNum(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    dn_.initialize(mesh.vertIndexNum(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    dt_.initialize(mesh.vertIndexNum(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    dx_.initialize(mesh.vertIndexNum(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    dy_.initialize(mesh.vertIndexNum(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
#pragma omp parallel for num_threads(thread_num)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);
        if (vert.onCutEdge[1])
        {
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, vert.nearestNodeLabel[1]);

            vector<double> _bself(6, 0.0);
            _bself[0] = 1.0;
            dself_.value[vertI] = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.selfCoeff, vert.inPhaseNum, _bself);

            vector<double> _bn(6, 0.0);
            _bn[1] = 1.0;
            dn_.value[vertI] = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.nCoeff, vert.inPhaseNum, _bn);

            vector<double> _bt(6, 0.0);
            _bt[2] = 1.0;
            dt_.value[vertI] = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.tCoeff, vert.inPhaseNum, _bt);

            vector<double> _bx(6, 0.0);
            _bx[1] = node.normal().x;
            _bx[2] = -node.normal().y;
            dx_.value[vertI] = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.xCoeff, vert.inPhaseNum, _bx);

            vector<double> _by(6, 0.0);
            _by[1] = node.normal().y;
            _by[2] = node.normal().x;
            dy_.value[vertI] = FiniteDifference::IIMDiscrete(iGeometry.levelSets[1], mesh, sharpInterface.geometry().allNodes(1), vert, node, vert.relatedVerts, vert.yCoeff, vert.inPhaseNum, _by);
        }
    }
}

template<typename Mesh>
void Electricity::computeInterfaceCoeff_RHS(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
#pragma omp parallel for num_threads(thread_num)
    for (int vertI = 0; vertI < mesh.vertIndexNum(); vertI++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[vertI]);
        if (vert.onCutEdge[1])
        {
            InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, vert.nearestNodeLabel[1]);

            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.selfCoeff, dself_.value[vertI]);

            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.nCoeff, dn_.value[vertI]);

            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.tCoeff, dt_.value[vertI]);

            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.xCoeff, dx_.value[vertI]);

            FiniteDifference::IIMDiscrete_RHS(iGeometry.levelSets[1], sharpInterface.geometry().allNodes(1), vert, node, vert.yCoeff, dy_.value[vertI]);
        }
    }
}

template<typename Mesh>
void Electricity::computeElectricField(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry)
{
    double _maxResDn = 0.0;
#pragma omp parallel for num_threads(thread_num) reduction(max: _maxResDn)
	for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[1].size(); i++)
	{
		InterfaceNode& node = sharpInterface.geometry().interfaceNode(1, sharpInterface.geometry().nodeIndexToLabel[1][i]);
		Vert& vert = mesh.verts(node.relatedNegVertLabel());

		// 界面电势求解
        iPhi[1].value[i] = 0.0;
		for (int j = 0; j < vert.relatedVertsNum; j++)
		{
			iPhi[1].value[i] += vert.selfCoeff.discOpr[j] * phi.value[mesh.verts(vert.relatedVerts[j]).index()];
		}
		iPhi[1].value[i] -= vert.selfCoeff.rhs;

		// 法向场强
		//		液相
        EnLiq[1].value[i] = 0.0;
		for (int j = 0; j < vert.relatedVertsNum; j++)
		{
			EnLiq[1].value[i] += vert.nCoeff.discOpr[j] * phi.value[mesh.verts(vert.relatedVerts[j]).index()];
		}
		EnLiq[1].value[i] -= vert.nCoeff.rhs;

        double _r = node.coord().x;
        double _z = node.coord().y;
        double _r2 = _r * _r;
        double _u1 = cos(_r2) * cos(_z);   // u^-
        double _u2 = cos(_r2) * sin(_z);   // u^+

        double _u1_dr = -2.0 * _r * sin(_r2) * cos(_z);
        double _u1_dz = -cos(_r2) * sin(_z);
        double _u2_dr = -2.0 * _r * sin(_r2) * sin(_z);
        double _u2_dz = cos(_r2) * cos(_z);

        double _u1_dn = _u1_dr * node.normal().x + _u1_dz * node.normal().y;
        double _u2_dn = _u2_dr * node.normal().x + _u2_dz * node.normal().y;

        _maxResDn = max(_maxResDn, fabs(_u1_dn - EnLiq[1].value[i]));

		// 		真空相
		EnVac[1].value[i] = -rhoe[1].value[i] + nonDim_epsilon * EnLiq[1].value[i];

		// 		取反
		EnLiq[1].value[i] = -EnLiq[1].value[i];
		EnVac[1].value[i] = -EnVac[1].value[i];

		// 切向场强
        Et[1].value[i] = 0.0;
		for (int j = 0; j < vert.relatedVertsNum; j++)
		{
			Et[1].value[i] += vert.tCoeff.discOpr[j] * phi.value[mesh.verts(vert.relatedVerts[j]).index()];
		}
		Et[1].value[i] -= vert.tCoeff.rhs;

		// 		取反
		Et[1].value[i] = -Et[1].value[i];
	}
    cout << "[Electricity] Max normal-derivative residual: " << _maxResDn << endl;

    // 空间电场求解
#pragma omp parallel for num_threads(thread_num)
	for(int i = 0; i < mesh.vertIndexNum(); i++)
	{
		Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);

        if (vert.onCutEdge[1])
        {
            Ex.value[i] = 0.0;
            for (int j = 0; j < vert.relatedVertsNum; j++)
            {
                Ex.value[i] += vert.xCoeff.discOpr[j] * phi.value[mesh.verts(vert.relatedVerts[j]).index()];
            }
            Ex.value[i] -= vert.xCoeff.rhs;

            Ey.value[i] = 0.0;
            for (int j = 0; j < vert.relatedVertsNum; j++)
            {
                Ey.value[i] += vert.yCoeff.discOpr[j] * phi.value[mesh.verts(vert.relatedVerts[j]).index()];
            }
            Ey.value[i] -= vert.yCoeff.rhs;

            Ex.value[i] = -Ex.value[i];
            Ey.value[i] = -Ey.value[i];
        }
        else
        {
            Ex.value[i] = -FiniteDifference::dx(phi, mesh, vert);
            Ey.value[i] = -FiniteDifference::dy(phi, mesh, vert);
        }
	}
}

template<typename Mesh>
void Electricity::computeSpaceChargeDensity(Mesh& mesh, vector<ParticleType>& particleTypes, ImmersedGeometry<Mesh>& iGeometry)
{
    spaceRhoe.initialize(mesh, 0.0);
    for (int ptI = 0; ptI < particleTypes.size(); ptI++)
    {
        vector<int> backgroundCells;
    
#pragma omp parallel for num_threads(thread_num)
        for(int p = 0; p < particleTypes[ptI].number; p++)
        {
            Particle& particle(particleTypes[ptI].particle[p]);
            
            if(particle.inBackgroundCellL == -1)
            {
                particle.inBackgroundCellL = mesh.findBackgroundCellLabel(particle.position);
            }
            Cell& inBackgroundCell = mesh.cells(particle.inBackgroundCellL);
    
            particle.updateWeight(mesh, inBackgroundCell);
    
            for(int i = 0; i < 4; i++)
            {
                Vert& vert = mesh.verts(inBackgroundCell.cellToVert[i]);
                if (vert.coord().x < 0.1 * inBackgroundCell.size().x)// 轴对称边界
                {
            #pragma omp atomic
                    spaceRhoe.value[vert.index()] += particleTypes[ptI].charge() * particle.weight[i] / (0.25 * Pi * inBackgroundCell.size().x * inBackgroundCell.size().x * inBackgroundCell.size().y);
                }
                else
                {
            #pragma omp atomic
                    spaceRhoe.value[vert.index()] += particleTypes[ptI].charge() * particle.weight[i] / (2.0 * Pi * vert.coord().x * inBackgroundCell.size().x * inBackgroundCell.size().y);
                }
            }

            #pragma omp critical
            {
                backgroundCells.push_back(particle.inBackgroundCellL);
            }
        }
    
        set<int> uniqueCells(backgroundCells.begin(), backgroundCells.end());
        backgroundCells.assign(uniqueCells.begin(), uniqueCells.end());
        for (int i =0; i< backgroundCells.size();i++)
        {
            Cell& cell = mesh.cells(backgroundCells[i]);
        }
    
        Field<bool> isTreated(mesh.vertIndexNum(), false);
    
#pragma omp parallel for num_threads(thread_num)
        for(int i = 0; i < backgroundCells.size(); i++)
        {
            Cell& bgCell = mesh.cells(backgroundCells[i]);

            interpolationForSpaceRhoe(mesh, bgCell, bgCell, isTreated);
        }
    }

#pragma omp parallel for num_threads(thread_num)
    for (int i = 0; i < mesh.vertIndexNum(); i++)
    {
        Vert& vert = mesh.verts(mesh.vertIndexToLabel[i]);
        for (int levelSetI = 0; levelSetI < iGeometry.levelSets.size(); levelSetI++)
        {
            if (iGeometry.levelSets[levelSetI].value[vert.index()] > 0.0)
            {
                continue;
            }
            else
            {
                spaceRhoe.value[i] = 0.0;
            }
        }
    }
}

template<typename Mesh>
void Electricity::interpolationForSpaceRhoe(Mesh& mesh, Cell& cell, const Cell& backgroundCell, Field<bool>& isTreated)
{
    Vert& bgVert0 = mesh.verts(backgroundCell.cellToVert[0]);
    Vert& bgVert1 = mesh.verts(backgroundCell.cellToVert[1]);
    Vert& bgVert2 = mesh.verts(backgroundCell.cellToVert[2]);
    Vert& bgVert3 = mesh.verts(backgroundCell.cellToVert[3]);

	for(int i = 0; i < 4; i++)
	{
		if(cell.son[i] != -1)
		{
            Cell& sonCell = mesh.cells(cell.son[i]);
			interpolationForSpaceRhoe(mesh, sonCell, backgroundCell, isTreated);
		}
		else
		{
			Vert& vert = mesh.verts(cell.cellToVert[i]);

			if(!isTreated.value[vert.index()])
			{
				double x = (vert.coord().x - bgVert0.coord().x) / backgroundCell.size().x;
				double y = (vert.coord().y - bgVert0.coord().y) / backgroundCell.size().y;

				spaceRhoe.value[vert.index()] = spaceRhoe.value[bgVert0.index()] * (1.0 - x) * (1.0 - y)
					+ spaceRhoe.value[bgVert1.index()] * x * (1.0 - y)
					+ spaceRhoe.value[bgVert2.index()] * x * y
					+ spaceRhoe.value[bgVert3.index()] * (1.0 - x) * y;

				isTreated.value[vert.index()] = true;
			}
		}
	}
}

template<typename Mesh>
double Electricity::emitParticle(Mesh& mesh, InterfaceNode& node, ParticleType& particleType, const int levelSetI, const double& dt)
{
    Cell& cell = mesh.cells(mesh.findCellLabel(node.coord()));

    if (node.coord().x < 0.1 * cell.size().x)
    {
        accCharge[levelSetI].value[node.index()] += evapCurrent[levelSetI].value[node.index()] * Pi * node.ds() * node.ds() * dt;
    }
    else
    {
        accCharge[levelSetI].value[node.index()] += evapCurrent[levelSetI].value[node.index()] * node.ds() * 2.0 * Pi * node.coord().x * dt;
    }

	int particleNum = static_cast<int>(fabs(accCharge[levelSetI].value[node.index()] / particleType.charge()));

	accCharge[levelSetI].value[node.index()] -= static_cast<double>(particleNum) * particleType.charge();
particleNum = 0;
	particleType.number += particleNum;

	const Vector2D nodePos = node.coord() + node.normal() * MINVAL;
	const double magV = cell.size().x / dt;

    uniform_real_distribution<double> direcDistri(-0.25 * Pi, 0.25 * Pi);
	uniform_real_distribution<double> chanceDistri(MINVAL, 1);

	for(int p = 0; p < particleNum; p++)
	{
		Vector2D pos = nodePos;
		if(pos.x < 0.0)
		{
			pos.x = -pos.x;
		}

		// 计算node处法向量与x正方向的夹角
		double direc = atan2(node.normal().y, node.normal().x);
		
		double theta = direc + direcDistri(generator_);

		double vx = magV * cos(theta);
		double vy = magV * sin(theta);

		Particle tempParticle(pos.x, pos.y, vx, vy, cell.label());
		
		tempParticle.emitChance = chanceDistri(generator_);

		particleType.particle.push_back(tempParticle);
	}

	return evapCurrent[levelSetI].value[node.index()] / particleType.q_m();
}

template<typename Mesh>
void Electricity::moveParticle(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, ParticleType& particleType, const double& dt, const bool& initial)
{
	vector<int> toDelete;

#pragma omp parallel for num_threads(thread_num)
	for(int p = 0; p < particleType.number; p++)
	{
		Particle& particle(particleType.particle[p]);

		particle.inCellL = mesh.findCellLabel(particle.position);
        Cell inCell = mesh.cells(particle.inCellL);

		double remainDt = dt;
		if(initial)
			remainDt = dt * particle.emitChance;

		while(remainDt > 0.0)
		{
			particle.updateWeight(mesh, inCell);
			particle.getDiscreteVert(mesh, iGeometry.levelSets[levelSetI], levelSetI);

			Vector2D force(0.0, 0.0);

			if(particle.discreteVertL == -1)
			{
				for (int i = 0; i < 4; i++)
				{
                    Vert& vert = mesh.verts(inCell.cellToVert[i]);
					force.x += Ex.value[vert.index()] * particle.weight[i];
					force.y += Ey.value[vert.index()] * particle.weight[i];
				}
			}
			else
			{
                Vert& vert = mesh.verts(particle.discreteVertL);
				InterfaceNode& node(sharpInterface.geometry().interfaceNode(levelSetI, vert.nearestNodeLabel[levelSetI]));

				double En = 0.0;
				for (int i = 0; i < vert.relatedVertsNum; i++)
				{
					En += vert.nCoeff.discOpr[i] * phi.value[mesh.verts(vert.relatedVerts[i]).index()];
				}
				En -= vert.nCoeff.rhs;

				double Enn = 0.0;
				for (int i = 0; i < vert.relatedVertsNum; i++)
				{
					Enn += vert.nnCoeff.discOpr[i] * phi.value[mesh.verts(vert.relatedVerts[i]).index()];
				}
				Enn -= vert.nnCoeff.rhs;

				double pEn = En + ((particle.position.x - node.coord().x) * node.normal().x 
					+ (particle.position.y - node.coord().y) * node.normal().y) * Enn;

				double Et = 0.0;
				for (int i = 0; i < vert.relatedVertsNum; i++)
				{
					Et += vert.tCoeff.discOpr[i] * phi.value[mesh.verts(vert.relatedVerts[i]).index()];
				}
				Et -= vert.tCoeff.rhs;

				double Ett = 0.0;
				for (int i = 0; i < vert.relatedVertsNum; i++)
				{
					Ett += vert.ttCoeff.discOpr[i] * phi.value[mesh.verts(vert.relatedVerts[i]).index()];
				}
				Ett -= vert.ttCoeff.rhs;

				double pEt = Et + (-(particle.position.x - node.coord().x) * node.normal().y
					+ (particle.position.y - node.coord().y) * node.normal().x) * Ett;

				force.x = pEn * node.normal().x - pEt * node.normal().y;
				force.y = pEn * node.normal().y + pEt * node.normal().x;

				force.x = -force.x;
				force.y = -force.y;
			}			

			Vector2D predictVel(particle.velocity + (force * (particleType.q_m() * remainDt)));

			double localDt = inCell.size().mod() / predictVel.mod();

			double dTime = localDt < remainDt ? localDt : remainDt;

			remainDt -= localDt;

			particle.velocity += (force * (particleType.q_m() * dTime));

			particle.position += (particle.velocity * dTime);

			if(particle.position.x < 0.0)
			{
				particle.position.x = -particle.position.x;
				particle.velocity.x = -particle.velocity.x;
			}
			else if(particle.position.x > mesh.maxCoord.x)
			{
				// 标记要删除的粒子
			#pragma omp critical
				toDelete.push_back(p);

				break;
			}

			if(particle.position.y < mesh.minCoord.y || particle.position.y > mesh.maxCoord.y)
			{
				// 标记要删除的粒子
			#pragma omp critical
				toDelete.push_back(p);

				break;
			}

			particle.inCellL = mesh.findCellLabel(particle.position);
            inCell = mesh.cells(particle.inCellL);

            bool needDelete = false;
			particle.ensureLevelSetSize(iGeometry.levelSets.size());
			for (int lsI = 0; lsI < iGeometry.levelSets.size(); lsI++)
            {
                particle.levelSetValues[lsI] = 0.0;
                for(int i = 0; i < 4; i++)
                {
                    Vert& cellVert = mesh.verts(inCell.cellToVert[i]);
                    particle.levelSetValues[lsI] += iGeometry.levelSets[lsI].value[cellVert.index()] * particle.weight[i];
                }

                if (particle.levelSetValues[lsI] < 0.0)
                {
                    needDelete = true;
                    break;
                }
            }

			if(needDelete)
			{
				// 标记要删除的粒子
			#pragma omp critical
				toDelete.push_back(p);

				break;
			}
		}

		particle.inBackgroundCellL = mesh.findBackgroundCellLabel(particle.position);
	}

	// 删除所有被标记的粒子
	sort(toDelete.begin(), toDelete.end(), greater<int>());
	for(int p : toDelete)
	{
		particleType.particle.erase(particleType.particle.begin() + p);
	}

	particleType.number = particleType.particle.size();
}

template<typename Mesh>
vector<int> Electricity::findRelatedVertLabel(Mesh& mesh, Vert& vert, const int levelSetI)
{
    vector<int> relatedVertLabel;

    relatedVertLabel.push_back(vert.label());
    if (vert.neiVertices[0] != -1)
    {
        relatedVertLabel.push_back(vert.neiVertices[0]);
    }
    else
    {
        for (int i = 0; i < 4; i++)
        {
            if (vert.vertToCell[i] != -1)
            {
                relatedVertLabel.push_back(mesh.cells(vert.vertToCell[i]).cellToVert[i]);
                break;
            }
        }
    }

    for (int i = 1; i < 4; i++)
    {
        if (vert.neiVertices[i] != -1)
        {
            relatedVertLabel.push_back(vert.neiVertices[i]);
        }
        else
        {
            int oppo = (i + 2) % 4;

            relatedVertLabel.push_back(mesh.verts(vert.neiVertices[oppo]).neiVertices[oppo]);
        }
    }

    for (int i = 3; i >= 0; i--)
    {
        if (vert.vertToCell[i] != -1)
        {
            relatedVertLabel.push_back(mesh.cells(vert.vertToCell[i]).cellToVert[i]);
            break;
        }
    }
    
    return relatedVertLabel;
}

template<typename Mesh>
int Electricity::sortRelatedVertLabel(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, Vert& vert, const int levelSetI, int relatedVertLabel[6])
{
    int inPhaseNum = 0;

    vector<int> inPhaseVertLabel;
    vector<int> outPhaseVertLabel;

    double LS0 = iGeometry.levelSets[levelSetI].value[vert.index()];
    for (int i = 1; i < vert.relatedVertsNum; i++)
    {
        Vert& relatedVert = mesh.verts(relatedVertLabel[i]);
        double LS = iGeometry.levelSets[levelSetI].value[relatedVert.index()];
        if (LS * LS0 > 0.0)
        {
            inPhaseVertLabel.push_back(relatedVert.label());
        }
        else
        {
            outPhaseVertLabel.push_back(relatedVert.label());
        }
    }

    relatedVertLabel[0] = vert.label();
    for (std::size_t i = 0; i < inPhaseVertLabel.size(); i++)
    {
        relatedVertLabel[i + 1] = inPhaseVertLabel[i];
    }
    for (std::size_t i = 0; i < outPhaseVertLabel.size(); i++)
    {
        relatedVertLabel[i + inPhaseVertLabel.size() + 1] = outPhaseVertLabel[i];
    }
    vert.relatedVertsNum = static_cast<int>(inPhaseVertLabel.size() + outPhaseVertLabel.size() + 1);

    return static_cast<int>(inPhaseVertLabel.size() + 1);
}

template<typename Mesh>
void Electricity::writeCache(const Mesh& mesh, const std::string& filePath) const
{
    cout << "[Electricity] Writing electric-field cache..." << endl;

    ostringstream fileName;
	fileName << filePath << "electricity_cache.bin";

    ofstream cache(fileName.str(), ios::binary | ios::out);
	if (cache.fail())
	{
		cout << "Writing cache \"" << fileName.str() << "\" failed!" << endl;
		exit(-1);
	}

    const char magic[4] = {'E', 'L', 'C', '2'};
    const std::uint32_t version = 2;
    cache.write(reinterpret_cast<const char*>(magic), 4);
    cache.write(reinterpret_cast<const char*>(&version), sizeof(std::uint32_t));
    const int levelSetCount = static_cast<int>(iPhi.size());
    cache.write(reinterpret_cast<const char*>(&levelSetCount), sizeof(int));

    int vertFieldSize = phi.size();
    cache.write(reinterpret_cast<const char*>(&vertFieldSize), sizeof(int));
    for (int i = 0; i < vertFieldSize; i++)
    {
        cache.write(reinterpret_cast<const char*>(&phi.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&Ex.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&Ey.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&spaceRhoe.value[i]), sizeof(double));
    }

    for (int i = 0; i < levelSetCount; i++)
    {
        int nodeFieldSize = iPhi[i].size();
        cache.write(reinterpret_cast<const char*>(&nodeFieldSize), sizeof(int));
        for (int j = 0; j < nodeFieldSize; j++)
        {
            cache.write(reinterpret_cast<const char*>(&iPhi[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&EnVac[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&EnLiq[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&Et[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&expansivity[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&rhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&nextRhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&preRhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&prePreRhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&tempRhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&Vt[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&DtVt[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&accCharge[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&preNormal[i].value[j]), sizeof(Vector2D));
            cache.write(reinterpret_cast<const char*>(&prePreNormal[i].value[j]), sizeof(Vector2D));
            cache.write(reinterpret_cast<const char*>(&condCurrent[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&convCurrent[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&evapCurrent[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&dtRhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&nStressMaxwell[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&tStressMaxwell[i].value[j]), sizeof(double));
        }

        int preNodeFieldSize = preInterfaceRhoe[i].size();
        cache.write(reinterpret_cast<const char*>(&preNodeFieldSize), sizeof(int));
        for (int j = 0; j < preNodeFieldSize; j++)
        {
            cache.write(reinterpret_cast<const char*>(&preInterfaceRhoe[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&preAccCharge[i].value[j]), sizeof(double));
        }

        int prePreNodeFieldSize = prePreInterfaceRhoe[i].size();
        cache.write(reinterpret_cast<const char*>(&prePreNodeFieldSize), sizeof(int));
        for (int j = 0; j < prePreNodeFieldSize; j++)
        {
            cache.write(reinterpret_cast<const char*>(&prePreInterfaceRhoe[i].value[j]), sizeof(double));
        }
    }
    cache.write(reinterpret_cast<const char*>(&generator_), sizeof(std::default_random_engine));

    cache.close();
    cout << "[Electricity] Electric-field cache written." << endl;
}

template<typename Mesh>
void Electricity::readCache(Mesh& mesh, const std::string& filePath)
{
    cout << "[Electricity] Reading electric-field cache..." << endl;

    ostringstream fileName;
    fileName << filePath << "electricity_cache.bin";

    ifstream cache(fileName.str(), ios::binary | ios::in);
    if (cache.fail())
    {
        cout << "Reading cache \"" << fileName.str() << "\" failed!" << endl;
        exit(-1);
    }

    char magic[4] = {0, 0, 0, 0};
    cache.read(magic, 4);
    const char expected[4] = {'E', 'L', 'C', '2'};
    if (std::memcmp(magic, expected, 4) != 0)
    {
        cout << "[Electricity] ERROR: incompatible cache format in \"" << fileName.str()
             << "\" (delete cache)." << endl;
        exit(-1);
    }
    std::uint32_t version = 0;
    cache.read(reinterpret_cast<char*>(&version), sizeof(std::uint32_t));
    if (version != 2u)
    {
        cout << "[Electricity] ERROR: unsupported cache version " << version
             << " in \"" << fileName.str() << "\"." << endl;
        exit(-1);
    }
    int levelSetCount = 0;
    cache.read(reinterpret_cast<char*>(&levelSetCount), sizeof(int));
    if (levelSetCount < 0)
        levelSetCount = 0;
    ensureLevelSetCount(levelSetCount);

    int vertFieldSize;
    cache.read(reinterpret_cast<char*>(&vertFieldSize), sizeof(int));
    phi.initialize(mesh, 0.0);
    Ex.initialize(mesh, 0.0);
    Ey.initialize(mesh, 0.0);
    spaceRhoe.initialize(mesh, 0.0);

    for (int i = 0; i < vertFieldSize; i++)
    {
        cache.read(reinterpret_cast<char*>(&phi.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&Ex.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&Ey.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&spaceRhoe.value[i]), sizeof(double));
    }

    for (int i = 0; i < levelSetCount; i++)
    {
        int nodeFieldSize;
        cache.read(reinterpret_cast<char*>(&nodeFieldSize), sizeof(int));
        iPhi[i].initialize(nodeFieldSize, 0.0);
        EnVac[i].initialize(nodeFieldSize, 0.0);
        EnLiq[i].initialize(nodeFieldSize, 0.0);
        Et[i].initialize(nodeFieldSize, 0.0);
        expansivity[i].initialize(nodeFieldSize, 0.0);
        rhoe[i].initialize(nodeFieldSize, 0.0);
        nextRhoe[i].initialize(nodeFieldSize, 0.0);
        preRhoe[i].initialize(nodeFieldSize, 0.0);
        prePreRhoe[i].initialize(nodeFieldSize, 0.0);
        tempRhoe[i].initialize(nodeFieldSize, 0.0);
        Vt[i].initialize(nodeFieldSize, 0.0);
        DtVt[i].initialize(nodeFieldSize, 0.0);
        accCharge[i].initialize(nodeFieldSize, 0.0);
        preNormal[i].initialize(nodeFieldSize, Vector2D(1.0, 1.0));
        prePreNormal[i].initialize(nodeFieldSize, Vector2D(1.0, 1.0));
        condCurrent[i].initialize(nodeFieldSize, 0.0);
        convCurrent[i].initialize(nodeFieldSize, 0.0);
        evapCurrent[i].initialize(nodeFieldSize, 0.0);
        dtRhoe[i].initialize(nodeFieldSize, 1.0);
        nStressMaxwell[i].initialize(nodeFieldSize, 1.0);
        tStressMaxwell[i].initialize(nodeFieldSize, 1.0);

        for (int j = 0; j < nodeFieldSize; j++)
        {
            cache.read(reinterpret_cast<char*>(&iPhi[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&EnVac[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&EnLiq[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&Et[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&expansivity[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&rhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&nextRhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&preRhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&prePreRhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&tempRhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&Vt[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&DtVt[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&accCharge[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&preNormal[i].value[j]), sizeof(Vector2D));
            cache.read(reinterpret_cast<char*>(&prePreNormal[i].value[j]), sizeof(Vector2D));
            cache.read(reinterpret_cast<char*>(&condCurrent[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&convCurrent[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&evapCurrent[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&dtRhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&nStressMaxwell[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&tStressMaxwell[i].value[j]), sizeof(double));
        }

        int preNodeFieldSize;
        cache.read(reinterpret_cast<char*>(&preNodeFieldSize), sizeof(int));
        preInterfaceRhoe[i].initialize(preNodeFieldSize, 0.0);
        preAccCharge[i].initialize(preNodeFieldSize, 0.0);
        for (int j = 0; j < preNodeFieldSize; j++)
        {
            cache.read(reinterpret_cast<char*>(&preInterfaceRhoe[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&preAccCharge[i].value[j]), sizeof(double));
        }

        int prePreNodeFieldSize;
        cache.read(reinterpret_cast<char*>(&prePreNodeFieldSize), sizeof(int));
        prePreInterfaceRhoe[i].initialize(prePreNodeFieldSize, 0.0);
        for (int j = 0; j < prePreNodeFieldSize; j++)
        {
            cache.read(reinterpret_cast<char*>(&prePreInterfaceRhoe[i].value[j]), sizeof(double));
        }
    }

    cache.read(reinterpret_cast<char*>(&generator_), sizeof(std::default_random_engine));

    cache.close();
    cout << "[Electricity] Electric-field cache read." << endl;
}
