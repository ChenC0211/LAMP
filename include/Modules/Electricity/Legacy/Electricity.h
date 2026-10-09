#pragma once

#include <algorithm>
#include <iostream>
#include <random>
#include <set>
#include <vector>
#include "Core/Base/Base.h"
#include "Modules/Hydrodynamics/LegacyNodalState.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/Vert.h"
#include "Core/Field/VertField.h"
#include "Core/Field/CellField.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

class Electricity
{
public:
    VertField<double> phi;
    VertField<double> Ex;
    VertField<double> Ey;
    VertField<double> spaceRhoe;
    VertField<double> res;

    std::vector<Field<double>> iPhi;
    std::vector<Field<double>> EnVac;
    std::vector<Field<double>> EnLiq;
    std::vector<Field<double>> Et;
    std::vector<Field<double>> expansivity;
    std::vector<Field<double>> rhoe;// 电荷密度
    std::vector<Field<double>> nextRhoe;// 下一步电荷密度
    std::vector<Field<double>> preRhoe;// 当前节点的上一步电荷密度，由前一步界面的电荷preInterfaceRhoe插值得到
    std::vector<Field<double>> prePreRhoe;// 当前节点的上上一步电荷密度，由上上一步界面的电荷prePreInterfaceRhoe插值得到
    std::vector<Field<double>> tempRhoe;
    std::vector<Field<double>> Vt;
    std::vector<Field<double>> DtVt;
    std::vector<Field<double>> accCharge;
    std::vector<Field<double>> preAccCharge;
    std::vector<Field<Vector2D>> preNormal;
    std::vector<Field<Vector2D>> prePreNormal;
    std::vector<Field<double>> preInterfaceRhoe;// 前一步界面的电荷
    std::vector<Field<double>> prePreInterfaceRhoe;// 上上一步界面的电荷
    std::vector<Field<double>> condCurrent;// 传导电流
    std::vector<Field<double>> convCurrent;// 对流电流
    std::vector<Field<double>> evapCurrent;// 蒸发电流
    std::vector<Field<double>> dtRhoe;// 表面电荷密度变化速率
    std::vector<Field<double>> nStressMaxwell;// Maxwell法向应力
    std::vector<Field<double>> tStressMaxwell;// Maxwell切向应力

    template<typename Mesh>
    Electricity(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    ~Electricity();

    template<typename Mesh>
    void initialize(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void solve(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry, std::vector<ParticleType>& particleTypes, VertField<double>& Vx, VertField<double>& Vy, VertField<double>& preVx, VertField<double>& preVy, const double& dt, const double& preDt, const int _loop);

    template<typename Mesh>
    void writeCache(const Mesh& mesh, const std::string& filePath) const;
    template<typename Mesh>
    void readCache(Mesh& mesh, const std::string& filePath);

private:
    std::vector<Eigen::Triplet<double>> phiTrip_;
	Eigen::VectorXd b_phi_;

    double k_ = 3.0;
    double beta_liq_ = 1.0;
    double beta_vac_ = 1.0;

    Field<array<double, 6>> dlap_;
    Field<array<double, 6>> dself_;
    Field<array<double, 6>> dn_;
    Field<array<double, 6>> dt_;
    Field<array<double, 6>> dx_;
    Field<array<double, 6>> dy_;

    std::default_random_engine generator_;

    double rhoeRelaxCoeff_ = 0.03;// 松弛系数

    int levelSetCount_ = 0;

    void ensureLevelSetCount(int n)
    {
        levelSetCount_ = std::max(0, n);
        iPhi.resize(levelSetCount_);
        EnVac.resize(levelSetCount_);
        EnLiq.resize(levelSetCount_);
        Et.resize(levelSetCount_);
        expansivity.resize(levelSetCount_);
        rhoe.resize(levelSetCount_);
        nextRhoe.resize(levelSetCount_);
        preRhoe.resize(levelSetCount_);
        prePreRhoe.resize(levelSetCount_);
        tempRhoe.resize(levelSetCount_);
        Vt.resize(levelSetCount_);
        DtVt.resize(levelSetCount_);
        accCharge.resize(levelSetCount_);
        preAccCharge.resize(levelSetCount_);
        preNormal.resize(levelSetCount_);
        prePreNormal.resize(levelSetCount_);
        preInterfaceRhoe.resize(levelSetCount_);
        prePreInterfaceRhoe.resize(levelSetCount_);
        condCurrent.resize(levelSetCount_);
        convCurrent.resize(levelSetCount_);
        evapCurrent.resize(levelSetCount_);
        dtRhoe.resize(levelSetCount_);
        nStressMaxwell.resize(levelSetCount_);
        tStressMaxwell.resize(levelSetCount_);
    }

    template<typename Mesh>
    void computeInterfaceVar(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, VertField<double>& Vx, VertField<double>& Vy, VertField<double>& preVx, VertField<double>& preVy, const double& dt, const double& preDt, const int _loop);

    void computeInterfaceRhoe(LegacyNodalState& sharpInterface, const int levelSetI, const double& dt, const double& preDt, const int _loop);
    template<typename Mesh>
    void buildEleMatrixRZ(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void buildEleRHSRZ(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    std::vector<int> findRelatedVertLabel(Mesh& mesh, Vert& vert, const int levelSetI);
    template<typename Mesh>
    int sortRelatedVertLabel(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, Vert& vert, const int levelSetI, int relatedVertLabel[6]);
    template<typename Mesh>
    void computeInterfaceCoeff(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void computeInterfaceCoeff_RHS(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void computeElectricField(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void computeSpaceChargeDensity(Mesh& mesh, std::vector<ParticleType>& particleTypes, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void interpolationForSpaceRhoe(Mesh& mesh, Cell& cell, const Cell& backgroundCell, Field<bool>& isTreated);
    template<typename Mesh>
    double emitParticle(Mesh& mesh, InterfaceNode& node, ParticleType& particleType, const int levelSetI, const double& dt);
    template<typename Mesh>
    void moveParticle(Mesh& mesh, LegacyNodalState& sharpInterface, ImmersedGeometry<Mesh>& iGeometry, const int levelSetI, ParticleType& particleType, const double& dt, const bool& initial);
};

#include "Modules/Electricity/Legacy/Electricity_impl.h"
