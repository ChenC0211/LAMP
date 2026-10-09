#pragma once

#include "Core/Lagrangian/LegacyNodalGeometry.h"
#include "Core/Field/Field.h"
#include "Core/Field/VertField.h"

// Compatibility owner for the legacy nodal Hydro/IBM state and SIC2 payload.
// Modern interface/traction services do not depend on this mutable state.
class LegacyNodalState
{
public:
    LegacyNodalGeometry& geometry() { return geometry_; }
    const LegacyNodalGeometry& geometry() const { return geometry_; }

    // Indexed by geometry node label, including inactive nodes, as in SIC2.
    const std::vector<double>& evaporationVelocity(int levelSetI, Time layer = Time::current) const;
    void setEvaporationVelocity(int levelSetI, int nodeLabel, double value);
    double noSlipEvaporationVelocity(int nodeLabel) const;
    void setNoSlipEvaporationVelocity(int nodeLabel, double value);

    std::vector<Field<double>> Uxs;
    std::vector<Field<double>> Uys;
    std::vector<Field<double>> preUxs;
    std::vector<Field<double>> preUys;
    std::vector<Field<double>> prePreUxs;
    std::vector<Field<double>> prePreUys;
    std::vector<Field<double>> Fxs;
    std::vector<Field<double>> Fys;
    Field<double> noSlipUx;
    Field<double> noSlipUy;
    Field<double> noSlipFx;
    Field<double> noSlipFy;

    LegacyNodalState();
    ~LegacyNodalState();

    // Preserve the old geometry-to-velocity update order; callers initialize loads.
    // Velocity history is advanced explicitly by legacy callers, as before.
    template<typename Mesh>
    void initializeGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry);
    template<typename Mesh>
    void updateGeometry(Mesh& mesh, ImmersedGeometry<Mesh>& iGeometry, const VertField<double>& ux, const VertField<double>& uy, const VertField<double>& preUx, const VertField<double>& preUy, const double& dt, const double& preDt, const int loop);
    template<typename Field>
    double interfaceInterpolation(const Vector2D& coord, Field& field, const int levelSetI, const bool isTimeDiff, const bool isVt, const Time& timeStep, const Vector2D& normal);

    void output(const std::string& _name, const int levelSetI, const std::vector<std::string>& valNames = {}, const std::vector<Field<double>>& fields = {}) const;
    void temporalOutput(const std::string& _name, const std::vector<std::string>& valNames, const std::vector<Field<double>>& fields, const int levelSetI, const InterfaceNode& node, double nodeEvaporationVelocity, const double& t, const int loop, const bool needInitialize) const;

    void writeCache(const std::string& filePath) const;
    void readCache(const std::string& filePath);

private:
    LegacyNodalGeometry geometry_;
    std::vector<std::vector<double>> evaporationVelocity_;
    std::vector<std::vector<double>> preEvaporationVelocity_;
    std::vector<std::vector<double>> prePreEvaporationVelocity_;
    std::vector<double> noSlipEvaporationVelocity_;
    void ensureLevelSetCount(int n)
    {
        const int count = std::max(0, n);
        geometry_.resize(count);
        evaporationVelocity_.resize(count);
        preEvaporationVelocity_.resize(count);
        prePreEvaporationVelocity_.resize(count);
        Uxs.resize(count);
        Uys.resize(count);
        preUxs.resize(count);
        preUys.resize(count);
        prePreUxs.resize(count);
        prePreUys.resize(count);
        Fxs.resize(count);
        Fys.resize(count);
    }

    template<typename Mesh>
    void computeVelocity(Mesh& mesh, const int levelSetI, const VertField<double>& ux, const VertField<double>& uy, const VertField<double>& preUx, const VertField<double>& preUy, const double& dt, const double& preDt, const int loop);
};

#include "Modules/Hydrodynamics/LegacyNodalState_impl.h"
