#pragma once

#include "Modules/Interface/VofSmoothGeometry.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/Base/Vector2D.h"
#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace LAMP::Electric
{
class ElectricVofSharpField;
class ElectricSurfaceChargeTransport;
// Independent signed surface moments, ordered by the actual smooth arcs.
// moment[k]=integral qs*P_k(2t-1) dS; moment[0] is integrated charge.
class ElectricVofSurfaceCharge
{
public:
    static ElectricVofSurfaceCharge fromMoments(const QuadtreeMesh&,
        const Interface::VofSmoothGeometry&,LAMP::TimeLayer,
        const std::vector<std::array<double,3>>&);
    static ElectricVofSurfaceCharge fromDensity(const QuadtreeMesh&,
        const Interface::VofSmoothGeometry&,
        const std::function<double(std::size_t,double)>&);
    const Interface::VofSmoothGeometry& geometry() const;
    LAMP::TimeLayer timeLayer() const;
    const std::vector<std::array<double,3>>& moments() const;
    const std::vector<double>& integratedCharge() const;
    bool sameInstance(const ElectricVofSurfaceCharge&) const noexcept;
    void validateFor(const QuadtreeMesh&,const Interface::VofSmoothGeometry&,LAMP::TimeLayer) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofSurfaceCharge(std::shared_ptr<const Data>);
    static ElectricVofSurfaceCharge fromPrepared(const QuadtreeMesh&,
        const Interface::VofSmoothGeometry&,LAMP::TimeLayer,
        std::vector<std::array<double,3>>);
    friend class ElectricVofSharpField;
    friend class ElectricSurfaceChargeTransport;
};
struct ElectricVofSharpBoundary
{
    enum class Type { Dirichlet,HomogeneousNeumann };
    Type type=Type::Dirichlet;
    std::function<double(Vector2D)> value;
};
struct ElectricVofSharpProblem
{
    std::array<double,2> epsilon{{1,1}},conductivity{{0,0}};
    std::array<ElectricVofSharpBoundary,4> boundaries; // W,S,E,N.
};
struct ElectricVofSharpControls
{
    double relativeLinearTolerance=1e-11;
    double relativeConservationTolerance=1e-10;
};
struct ElectricVofSharpInterfaceField
{
    Vector2D normal,tangent; // n: liquid->gas; t: increasing actual arc parameter.
    std::array<double,2> normalElectricField{};
    double tangentialElectricField=0,surfaceChargeDensity=0;
    std::array<Vector2D,2> electricField{};
};
struct ElectricVofSharpReport
{
    std::size_t elements=0,sharedFaces=0,interfaceFaces=0,skeletonDofs=0,nonzeros=0;
    double relativeLinearResidual=0,relativeSymmetryDefect=0;
    double maximumRelativeGaussResidual=0,maximumRelativeChargeUpdateResidual=0;
    double maximumRelativeLocalResidual=0,maximumRelativeNativeBalance=0;
    double globalGaussResidual=0,relativeGlobalGaussResidual=0;
    double boundaryDisplacement=0,boundaryCurrent=0;
    double globalChargeUpdateResidual=0,relativeGlobalChargeUpdateResidual=0;
    double assembleSeconds=0,factorSeconds=0,solveSeconds=0,recoverySeconds=0;
    unsigned defectCorrections=0;
    // Correction branch wall time is a subset of solve/recovery, not additive.
    double defectCorrectionSeconds=0,initialMaximumRelativeGaussResidual=0;
    // Complete stored-C cell/free-trace residual; unavailable unless measured.
    bool defectResidualMeasured=false;
    double uncondensedResidualBefore=0,uncondensedResidualAfter=0;
    bool usedSuperLU=false;
};
// Real two-phase k=2 HHO on frozen B1 smooth geometry, no bulk charge.
// No mesh mutation, clock advance, transport, IO or fluid commit occurs here.
class ElectricVofSharpField
{
public:
    static ElectricVofSharpField prepareElectrostatic(const QuadtreeMesh&,
        const Interface::VofSmoothGeometry&,const ElectricVofSurfaceCharge&,
        const ElectricVofSharpProblem&,const ElectricVofSharpControls& = {});
    static ElectricVofSharpField prepareOhmic(const QuadtreeMesh&,
        const Interface::VofSmoothGeometry&,const ElectricVofSurfaceCharge&,
        LAMP::TimeLayer candidateClock,const ElectricVofSharpProblem&,
        const ElectricVofSharpControls& = {});
    static ElectricVofSharpField prepareOhmic(const QuadtreeMesh&,
        const ElectricSurfaceChargeTransport&,const ElectricVofSharpProblem&,
        const ElectricVofSharpControls& = {});
    const Interface::VofSmoothGeometry& geometry() const;
    const ElectricVofSurfaceCharge& sourceCharge() const;
    const ElectricVofSurfaceCharge& chargeCandidate() const;
    const ElectricVofSurfaceCharge& implicitInputCharge() const;
    const ElectricSurfaceChargeTransport* advection() const;
    LAMP::TimeLayer candidateTime() const;
    double effectiveDt() const;
    const ElectricVofSharpProblem& materials() const;
    const ElectricVofSharpReport& report() const;
    Vector2D sampleElectricField(int nativeCell,unsigned phase,Vector2D) const;
    // Query site is the exact sum of its retained anchor/offset. Cell bounds
    // and the full site-minus-basis-center difference precede native RN.
    Vector2D sampleElectricFieldAt(int nativeCell,unsigned phase,
        Geometry::PlanarCellCut::AnchoredPoint) const;
    ElectricVofSharpInterfaceField interfaceField(std::size_t arc,double t) const;
    bool sameInstance(const ElectricVofSharpField&) const noexcept;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit ElectricVofSharpField(std::shared_ptr<const Data>);
    static ElectricVofSharpField prepare(const QuadtreeMesh&,
        const Interface::VofSmoothGeometry&,const ElectricVofSurfaceCharge&,
        LAMP::TimeLayer,const ElectricVofSharpProblem&,const ElectricVofSharpControls&,bool,
        const ElectricSurfaceChargeTransport* = nullptr);
};
}
