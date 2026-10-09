#pragma once

#include "Modules/Hydrodynamics/CchCharacteristic.h"
#include "Modules/Interface/HermiteGalsFirstJetStep.h"
#include <array>
#include <memory>

namespace Quantum::Hydrodynamics
{
struct CchMaterialSources
{
    Interface::HermiteGalsIdentity begin, midpoint, end;
    CchAdvectorIdentity advector; // Actual spatial owner; zero only for temporary legacy consumers.
};
struct CchMaterialCell
{
    long double massBegin=0, massMidpoint=0, massEnd=0;
    long double massChange=0, outwardIntegral=0, integratedDefect=0, defectRate=0;
    long double relativeDefect=0; // |defect|/(|massChange|+sum |face integral|).
};
struct CchMaterialReport
{
    std::array<long double,3> geometryMass{{0,0,0}};
    long double massChange=0, boundaryIntegral=0, integratedDefect=0, l1IntegratedDefect=0;
    long double maxAbsDefectRate=0, maxRelativeCellDefect=0, balanceIdentity=0;
};
struct CchSpaceTimeConfig
{
    CchCharacteristicConfig characteristic;
    double relativeTolerance=1e-7;
    int initialPanels=4,maxDepth=24,maxGeometrySamples=8192;
};
struct CchSpaceTimeReport
{
    std::size_t geometrySamples=0,acceptedPanels=0;
    int deepestPanel=0;
    long double absolutePhaseErrorEstimate=0,maxRelativePhaseErrorEstimate=0;
};
// UNIQUE instantaneous density-weighted spatial flux. Valid at any declared
// Initial/Midpoint/End layer; no dt, inventory, mass update or accepted state.
class CchMaterialFlux
{
public:
    static CchMaterialFlux prepare(const QuadtreeMesh&, CchPhaseRate,
        double liquidDensity, double gasDensity);
    const CchPhaseRate& phaseRate() const noexcept;
    const CchMomentumGeometry& momentumGeometry() const;
    const CchDofLayout& layout() const noexcept;
    CchAdvectorIdentity advectorIdentity() const;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    double liquidDensity() const noexcept;
    double gasDensity() const noexcept;
    const std::vector<long double>& atomRates() const noexcept;
    const std::vector<long double>& faceRates() const noexcept;
    const std::vector<long double>& dualFaceRates() const noexcept;
    // One-sided point view of this SAME density-weighted physical flux.
    // Phase is supplied by the retained publication/cubature.
    Vector2D physicalFlux(int pressure,Vector2D point,bool liquid) const;
    void validateFor(const QuadtreeMesh&, const CchMomentumGeometry& expected,
        Interface::HermiteGalsIdentity, CchAdvectorIdentity,
        Interface::HermiteGalsTimeLayer) const;
private:
    struct Data;
    explicit CchMaterialFlux(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
// One fixed-layout interval: measured geometry masses + UNIQUE material
// integral + explicit nonconservative defect. No update/commit or repair.
// prepare keeps the original midpoint rule; prepareSpaceTime integrates the
// actual moving geometry and endpoint source over the whole interval.
class CchMaterialTransport
{
public:
    static CchMaterialTransport prepare(const QuadtreeMesh&, CchPhaseGeometry begin,
        CchPhaseRate midpoint, CchPhaseGeometry end, double liquidDensity, double gasDensity);
    static CchMaterialTransport prepare(const QuadtreeMesh&, CchPhaseGeometry begin,
        CchMaterialFlux midpoint, CchPhaseGeometry end);
    // True spacetime integration along ONE original FirstJet/endpoint source.
    // Every sampled geometry uses that original state, not interpolated snapshots.
    // If supplied, boundary prescribes the full outer arrival FirstJet (including
    // tangent/outflow sides and corners); it is not an exterior scalar extension.
    static CchMaterialTransport prepareSpaceTime(const QuadtreeMesh&,
        Interface::HermiteGalsPublication original,CchCharacteristicSource,
        double liquidDensity,double gasDensity,CchSpaceTimeConfig,
        const Interface::HermiteGalsBoundaryFirstJetQuery& boundary={},
        const Interface::HermiteGalsPublicationConfig& publication={});
    bool isSpaceTime() const noexcept;
    const CchCharacteristicSource& characteristicSource() const;
    const CchSpaceTimeReport& spaceTimeReport() const;
    const std::vector<CchPhaseMeasure>& phaseAtomIntegrals() const;
    const std::vector<long double>& atomIntegrals() const noexcept;
    const std::vector<long double>& dualFaceIntegrals() const noexcept;
    const CchPhaseGeometry& beginGeometry() const noexcept;
    const CchPhaseRate& midpointRate() const noexcept;
    const CchPhaseGeometry& endGeometry() const noexcept;
    const CchDofLayout& layout() const noexcept;
    CchMaterialSources sources() const;
    Interface::HermiteGalsTimeLayer timeLayer() const; // Full Midpoint clock.
    double deltaTime() const noexcept;
    double liquidDensity() const noexcept;
    double gasDensity() const noexcept;
    bool isSpatial() const noexcept;
    // Original spatial geometry and single density-weighted atom owner.
    // These ports reject temporary legacy preparations.
    const CchMomentumGeometry& momentumGeometry() const;
    const CchMaterialFlux& midpointFlux() const;
    const std::vector<long double>& atomRates() const;
    const std::vector<long double>& dualFaceRates() const;
    const std::vector<long double>& faceRates() const noexcept; // interval average for spacetime; midpoint otherwise.
    const std::vector<long double>& faceIntegrals() const noexcept;
    const std::vector<CchMaterialCell>& cells() const noexcept; // pressure DOF order.
    const CchMaterialReport& report() const noexcept;
    void validateFor(const QuadtreeMesh&, CchMaterialSources,
        Interface::HermiteGalsTimeLayer midpointClock) const;
    void validateFor(const QuadtreeMesh&, const CchMomentumGeometry& expected,
        CchMaterialSources, Interface::HermiteGalsTimeLayer midpointClock) const;
private:
    static CchMaterialTransport assemble(const QuadtreeMesh&, CchPhaseGeometry,
        CchPhaseRate, CchPhaseGeometry, double, double, const CchMaterialFlux*,
        const std::vector<CchPhaseMeasure>* integrated=nullptr,const CchCharacteristicSource* motion=nullptr,
        const CchSpaceTimeReport* report=nullptr);
    struct Data;
    explicit CchMaterialTransport(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
