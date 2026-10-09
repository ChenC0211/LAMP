#pragma once

#include "Modules/Hydrodynamics/CchPhaseGeometry.h"
#include <limits>

namespace Quantum::Hydrodynamics
{
struct CchRegularizedColorConfig
{
    double epsilon=0; // ONE explicitly supplied physical length, not a cellwise h factor.
    double relativeTolerance=1e-9;
    double relativeGradientResolution=256*std::numeric_limits<double>::epsilon();
    int maxDepth=20;
    std::size_t maxQueriesPerCell=262144; // Total original jet and affine-enclosure queries.
};
struct CchRegularizedColorCellReport
{
    std::size_t queries=0,enclosureQueries=0,acceptedPanels=0;
    int deepestPanel=0;
    long double absoluteErrorEstimate=0,transitionVolumeEstimate=0;
    long double transitionVolumeLower=0,transitionVolumeUpper=0;
    double gradientScale=0,minimumResolvedGradientRatio=std::numeric_limits<double>::infinity();
    double minimumCertifiedGradientRatio=std::numeric_limits<double>::infinity();
    bool certifiedPure=false;
};
struct CchRegularizedColorReport
{
    std::size_t sourceQueries=0,enclosureQueries=0,pureCells=0;
    long double physicalVolume=0,absoluteErrorEstimate=0;
    long double transitionVolumeEstimate=0,transitionVolumeLower=0,transitionVolumeUpper=0;
    long double colorMinusSharpL1=0;
    double colorMinusSharpMax=0;
};
// ONE read-only force color derived from the retained FirstJet publication.
// No scalar/material inventory, transport, pressure/velocity, IO or main clock.
class CchRegularizedColor
{
public:
    static CchRegularizedColor build(const QuadtreeMesh&,CchPhaseGeometry,CchRegularizedColorConfig);
    static double heaviside(double signedLength,double epsilon);
    static double delta(double signedLength,double epsilon);
    const CchPhaseGeometry& phaseGeometry() const noexcept;
    const CchRegularizedColorConfig& config() const noexcept;
    const std::vector<double>& colors() const noexcept; // native pressure order
    const std::vector<CchRegularizedColorCellReport>& cellReports() const noexcept;
    const CchRegularizedColorReport& report() const noexcept;
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    void validateFor(const QuadtreeMesh&,const CchPhaseGeometry& expected) const;
    void validateFor(const QuadtreeMesh&,Interface::HermiteGalsIdentity,Interface::HermiteGalsTimeLayer) const;
    double pointColor(int pressure,Vector2D point) const;
private:
    struct Data;
    explicit CchRegularizedColor(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
