#pragma once

#include "Modules/Interface/HermiteGalsPhaseFlux.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Interface
{
class HermiteGalsIntegratedInventoryResult;
struct HermiteGalsSideFlux
{
    int cell = -1; // native cell index
    FaceSide side = FaceSide::West;
    int outwardSign = -1;
    std::size_t subfaceCount = 0;
    // Indices into source().faces(), ordered along the positive tangent.
    // One actual face on a conforming/fine side, two on a 2:1 coarse side.
    std::array<std::size_t, 2> subfaces{};
    // All rates retain GLOBAL +x/+y orientation. Use outwardSign for the
    // cell-outward signed rate; swap +/- for outward gross rates when -1.
    long double totalPositiveRate = 0, totalNegativeRate = 0;
    long double liquidPositiveRate = 0, liquidNegativeRate = 0;
    long double totalRate = 0, liquidRate = 0;
    long double liquidVolume = 0; // sum of cached dt*Q_l, never dt*sum(Q_l)
    long double integrationErrorEstimate = 0;
};
struct HermiteGalsRefluxReport
{
    std::size_t cellSides = 0, coarseSides = 0, interiorCoarseSides = 0;
    std::size_t coarseFinePairs = 0;
    // Max |coarse outward + sum(fine outward)|, unnormalized absolute units.
    long double maxTotalRateResidual = 0, maxLiquidRateResidual = 0;
    long double maxLiquidVolumeResidual = 0;
};
class HermiteGalsRefluxResult
{
public:
    // Same immutable authority: faces, M, geometry, identities and time.
    const HermiteGalsPhaseFluxResult& source() const;
    const std::vector<HermiteGalsSideFlux>& sides() const; // cell ascending, W/S/E/N
    const HermiteGalsRefluxReport& report() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsRefluxResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsRefluxResult assembleHermiteGalsReflux(
        const QuadtreeMesh&, const HermiteGalsPhaseFluxResult&);
};
// Single shared midpoint interval, including I08's final roundoff repair.
// No coarse-face quadrature or second correction of the already refluxed M.
HermiteGalsRefluxResult assembleHermiteGalsReflux(
    const QuadtreeMesh& mesh, const HermiteGalsPhaseFluxResult& flux);

struct HermiteGalsIntegratedSideTransfer
{
    int cell=-1;
    FaceSide side=FaceSide::West;
    int outwardSign=-1;
    std::size_t subfaceCount=0;
    std::array<std::size_t,2> subfaces{}; // indices into source().faces()
    // Global +x/+y orientation. Derived sums of FINAL actual subfaces.
    long double totalPositiveVolume=0,totalNegativeVolume=0;
    long double liquidPositiveVolume=0,liquidNegativeVolume=0;
    // Sum each native signed volume in tangent order, preserving the ledger
    // grouping. Never recompute a conserved volume from an averaged rate.
    long double totalVolume=0,liquidVolume=0;
};
struct HermiteGalsIntegratedRefluxReport
{
    std::size_t cellSides=0,coarseSides=0,interiorCoarseSides=0,coarseFinePairs=0;
    long double maxTotalVolumeResidual=0,maxLiquidVolumeResidual=0;
};
class HermiteGalsIntegratedRefluxResult
{
public:
    const HermiteGalsIntegratedInventoryResult& source() const;
    const std::vector<HermiteGalsIntegratedSideTransfer>& sides() const;
    const HermiteGalsIntegratedRefluxReport& report() const;
    long double averageTotalFlux(std::size_t sideRow) const;
    long double averageLiquidFlux(std::size_t sideRow) const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsIntegratedRefluxResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsIntegratedRefluxResult assembleHermiteGalsReflux(
        const QuadtreeMesh&,const HermiteGalsIntegratedInventoryResult&);
};
// Same topology/ownership assembly as I09; no integration or M correction.
HermiteGalsIntegratedRefluxResult assembleHermiteGalsReflux(
    const QuadtreeMesh& mesh,const HermiteGalsIntegratedInventoryResult& flux);
}
