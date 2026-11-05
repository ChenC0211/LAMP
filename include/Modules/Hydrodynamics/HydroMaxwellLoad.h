#pragma once
#include "Modules/Hydrodynamics/HydroPressureLoad.h"
#include "Modules/Hydrodynamics/HydroPhaseGeometry.h"
#include "Modules/Hydrodynamics/HydroEnergyViscosity.h"
#include "Modules/Electric/ElectricVofMaxwellStress.h"
#include "Modules/Electric/ElectricVofSharpMaxwellStress.h"
#include <memory>
#include <vector>

namespace LAMP::Hydrodynamics
{
enum class HydroMaxwellPairing { CellMeanStrain, Rt0Green };
struct HydroMaxwellLoadReport
{
    // A: v_raw div(T) - sum_sides v_shared T n + physical exterior T n.
    // B1 retains the partitioned sharp cell-mean strain representation.
    HydroMaxwellPairing pairing=HydroMaxwellPairing::CellMeanStrain;
    std::size_t cellSidePoints=0;
    std::size_t volumePoints=0,boundaryPoints=0,axisPoints=0;
    // This load's successful exact compiled contraction, not cached field data.
    bool compiledContraction=false,compiledGeometryHit=false,compiledTimingAvailable=false;
    std::size_t compiledVolumeTerms=0,compiledBoundaryTerms=0,compiledRetainedBytes=0;
    // Nested owner wall intervals; sampling and contraction share one original
    // team interval (no extra worker barriers). Do not add them to owner total.
    double validationSeconds=0,compiledPlanSeconds=0,samplingAndContractionSeconds=0,finalRneSeconds=0;
};
// Actual typed Electric stress -> native force. Does not solve, advance,
// accept, alter viscosity, or manufacture a caller-provided RHS array.
class HydroMaxwellLoad
{
public:
    static HydroMaxwellLoad prepare(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,Electric::ElectricVofMaxwellStress);
    // A-only request: retain actual RT0 volume/face velocity geometry without
    // constructing the unused full strain representation.
    static HydroMaxwellLoad prepareRt0Green(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,Electric::ElectricVofMaxwellStress);
    // Reuse only immutable RT0 Green-capable geometry retaining the actual
    // mesh/layout/divergence owner. Phase, stress and load are freshly qualified.
    static HydroMaxwellLoad prepare(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,Electric::ElectricVofMaxwellStress,HydroStrainOperator);
    static HydroMaxwellLoad prepareSharp(const QuadtreeMesh&,HydroPressureLoad,
        HydroPhaseGeometry,Electric::ElectricVofSharpMaxwellStress);
    const HydroPressureLoad& pressureLoad() const;
    const HydroPhaseGeometry& phaseGeometry() const;
    const Electric::ElectricVofMaxwellStress& stress() const;
    bool isSharp() const noexcept;
    const Electric::ElectricVofSharpMaxwellStress& sharpStress() const;
    const HydroStrainOperator& geometry() const;
    const std::vector<double>& momentumRhsLoad() const;
    const HydroMaxwellLoadReport& report() const;
    void validateFor(const QuadtreeMesh&,const HydroPressureLoad&,const HydroPhaseGeometry&,
        const Electric::ElectricVofMaxwellStress&,LAMP::TimeLayer) const;
    void validateSharpFor(const QuadtreeMesh&,const HydroPressureLoad&,const HydroPhaseGeometry&,
        const Electric::ElectricVofSharpMaxwellStress&,LAMP::TimeLayer) const;
private:
    struct Data;
    static HydroMaxwellLoad assemble(std::shared_ptr<Data>,
        const std::vector<unsigned>& volumePhases,const std::vector<unsigned>& facePhases);
    std::shared_ptr<const Data> data_;
    explicit HydroMaxwellLoad(std::shared_ptr<const Data>);
};
}
