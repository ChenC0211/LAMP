#pragma once

#include <optional>

#include "Modules/Hydrodynamics/HydroPhysicalDivergence.h"
#include <Eigen/SparseCore>
#include <memory>
#include <vector>

namespace LAMP::Mesh { class NativeLeafOverlap; }

namespace LAMP::Hydrodynamics
{
class HydroMomentumGeometry;
class HydroPhaseViscosity;
// Current owned raw preparation only; never a duration retained in strain Data.
struct HydroNoSlipPreparationReport
{
    bool available=false,rawCacheHit=false,coldBuilt=false;
    bool parallelUsed=false,parallelFallback=false,timingAvailable=false;
    int maximumWorkers=0;
    double coldSeconds=0;
    std::size_t fitRequests=0,fitHits=0,qrFactorizations=0,svdFactorizations=0;
    std::size_t bankPeakBytes=0,bankLimitBytes=32u*1024u*1024u;
};
// A volume quadrature sample, not a material/phase classifier. Physical weight
// includes 2*pi*r in RZ. H07 must supply cut samples from the unique geometry.
struct HydroStrainPoint
{
    int cellDof = -1;
    Vector2D point{0,0};
    double physicalWeight = 0;
};

// Actual shared native-face quadrature for a virtual velocity, not a force.
// The normal trace is the native mean; the tangent is one shared equal mean
// of the adjacent original P2 values (one-sided at a physical boundary).
struct HydroStressFacePoint
{
    int velocityDof=-1,nativeFace=-1,negativeCellDof=-1,positiveCellDof=-1;
    Vector2D point{0,0};
    double physicalWeight=0; // Includes 2*pi*r once in RZ; true axis is zero.
};
struct HydroStressCompatibilityReport
{
    std::vector<double> cellQuadratureWeightSums,cellQuadratureWeightDefects;
    double maximumRelativeQuadratureWeightDefect=0;
    double maximumRelativeFaceQuadratureWeightDefect=0;
};

// Unique target cells/face rows certified for local same-End AMR geometry reuse.
// No material values, physical step counts or old owner handles are included.
struct HydroStressGeometryReuseReport
{
    bool requested=false,enabled=false;
    std::size_t targetCells=0,reusedRawCells=0,reusedStressCells=0,recomputedFitCells=0;
    std::size_t reusedFaces=0,rebuiltFaces=0,stencilRejectedCells=0;
    std::size_t quadratureRejectedCells=0,columnRejectedCells=0;
    // Original direct/literal counters: shared face, cell integral and Ehat only.
    // Raw/volume have separate counters below; scratch is the total construction peak.
    std::size_t directReboundRows=0,literalReboundRows=0,csrRebindScratchBytes=0;
    // Separate Task49 actual scalar source-copy rows; construction snapshot.
    std::size_t rawDirectReboundRows=0,volumeDirectReboundRows=0;
};
// Derived local P2 coefficients are eliminated into this native-velocity map.
// Rows per point: XY (xx,yy,sqrt(2)xy), RZ (rr,zz,hoop,sqrt(2)rz).
class HydroStrainOperator
{
public:
    using Matrix = Eigen::SparseMatrix<double,Eigen::RowMajor,Eigen::Index>;
    static std::vector<HydroStrainPoint> volumeQuadrature(const QuadtreeMesh& mesh,
        const HydroPhysicalDivergence& metrics);
    static HydroStrainOperator build(const QuadtreeMesh& mesh, HydroPhysicalDivergence metrics,
        const std::vector<HydroStrainPoint>& points);
    // Explicit stationary no-slip on actual positive-area boundary DOFs.
    // P2 fits consume physical native face means, and homogeneous full traces
    // on one coplanar wall per adjacent cell. Unsupported multiwall corners
    // or insufficient support reject. No boundary penalty/viscosity floor.
    static HydroStrainOperator buildStationaryNoSlip(const QuadtreeMesh&,
        HydroPhysicalDivergence,const std::vector<HydroStrainPoint>&,
        const std::vector<int>& wallVelocityDofs);
    const std::vector<int>& stationaryNoSlipVelocityDofs() const noexcept;
    // Geometry-only conservative stress virtual strain. The original build
    // and viscosity keep their original raw P2 strain arithmetic.
    static HydroStrainOperator buildStressCompatible(const QuadtreeMesh&,
        HydroPhysicalDivergence metrics);
    // Same native P2/Ehat operator on an explicitly partitioned volume/face
    // quadrature. Each actual cell and face must have its full physical weight.
    static HydroStrainOperator buildStressCompatible(const QuadtreeMesh&,
        HydroPhysicalDivergence metrics,const std::vector<HydroStrainPoint>&,
        const std::vector<HydroStressFacePoint>&);
    // Fresh target metrics/CSR, copying only exactly certified local scalar
    // blocks of the builtin canonical stress factory. Unsupported factories
    // keep the original full fresh path; both live meshes/overlap are checked.
    static HydroStrainOperator buildStressCompatibleRemapped(const QuadtreeMesh& sourceMesh,
        const HydroStrainOperator& sourceGeometry,const QuadtreeMesh& targetMesh,
        HydroPhysicalDivergence targetMetrics,const LAMP::Mesh::NativeLeafOverlap&);
    const HydroStressGeometryReuseReport& geometryReuseReport() const noexcept;
    bool isStressCompatible() const noexcept;
    // P2 velocity maps for the actual RT0 tensor Green pairing. Full also
    // satisfies this capability; Green-only has no strain/moment operator.
    static HydroStrainOperator buildRt0Green(const QuadtreeMesh&,HydroPhysicalDivergence);
    static HydroStrainOperator buildRt0GreenRemapped(const QuadtreeMesh& sourceMesh,
        const HydroStrainOperator& sourceGeometry,const QuadtreeMesh& targetMesh,
        HydroPhysicalDivergence targetMetrics,const LAMP::Mesh::NativeLeafOverlap&);
    bool isRt0GreenCompatible() const noexcept;
    // Only the builtin complete, original-order volumeQuadrature stress factory
    // retains its actual pre-correction raw CSR. Partitioned/sharp factories do not.
    bool hasCanonicalRawStrain() const noexcept;
    HydroStrainOperator canonicalRawStrain() const;
    const std::vector<HydroStressFacePoint>& stressFacePoints() const;
    const Matrix& stressFaceVelocityMatrix() const; // 2 rows/face point: x/r,y/z.
    const Matrix& stressVolumeVelocityMatrix() const; // 2 rows/volume point, raw P2.
    const Matrix& stressComponentIntegralMatrix() const; // cellDof*components()+k.
    const HydroStressCompatibilityReport& stressCompatibilityReport() const;
    const HydroPhysicalDivergence& metrics() const noexcept;
    const HydroDofLayout& layout() const noexcept;
    const std::vector<HydroStrainPoint>& points() const noexcept;
    int components() const noexcept;
    const Matrix& matrix() const; // Green-only rejects unrequested strain.
    std::vector<double> apply(const QuadtreeMesh& mesh, const HydroNativeValues& values) const;
private:
    // Only the explicit zero bilinear form creates an inactive empty strain
    // map; the ordinary sampled-strain factory still rejects empty points.
    friend class HydroEnergyViscosity;
    struct Data;
    static HydroStrainOperator buildStationaryNoSlipImpl(const QuadtreeMesh&,
        HydroPhysicalDivergence,const std::vector<HydroStrainPoint>&,
        const std::vector<int>&,HydroNoSlipPreparationReport*);
    static HydroStrainOperator buildImpl(const QuadtreeMesh&, HydroPhysicalDivergence,
        const std::vector<HydroStrainPoint>&,bool stressCompatible,
        const std::vector<HydroStressFacePoint>* partitionedFaces=nullptr,
        bool retainCanonicalRaw=false,const QuadtreeMesh* remapSourceMesh=nullptr,
        const HydroStrainOperator* remapSource=nullptr,
        const LAMP::Mesh::NativeLeafOverlap* remapOverlaps=nullptr,bool rt0GreenOnly=false);
    explicit HydroStrainOperator(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};

// Actual K-product work of one build, never an estimate of unchanged cells.
// reason: 0 qualified row work; 1 first Eigen/no row receipt; 2 unsupported FP;
// 3 optional memory scope; 4 concurrent slot busy; 5 optional failure.
struct HydroViscosityProductReport
{
    bool requested=false,eligible=false,sourceAvailable=false,sourceFpReceipt=false;
    bool currentFpQualified=false,localUsed=false,fullAssembly=false,memoPublished=false;
    int fallbackReason=0;
    std::size_t weightRows=0,changedWeightRows=0,outputRows=0,recomputedRows=0,reusedRows=0;
    std::size_t symbolicStoredValues=0,retainedBytes=0,peakOptionalBytes=0;
    std::size_t byteLimit=128u*1024u*1024u;
    // Passive current-build timings. Zero symbolic workers means no new
    // structure work; it is independent of timing availability.
    bool timingAvailable=false,symbolicBuilt=false;
    int symbolicWorkers=0;
    double columnRoutingSeconds=0,symbolicGraphSeconds=0,numericAssemblySeconds=0;
    HydroNoSlipPreparationReport noSlip;
};

// K_mu = E^T diag(2*mu_q*dV_q) E; momentum uses -K_mu*u.
// No minimum viscosity, boundary penalty, phase averaging or new primary DOF.
class HydroEnergyViscosity
{
public:
    using Matrix = HydroStrainOperator::Matrix;
    // Exact inactive form retaining the actual native metrics and axis law.
    // Strain has zero rows and no samples; K is velocityCount x velocityCount.
    static HydroEnergyViscosity zero(const QuadtreeMesh& mesh, HydroPhysicalDivergence metrics);
    static HydroEnergyViscosity build(HydroStrainOperator strain, const std::vector<double>& viscosity);
    const HydroStrainOperator& strain() const noexcept;
    const std::vector<double>& viscosity() const noexcept;
    const std::vector<double>& rowWeights() const noexcept;
    const Matrix& matrix() const noexcept;
    const HydroViscosityProductReport& productReport() const noexcept;
    // Certified absolute row-sum bounds for diag(s)*(L(K)-G)*diag(s),
    // L(K): symmetric completion of the actual stored lower triangle;
    // G: exact-real Gram of stored strain and actual binary64 rowWeights.
    // Source build and caller must be strict RNE/masked/gradual binary64.
    // Unsupported/invalid/overlarge input returns no certificate. The query
    // preserves the complete caller fenv/CSR/errno and never changes K.
    // Coarse scalar upper bound valid for every scaled row when 0<s_i<=maximumScale.
    // Same exact Lower/Gram and source/current FP contract as the detailed query.
    std::optional<double> coarseScaledGramLowerErrorBound(double maximumScale) const;
    std::optional<std::vector<double>> scaledGramLowerErrorBounds(
        const std::vector<double>& positiveScale) const;
    std::vector<double> apply(const QuadtreeMesh& mesh, const HydroNativeValues& values) const;
    double dissipation(const QuadtreeMesh& mesh, const HydroNativeValues& values) const;
private:
    friend class HydroPhaseViscosity;
    // Internal optional same-Momentum raw geometry; current viscosity/K stay fresh.
    static HydroEnergyViscosity buildStationaryNoSlipOwned(const QuadtreeMesh&,
        const HydroMomentumGeometry&,HydroPhysicalDivergence,const std::vector<HydroStrainPoint>&,
        const std::vector<int>& stationaryNoSlipVelocityDofs,const std::vector<double>& viscosity);
    static HydroEnergyViscosity buildImpl(HydroStrainOperator,const std::vector<double>&,
        const HydroNoSlipPreparationReport*);
    struct Data;
    explicit HydroEnergyViscosity(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
}
