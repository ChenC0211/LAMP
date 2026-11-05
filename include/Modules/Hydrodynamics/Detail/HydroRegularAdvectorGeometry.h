#pragma once
#include "Modules/Hydrodynamics/HydroMomentumGeometry.h"
#include <array>
#include "Core/Geometry/LevelSetBoundaryTrace.h"

namespace LAMP::Hydrodynamics
{
class HydroAdvectorBoundary;
class HydroAdvectorReconstruction;
struct HydroAdvectorCellJet;
struct HydroAdvectorCellJetPair;
struct HydroNormalPolynomial;
struct HydroAdvectorPiece;
struct HydroAdvectorCurveIntegral;
struct HydroAdvectorCurveIntegralPair;
namespace Detail
{
// Implementation data of HydroAdvectorReconstruction, not a second field owner.
struct HydroRegularCell
{
    int nativeCell=-1;
    Vector2D lower{0,0}, upper{0,0};
    std::array<int,4> vertices{}, edges{};
    std::vector<int> faces;
    std::array<std::vector<double>,2> supportKnots;
};
struct HydroRegularEdge
{
    std::array<int,2> vertices{};
    int normalAxis=0, parent=-1;
    Vector2D lower{0,0}, upper{0,0};
};
struct HydroRegularFace
{
    int normalAxis=0, negativePressure=-1, positivePressure=-1;
    Vector2D lower{0,0}, upper{0,0};
    std::array<double,2> supportWidth{};
    long double area=0;
    bool external=false;
    std::array<double,25> meanWeights{};
};
struct HydroRegularRow
{
    enum class Kind { Vertex, Edge, Center };
    Kind kind=Kind::Vertex;
    int index=-1, component=0, rank=0, rings=0, degree=2;
    int p2Rank=0, p2ColumnCount=0; // Independent old-degree observation eligibility.
    double p2Deficit=0; // Actual inverse preserves the original identifiable P2 representative.
    bool zeroNormal=false;
    Vector2D point{0,0}, scale{0,0};
    std::vector<int> columns, support;
    std::array<std::vector<double>,4> weights; // value, dx, dy, dxy
    double coefficientNorm=0, deficit=0, nullResidual=0;
    int priorityDimension=0, priorityRank=0;
    double priorityDeficit=0;
    std::array<double,10> nullPolynomial{}; // Actual selected-degree unit certificate when deficient
};
struct HydroRegularPreparationReport
{
    std::size_t svdCount=0, weights=0, maxSupport=0, rankDeficient=0;
    double maxCoefficientNorm=0;
    std::size_t normalPorts=0, normalPortWeights=0;
    std::size_t virtualCells=0,virtualFaces=0,prolongationWeights=0;
    // Full normalized local operators, never reused state/phase/velocity data.
    // svdCount above still counts every actual main/port/lower solve, including
    // unsuccessful template attempts followed by the original literal path.
    std::size_t templateLookups=0,templateRowQueries=0,templateRowHits=0,templateRowBuilds=0;
    std::size_t templateRoutingHits=0,templatePortQueries=0,templatePortHits=0,templatePortBuilds=0;
    std::size_t templateFallbacks=0,templatePublishedEntries=0,templateEntryCount=0,templateRetainedBytes=0;
    std::size_t templateRowEntries=0,templateRoutingEntries=0;
    std::size_t templateEntryLimit=0,templateByteLimit=0; // accounted payload plus conservative entry overhead
};
// Final normal value and tangential derivative at a shared master vertex.
// Observations are physical native face means on the same coordinate line.
struct HydroRegularNormalPort
{
    int vertex=-1, component=0, observedModes=0;
    std::vector<int> support;
    std::array<std::vector<double>,2> weights;
};
struct HydroRegularCoefficients
{
    // General cells: Q4 Hermite data, tensor index 5*y+x.
    // Compatible XY: sixteen quarter-cell power patches per native uniform
    // cell, or per derived virtual cell of the same native U on one-level AMR.
    // The accepted mesh/layout/field are unchanged. Uniform patches are
    // ordered 16*pressure+4*quarterY+quarterX; the same 25-slot tensor storage.
    std::vector<std::array<std::array<double,25>,2>> cellDofs;
    std::vector<double> faceCorrections;
    // Same applied native means retained for exact cardinal face restrictions.
    // These are coefficient auxiliaries, with no separate field/provenance owner.
    std::vector<double> compatibleNativeMeans;
private:
    std::shared_ptr<const void> geometry_; // Same immutable preparation; no velocity/clock identity.
    friend class HydroRegularAdvectorGeometry;
};
// Same endpoint field, in this actual support piece's unit-square chart.
// Prepared once from cell DOFs and the active native face corrections.
struct HydroRegularPiecePolynomial
{
    int pressure=-1;
    Vector2D lower{0,0},upper{0,0};
    std::array<std::array<long double,25>,2> powers{}; // 5*y+x
private:
    std::shared_ptr<const void> geometry_;
    friend class HydroRegularAdvectorGeometry;
};
struct HydroRegularScalarWeight { int velocityDof=-1; double weight=0; };
struct HydroRegularJetWeight
{
    int velocityDof=-1;
    // value x/y, dx x, dy x, dx y, dy y of the same source.
    std::array<double,6> jet{};
};
struct HydroRegularCellStencil
{
    int pressure=-1;
    std::array<std::array<std::vector<HydroRegularScalarWeight>,25>,2> cellDofs;
    std::vector<std::vector<HydroRegularScalarWeight>> faceCorrections;
private:
    std::shared_ptr<const void> geometry_;
    std::vector<int> excludedVelocityDofs_; // local prescribed means, compatible source only
    friend class HydroRegularAdvectorGeometry;
};
class HydroRegularAdvectorGeometry
{
public:
    const HydroMomentumGeometry& momentumGeometry() const noexcept;
    const std::vector<HydroRegularCell>& cells() const noexcept;
    const std::vector<HydroRegularEdge>& edges() const noexcept;
    const std::vector<HydroRegularFace>& faces() const noexcept;
    const std::vector<HydroRegularRow>& rows() const noexcept;
    const std::vector<HydroRegularNormalPort>& normalPorts() const noexcept;
    const HydroRegularPreparationReport& report() const noexcept;
    bool usesCompatibleTensor() const noexcept; // native uniform or eligible one-level XY AMR compact L/E tensor
    bool usesVirtualCompatibleTensor() const noexcept; // derived sparse prolongation of the original native U
    // Pure topology scatter plans successfully prepared on this immutable
    // regular owner. Counts tensor-cell/mask plans, never velocity evaluations.
    std::size_t cachedVirtualJetPlans() const noexcept;
    bool sameIdentity(const HydroRegularAdvectorGeometry&) const noexcept;
    void validateFor(const QuadtreeMesh&, const HydroMomentumGeometry&) const;
    // Only sparse weight application and shared Q4 algebra; never a fit/SVD.
    HydroRegularCoefficients apply(const QuadtreeMesh&, const HydroNativeValues&,
        const HydroAdvectorBoundary* =nullptr) const;
    // Local symbolic application of the existing prepared source. No global
    // unit-velocity probes, new field, or runtime moment solve.
    HydroRegularCellStencil cellStencil(const QuadtreeMesh&,int pressure,
        const HydroAdvectorBoundary* =nullptr) const;
    // Same private geometry identity and symbolic rows, with mesh/boundary
    // validation performed once for this requested compatible-cell batch.
    std::vector<HydroRegularCellStencil> cellStencils(const QuadtreeMesh&,
        const std::vector<int>& pressures,const HydroAdvectorBoundary* =nullptr) const;
    std::vector<HydroRegularJetWeight> cellJetWeights(const HydroRegularCellStencil&,
        Vector2D point) const;
    // Analytic queries of these same coefficients, not a separately published
    // field. Interior support knots select the positive-side closed piece.
    HydroAdvectorCellJet cellJet(const HydroRegularCoefficients&,int pressure,Vector2D point) const;
    // One known actual piece: share basis/bubble geometry while contracting
    // each endpoint's original coefficients in its original arithmetic order.
    HydroAdvectorCellJetPair pieceJetPair(const HydroRegularCoefficients& previous,
        const HydroRegularCoefficients& arrival,const HydroAdvectorPiece&,Vector2D point) const;
    // Algebraic conversion of the actual endpoint field; never a fit.
    std::vector<HydroRegularPiecePolynomial> piecePolynomials(const HydroRegularCoefficients&,
        const std::vector<HydroAdvectorPiece>&) const;
    HydroAdvectorCurveIntegral curveIntegral(const HydroRegularPiecePolynomial&,
        Vector2D r0,Vector2D r1,Vector2D r2,
        LevelSetBoundaryTrace::Jacobian P0,LevelSetBoundaryTrace::Jacobian P1,
        LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const;
    // Common geometric monomial moments of one actual quadratic span, then
    // linear contraction with the two retained endpoint polynomials.
    HydroAdvectorCurveIntegralPair curveIntegralPair(const HydroRegularPiecePolynomial& previous,
        const HydroRegularPiecePolynomial& arrival,
        Vector2D r0,Vector2D r1,Vector2D r2,
        LevelSetBoundaryTrace::Jacobian P0,LevelSetBoundaryTrace::Jacobian P1,
        LevelSetBoundaryTrace::Jacobian P2,double lo,double hi) const;
    // Canonical axis-aligned normal trace inside one cell/support piece.
    HydroNormalPolynomial normalPolynomial(const HydroRegularCoefficients&,int pressure,
        int normalAxis,Vector2D begin,Vector2D end) const;
private:
    // Only a same-owner synchronous Motion preparation that has already
    // validated mesh/source/values/walls may use this local derivation.
    HydroRegularCellStencil cellStencilPrepared(int pressure,const HydroAdvectorBoundary*) const;
    HydroRegularCellStencil compatibleCellStencil(int pressure,const std::vector<int>* boundaryDofs) const;
    void validateCoefficients(const HydroRegularCoefficients&,int pressure) const;
    // The sole Reconstruction caller has validated this live mesh and actual
    // boundary owner in the same callback-free preparation.
    static HydroRegularAdvectorGeometry prepare(const QuadtreeMesh&, HydroMomentumGeometry,
        const HydroAdvectorBoundary* =nullptr);
    static std::shared_ptr<const HydroRegularRemapPreparation> makeRemapPreparation(
        const QuadtreeMesh& sourceMesh,const HydroRegularAdvectorGeometry& source,
        const QuadtreeMesh& targetMesh,const Mesh::NativeLeafOverlap& actualOverlaps);
    friend class LAMP::Hydrodynamics::HydroMomentumGeometry;
    struct Data;
    explicit HydroRegularAdvectorGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class LAMP::Hydrodynamics::HydroAdvectorReconstruction;
};
}
}
