#pragma once
#include "Modules/Hydrodynamics/CchMomentumGeometry.h"
#include <array>

namespace Quantum::Hydrodynamics
{
class CchAdvectorBoundary;
class CchAdvectorReconstruction;
struct CchAdvectorCellJet;
struct CchNormalPolynomial;
namespace Detail
{
// Implementation data of CchAdvectorReconstruction, not a second field owner.
struct CchRegularCell
{
    int nativeCell=-1;
    Vector2D lower{0,0}, upper{0,0};
    std::array<int,4> vertices{}, edges{};
    std::vector<int> faces;
    std::array<std::vector<double>,2> supportKnots;
};
struct CchRegularEdge
{
    std::array<int,2> vertices{};
    int normalAxis=0, parent=-1;
    Vector2D lower{0,0}, upper{0,0};
};
struct CchRegularFace
{
    int normalAxis=0, negativePressure=-1, positivePressure=-1;
    Vector2D lower{0,0}, upper{0,0};
    std::array<double,2> supportWidth{};
    long double area=0;
    bool external=false;
    std::array<double,25> meanWeights{};
};
struct CchRegularRow
{
    enum class Kind { Vertex, Edge, Center };
    Kind kind=Kind::Vertex;
    int index=-1, component=0, rank=0, rings=0;
    bool zeroNormal=false;
    Vector2D point{0,0}, scale{0,0};
    std::vector<int> columns, support;
    std::array<std::vector<double>,4> weights; // value, dx, dy, dxy
    double coefficientNorm=0, deficit=0, nullResidual=0;
    int priorityDimension=0, priorityRank=0;
    double priorityDeficit=0;
    std::array<double,6> nullPolynomial{}; // unit certificate when deficient
};
struct CchRegularPreparationReport
{
    std::size_t svdCount=0, weights=0, maxSupport=0, rankDeficient=0;
    double maxCoefficientNorm=0;
    std::size_t normalPorts=0, normalPortWeights=0;
};
// Final normal value and tangential derivative at a shared master vertex.
// Observations are physical native face means on the same coordinate line.
struct CchRegularNormalPort
{
    int vertex=-1, component=0, observedModes=0;
    std::vector<int> support;
    std::array<std::vector<double>,2> weights;
};
struct CchRegularCoefficients
{
    // Tensor index 5*y+x, axes [value0,value1,derivative0,derivative1,mid].
    std::vector<std::array<std::array<double,25>,2>> cellDofs;
    std::vector<double> faceCorrections;
private:
    std::shared_ptr<const void> geometry_; // Same immutable preparation; no velocity/clock identity.
    friend class CchRegularAdvectorGeometry;
};
struct CchRegularScalarWeight { int velocityDof=-1; double weight=0; };
struct CchRegularJetWeight
{
    int velocityDof=-1;
    // value x/y, dx x, dy x, dx y, dy y of the same source.
    std::array<double,6> jet{};
};
struct CchRegularCellStencil
{
    int pressure=-1;
    std::array<std::array<std::vector<CchRegularScalarWeight>,25>,2> cellDofs;
    std::vector<std::vector<CchRegularScalarWeight>> faceCorrections;
private:
    std::shared_ptr<const void> geometry_;
    friend class CchRegularAdvectorGeometry;
};
class CchRegularAdvectorGeometry
{
public:
    const CchMomentumGeometry& momentumGeometry() const noexcept;
    const std::vector<CchRegularCell>& cells() const noexcept;
    const std::vector<CchRegularEdge>& edges() const noexcept;
    const std::vector<CchRegularFace>& faces() const noexcept;
    const std::vector<CchRegularRow>& rows() const noexcept;
    const std::vector<CchRegularNormalPort>& normalPorts() const noexcept;
    const CchRegularPreparationReport& report() const noexcept;
    bool sameIdentity(const CchRegularAdvectorGeometry&) const noexcept;
    void validateFor(const QuadtreeMesh&, const CchMomentumGeometry&) const;
    // Only sparse weight application and shared Q4 algebra; never a fit/SVD.
    CchRegularCoefficients apply(const QuadtreeMesh&, const CchNativeValues&,
        const CchAdvectorBoundary* =nullptr) const;
    // Local symbolic application of the existing prepared source. No global
    // unit-velocity probes, new field, or runtime moment solve.
    CchRegularCellStencil cellStencil(const QuadtreeMesh&,int pressure,
        const CchAdvectorBoundary* =nullptr) const;
    std::vector<CchRegularJetWeight> cellJetWeights(const CchRegularCellStencil&,
        Vector2D point) const;
    // Analytic queries of these same coefficients, not a separately published
    // field. Interior support knots select the positive-side closed piece.
    CchAdvectorCellJet cellJet(const CchRegularCoefficients&,int pressure,Vector2D point) const;
    // Canonical axis-aligned normal trace inside one cell/support piece.
    CchNormalPolynomial normalPolynomial(const CchRegularCoefficients&,int pressure,
        int normalAxis,Vector2D begin,Vector2D end) const;
private:
    void validateCoefficients(const CchRegularCoefficients&,int pressure) const;
    static CchRegularAdvectorGeometry prepare(const QuadtreeMesh&, CchMomentumGeometry);
    struct Data;
    explicit CchRegularAdvectorGeometry(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
    friend class Quantum::Hydrodynamics::CchAdvectorReconstruction;
};
}
}
