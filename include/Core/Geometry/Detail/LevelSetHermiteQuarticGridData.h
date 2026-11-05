#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Core/Field/VertField.h"
#include "Core/Geometry/Detail/LevelSetHermiteQuarticCell.h"
#include "Core/Geometry/LevelSetCurvature.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

namespace LAMP::Geometry::Detail
{
namespace
{
namespace LevelSetHermiteScalar
{
// Shared data generation for the existing complete equal-level grid port.
// Reuses its validated native labels, original corner phi and existing Q4 fit.
// Each edge midpoint is assembled ONCE from equally nearest admitted Q4 fits
// and supplies both neighboring cells; no positive-coordinate tie bias.
// No root/projection/certificate or production Snapshot switch is implied.
// Native hanging edges and explicit corner sectors require their own retained
// constraints; they are NOT accepted through a uniform background substitute.
class QuarticGridData
{
  public:
    using ExpandedPolynomial = QuarticCell::ExpandedPolynomial;
    using BernsteinControls = std::array<long double, 25>;
    using RationalSegmentCoordinate = LevelSetHermiteScalar::RationalSegmentCoordinate;
    QuarticGridData(const QuadtreeMesh &mesh, const VertField<double> &phi,
                    bool radialAxisSymmetry = false, CartesianGridScope scope = CartesianGridScope::RootOnly)
        : QuarticGridData(mesh, phi, reconstruction(mesh, phi, radialAxisSymmetry, scope))
    {
    }

    // Same retained quartic cells; only their data reconstruction is enriched
    // by recursive extension of the existing Q4 fits. Keep the symmetric-Q4
    // constructor as an explicit baseline; never change an existing snapshot.
    static QuarticGridData fromEnrichedQ4(const QuadtreeMesh& mesh,const VertField<double>& phi,
        bool radialAxisSymmetry=false,CartesianGridScope scope=CartesianGridScope::RootOnly)
    {
        return QuarticGridData(mesh,phi,reconstruction(mesh,phi,radialAxisSymmetry,scope),true);
    }

    const UniformLeafGrid &grid() const noexcept { return grid_; }
    std::size_t cellCount() const noexcept { return cells_.size(); }
    int nx() const noexcept { return grid_.nx(); }
    int ny() const noexcept { return grid_.ny(); }
    Vector2D minCoord() const { return grid_.minCoord(); }
    Vector2D maxCoord() const { return grid_.maxCoord(); }
    Vector2D cellSize() const { return grid_.spacing(); }
    std::size_t vertexCount() const { return grid_.vertices().size(); }
    bool radialAxisSymmetry() const noexcept { return radial_; }
    bool restrictedDomain() const noexcept { return false; }
    bool hasCornerSectors() const noexcept { return false; }
    bool cellActive(int i,int j) const { cell(i,j); return true; }

    // Scalar adapter for the SAME geometry engine as cubic snapshots. These
    // queries never refit or borrow another scalar; native/active patch and
    // corner-sector construction remain unsupported by this data owner.
    BernsteinControls cellBernsteinApproximate(int i,int j) const
    { return cell(i,j).bernsteinControlsApproximate(); }
    Sample evaluateCell(int i,int j,double x,double y) const
    { return evaluateCellChart(i,j,{{x,false},{y,false}}); }
    Sample evaluateCellChart(int i,int j,const CellChartPoint& p) const
    { return cell(i,j).approximateJetAtChart(p); }
    ThirdDerivatives thirdDerivativesCell(int i,int j,double x,double y) const
    { return thirdDerivativesCellChart(i,j,{{x,false},{y,false}}); }
    ThirdDerivatives thirdDerivativesCellChart(int i,int j,const CellChartPoint& p) const
    { return cell(i,j).approximateThirdDerivativesAtChart(p); }
    ThirdDerivatives reconstructedThirdDerivativesCell(int,int,double,double) const
    {
        throw std::runtime_error("quartic data do not expose a separate reconstructed D3 field; use same-cell derivatives");
    }
    ExactExpansion::Interval cellBounds(int i,int j,Vector2D lo={0,0},Vector2D hi={1,1}) const
    { return cell(i,j).cellBounds(lo,hi); }
    ExactExpansion::Interval evaluateEnclosureCell(int i,int j,double x,double y) const
    { return cell(i,j).evaluateEnclosureAtChart({{x,false},{y,false}}); }
    ExactExpansion::Expansion evaluateExactCell(int i,int j,double x,double y) const
    { return cell(i,j).localDerivativeExactAtChart({{x,false},{y,false}}); }
    int evaluateSignCell(int i,int j,double x,double y) const
    { return ExactExpansion::sign(evaluateExactCell(i,j,x,y)); }
    int cellSignExceptVertices(int i,int j,unsigned mask) const
    { return cell(i,j).cellSignExceptVertices(mask); }
    ExpandedPolynomial segmentPolynomialCleared(int i,int j,const RationalSegmentCoordinate& x,
        const RationalSegmentCoordinate& y) const
    { return cell(i,j).homogeneousRestriction(x,y).numerator; }
    ExpandedPolynomial segmentPolynomialExact(int i,int j,Vector2D a,Vector2D b) const
    {
        const auto coordinate=[](double first,double last)
        { return RationalSegmentCoordinate{ExactExpansion::fromDouble(first),ExactExpansion::fromDouble(last),
            ExactExpansion::fromDouble(1.0)}; };
        return segmentPolynomialCleared(i,j,coordinate(a.x,b.x),coordinate(a.y,b.y));
    }
    std::array<double,4> cellNodeDofs(int i,int j,int x,int y) const
    {
        Detail::require((x==0 || x==1) && (y==0 || y==1),"quartic corner indices must be zero or one");
        const auto& data=cell(i,j).localDofs();
        return {{data[5*y+x],data[5*y+x+2],data[5*(y+2)+x],data[5*(y+2)+x+2]}};
    }
    std::array<double,4> nodeDofs(int i,int j) const
    {
        Detail::require(i>=0 && j>=0 && i<=nx() && j<=ny(),"quartic vertex index outside grid");
        const int ci=std::min(i,nx()-1),cj=std::min(j,ny()-1);
        return cellNodeDofs(ci,cj,i-ci,j-cj);
    }
    const QuarticCell &cell(int i, int j) const
    {
        ExactExpansion::validateEnvironment();
        Detail::require(i >= 0 && i < grid_.nx() && j >= 0 && j < grid_.ny(),
                        "quartic grid cell index outside retained layout");
        return cells_[static_cast<std::size_t>(j) * grid_.nx() + i];
    }
    QuarticCell::Dofs cellMappedDofs(int i,int j) const { return cell(i,j).localDofs(); }

  private:
    using Reconstruction = LevelSetCurvature::Reconstruction<QuadtreeMesh>;
    UniformLeafGrid grid_;
    std::vector<QuarticCell> cells_;
    bool radial_ = false;

    static Reconstruction reconstruction(const QuadtreeMesh &mesh, const VertField<double> &phi,
                                          bool radial, CartesianGridScope scope)
    {
        ExactExpansion::validateEnvironment();
        return Reconstruction(mesh, phi, 4, LevelSetCurvature::Basis::TensorProduct, radial, scope);
    }

    static double scaled(double derivative, double spacing)
    {
        const double result = Detail::checked(static_cast<long double>(derivative) * spacing,
                                               "unrepresentable scaled quartic edge derivative");
        Detail::require(derivative == 0.0 || result != 0.0, "quartic edge derivative scaling underflow");
        return result;
    }

    // q4 is a constructor-only borrow of a temporary; all retained data below
    // are owned values. There is no saved reconstruction, mesh or field pointer.
    QuarticGridData(const QuadtreeMesh &mesh, const VertField<double> &phi, const Reconstruction &q4,bool enriched=false)
        : grid_(q4.grid()), radial_(q4.radialAxisSymmetry())
    {
        const auto h = grid_.spacing();
        Detail::require(grid_.nx() <= std::numeric_limits<int>::max()/2 && grid_.ny() <= std::numeric_limits<int>::max()/2,
                        "quartic midpoint index range exceeds supported doubled indices");
        const auto nx = static_cast<std::size_t>(grid_.nx());
        std::vector<Node> corners(grid_.vertices().size());
        for (std::size_t k = 0; k < corners.size(); ++k)
        {
            const auto &vertex = mesh.verts(grid_.vertices()[k]);
            const auto sample = enriched
                ? q4.evaluateHermiteDataJet(2*static_cast<int>(k%(nx+1)),2*static_cast<int>(k/(nx+1))).derivatives
                : q4.evaluateDerivatives(vertex.coord());
            auto &node = corners[k];
            // Preserve the old nodal VALUE and scaled derivative operation
            // order; the named enriched factory intentionally improves the
            // numerical derivatives. QR values never replace corner phi.
            node.value = phi.value[vertex.index()];
            node.x = scaled(sample.dx, h.x);
            node.y = scaled(sample.dy, h.y);
            node.xy = Detail::checked(static_cast<long double>(sample.dxy) * h.x * h.y,
                                       "unrepresentable quartic nodal mixed derivative");
            Detail::require(sample.dxy == 0.0 || node.xy != 0.0, "quartic mixed derivative scaling underflow");
        }
        // Entries are (value, scaled physical normal derivative) in actual
        // validated edge order. No per-cell averaging or independent refit.
        std::vector<std::array<double, 2>> vertical(grid_.vertical().size()), horizontal(grid_.horizontal().size());
        const auto midpoint=[&](int x,int y)
            { return enriched?q4.evaluateHermiteDataJet(x,y):q4.evaluateMidpointJet(x,y); };
        for (std::size_t k = 0; k < vertical.size(); ++k)
        {
            const auto sample = midpoint(2*static_cast<int>(k%(nx+1)), 2*static_cast<int>(k/(nx+1))+1);
            vertical[k] = {{sample.value, scaled(sample.derivatives.dx, h.x)}};
        }
        for (std::size_t k = 0; k < horizontal.size(); ++k)
        {
            const auto sample = midpoint(2*static_cast<int>(k%nx)+1, 2*static_cast<int>(k/nx));
            horizontal[k] = {{sample.value, scaled(sample.derivatives.dy, h.y)}};
        }
        cells_.reserve(grid_.cells().size());
        for (int j = 0; j < grid_.ny(); ++j)
            for (int i = 0; i < grid_.nx(); ++i)
            {
                QuarticCell::Dofs data{};
                for (int y = 0; y < 2; ++y)
                    for (int x = 0; x < 2; ++x)
                    {
                        const auto &node = corners[static_cast<std::size_t>(j + y) * (nx + 1) + i + x];
                        data[5 * y + x] = node.value;
                        data[5 * y + x + 2] = node.x;
                        data[5 * (y + 2) + x] = node.y;
                        data[5 * (y + 2) + x + 2] = node.xy;
                    }
                for (int side = 0; side < 2; ++side)
                {
                    const auto &v = vertical[static_cast<std::size_t>(j) * (nx + 1) + i + side];
                    const auto &hEdge = horizontal[static_cast<std::size_t>(j + side) * nx + i];
                    data[20 + side] = v[0];
                    data[22 + side] = v[1];
                    data[5 * side + 4] = hEdge[0];
                    data[5 * (side + 2) + 4] = hEdge[1];
                }
                data[24] = midpoint(2*i+1,2*j+1).value;
                cells_.emplace_back(data, h.x, h.y);
            }
    }
};
} // namespace LevelSetHermiteScalar
} // namespace
} // namespace LAMP::Geometry::Detail
