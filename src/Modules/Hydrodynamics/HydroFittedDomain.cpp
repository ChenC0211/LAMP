#include "Modules/Hydrodynamics/HydroFittedDomain.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Modules/Hydrodynamics/HydroNativePressureStencil.h"


namespace
{
using Type = HydroFittedDomain::BoundaryType;
using BoundaryData = HydroFittedDomain::Boundary;

void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[HydroFittedDomain] ") + message);
}

bool near(double a, double b, double scale)
{
    return std::isfinite(a) && std::isfinite(b) && std::isfinite(scale) &&
        std::abs(a - b) <= 256.0 * std::numeric_limits<double>::epsilon() * scale;
}

void validateBoundary(const BoundaryData& boundary, bool interior)
{
    require(boundary.id >= -1 && std::isfinite(boundary.pressure) &&
            std::isfinite(boundary.incomingDensity) && std::isfinite(boundary.incomingViscosity),
            "invalid boundary ID or non-finite boundary data");
    if (interior)
    {
        require(boundary.type == Type::Interior && boundary.id == -1 && boundary.pressure == 0.0 &&
                boundary.incomingDensity == 0.0 && boundary.incomingViscosity == 0.0,
                "active-active segment requires default Interior data");
        return;
    }
    if (boundary.type == Type::PressureReservoir)
        require(boundary.incomingDensity > 0.0 && boundary.incomingViscosity > 0.0,
                "pressure reservoir requires explicitly positive incoming liquid density/viscosity");
    else
    {
        require(boundary.type == Type::PressureOutlet || boundary.type == Type::NoSlipWall ||
                boundary.type == Type::Symmetry, "every fluid boundary must have an explicit supported condition");
        require(boundary.incomingDensity == 0.0 && boundary.incomingViscosity == 0.0,
                "only pressure reservoirs accept incoming material data");
        if (boundary.type != Type::PressureOutlet)
            require(boundary.pressure == 0.0, "wall/symmetry pressure field must be unused zero");
    }
}

// The virtual-child gradient has a mathematical zero row sum. Scaling its
// coefficients by the dual volume and merging triplets can lose that identity
// in floating point, so C*1 activates spurious phase-plateau loads. Restore it
// in the FINAL matrix consumed by pressure, capillarity and electric loads,
// not by clipping any evaluated force or changing the interpolation stencil.
// The last stored coefficient closes the sum in sparse column-product order.
// This guarantees C*1=0 on interior rows; arbitrary constant amplitudes may
// still incur rounding in the subsequent coefficient/value multiplications.
void preserveInteriorConstantMode(Eigen::SparseMatrix<double>& matrix,
    const std::vector<HydroFittedDomain::NativeFace>& faces, const std::vector<int>& nativeByRow)
{
    std::vector<int> lastColumn(matrix.rows(), -1);
    std::vector<double> prefix(matrix.rows(), 0.0), magnitude(matrix.rows(), 0.0);
    const auto interior = [&](Eigen::Index row) {
        return faces[nativeByRow[row]].boundary.type == Type::Interior;
    };
    for (int column = 0; column < matrix.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column); entry; ++entry)
        {
            require(std::isfinite(entry.value()), "non-finite pressure gradient coefficient");
            if (!interior(entry.row())) continue; // Physical pressure traction is unchanged.
            if (entry.value() != 0.0) lastColumn[entry.row()] = column;
            magnitude[entry.row()] += std::abs(entry.value());
        }
    for (Eigen::Index row = 0; row < matrix.rows(); ++row)
        if (interior(row))
            require(lastColumn[row] >= 0 && std::isfinite(magnitude[row]) && magnitude[row] > 0,
                "empty or non-finite interior pressure gradient row");
    for (int column = 0; column < matrix.outerSize(); ++column)
        for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column); entry; ++entry)
        {
            const auto row = entry.row();
            if (!interior(row)) continue;
            if (column == lastColumn[row])
            {
                const double replacement = -prefix[row];
                // Reject a broken stencil rather than hiding a non-roundoff
                // row defect. No new column, mask or physical tolerance.
                require(std::isfinite(replacement) &&
                    std::abs(replacement - entry.value()) <=
                        64 * std::numeric_limits<double>::epsilon() * magnitude[row],
                    "native pressure gradient does not preserve constants within roundoff");
                entry.valueRef() = replacement;
            }
            prefix[row] += entry.value();
        }
    for (Eigen::Index row = 0; row < matrix.rows(); ++row)
        if (interior(row)) require(prefix[row] == 0.0, "interior pressure gradient row sum is not zero");
}
} // namespace

HydroFittedDomain HydroFittedDomain::build(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<Boundary>& bySegment, bool axisymmetric)
{
    HydroFittedDomain result;
    require(bySegment.size() == graph.segments().size(), "graph offsets/boundary count mismatch");
    result.geometry_ = NativeLeafDomain::build(mesh, graph, axisymmetric);
    const int nv = mesh.vEdgeIndexNum();
    const auto& segments = graph.segments();
    result.faces_.resize(result.geometry_.nativeFaceCount());
    result.velocityRowByNative_.assign(result.faces_.size(), -1);
    for (int row = 0; row < result.geometry_.nativeFaceCount(); ++row)
        static_cast<NativeLeafFace&>(result.faces_[row]) = result.geometry_.face(mesh, row);
    for (std::size_t sid = 0; sid < segments.size(); ++sid)
    {
        const auto& segment = segments[sid];
        const auto& edge = segment.isHorizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
        const int row = edge.index() + (segment.isHorizontal ? nv : 0);
        auto& face = result.faces_[row];
        const bool interior = face.negativeCell >= 0 && face.positiveCell >= 0;
        validateBoundary(bySegment[sid], interior);
        face.boundary = bySegment[sid]; // Interior entries are both exactly default.
        if (!interior)
        {
            if (axisymmetric && face.center.x == 0.0)
                require(!face.horizontal && mesh.minCoord.x == 0.0 && segment.L == -1 &&
                        segment.side == FaceSide::West && face.boundary.type == Type::Symmetry,
                        "the true zero-area radial axis must be a Symmetry boundary");
        }
    }

    const int activeCount = static_cast<int>(result.cellIndexByRow().size());
    std::vector<std::vector<int>> neighbors(activeCount);
    std::vector<unsigned char> anchored(activeCount, 0), visited(activeCount, 0);
    std::vector<Eigen::Triplet<double>> divergenceEntries;
    std::vector<Eigen::Triplet<double>> pressureEntries;
    std::map<std::pair<int,int>,supCoeff> pressureCenters;
    const auto centers = [&](const Cell& cell,int child) -> const supCoeff& {
        const auto key=std::make_pair(cell.index(),child);
        auto found=pressureCenters.find(key);
        if (found==pressureCenters.end())
            found=pressureCenters.emplace(key,HydroNativePressureStencil::childCenter(mesh,cell,child)).first;
        for (int index:found->second.ind)
            require(result.cellRowByIndex()[index]>=0, "coarse/fine pressure interpolation leaves the active fluid");
        return found->second;
    };
    divergenceEntries.reserve(segments.size());
    for (auto& face : result.faces_)
    {
        const bool negative = face.negativeCell >= 0, positive = face.positiveCell >= 0;
        const bool interior = negative && positive;
        if (!negative && !positive) continue; // Solid-only storage, never an unknown.
        const bool pressure = face.boundary.type == Type::PressureReservoir || face.boundary.type == Type::PressureOutlet;
        const int owner = negative ? face.negativeCell : face.positiveCell;
        const auto& ownerCell = mesh.cells(mesh.cellIndexToLabel[owner]);
        auto low = negative ? mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).centerCoord() : face.center;
        auto high = positive ? mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).centerCoord() : face.center;
        bool coarseFine=false;
        if (interior)
        {
            const auto& a=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]);
            const auto& b=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]);
            coarseFine=a.level()!=b.level();
            if (coarseFine)
            {
                // The old momentum/pressure stencil uses a virtual child on
                // the coarse side, so its dual is fine-width, not the distance
                // between tangentially offset real coarse/fine centers.
                if (a.level()<b.level())
                {
                    low=face.center;
                    if (face.horizontal) low.y-=.5*face.normalSpacing;
                    else low.x-=.5*face.normalSpacing;
                }
                else
                {
                    high=face.center;
                    if (face.horizontal) high.y+=.5*face.normalSpacing;
                    else high.x+=.5*face.normalSpacing;
                }
            }
        }
        else require(near(face.tangentSpacing,face.horizontal?ownerCell.size().x:ownerCell.size().y,
            face.tangentSpacing), "boundary half-dual requires an unsplit fluid-owner side");
        const double distance = face.horizontal ? high.y - low.y : high.x - low.x;
        const double normalCenter = face.horizontal ? ownerCell.centerCoord().y : ownerCell.centerCoord().x;
        require(distance > 0.0 && near(distance, (interior ? 1.0 : 0.5) * face.normalSpacing,
                std::max(std::abs(normalCenter), face.normalSpacing)), "invalid actual normal dual extent");
        const auto& edge = face.horizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
        if (axisymmetric && !face.horizontal)
        {
            // Exact real integral between actual centers (or owner and wall).
            // Factored difference of squares also handles an axis half dual.
            face.dualVolume = 0.5 * CartesianMeasure::twoPi * (high.x + low.x) * distance * edge.length();
        }
        else
            face.dualVolume = face.area * distance;
        require(std::isfinite(face.dualVolume) && face.dualVolume > 0.0 &&
                std::isfinite(face.area) && face.area >= 0.0, "unrepresentable physical face/dual measure");
        if (face.area == 0.0)
            require(axisymmetric && !face.horizontal && face.center.x == 0.0 &&
                    !interior && face.boundary.type == Type::Symmetry, "unexpected zero-area physical face");
        if (!interior && !pressure) continue; // u_n=0 eliminated, never pinned.
        require(face.area > 0.0, "free velocity requires positive physical area");
        face.velocityRow = static_cast<int>(result.nativeByVelocityRow_.size());
        result.nativeByVelocityRow_.push_back(face.nativeRow);
        result.velocityRowByNative_[face.nativeRow] = face.velocityRow;
        const int nrow = negative ? result.cellRowByIndex()[face.negativeCell] : -1;
        const int prow = positive ? result.cellRowByIndex()[face.positiveCell] : -1;
        if (negative) divergenceEntries.emplace_back(nrow, face.velocityRow, face.area);
        if (positive) divergenceEntries.emplace_back(prow, face.velocityRow, -face.area);
        if (!result.uniformLeaves())
        {
            if (coarseFine)
                HydroNativePressureStencil::visitGradient(mesh,edge,face.horizontal,centers,
                    [&](int index,double gradient) {
                        require(result.cellRowByIndex()[index]>=0, "pressure coefficient references solid-only storage");
                        pressureEntries.emplace_back(face.velocityRow,result.cellRowByIndex()[index],face.dualVolume*gradient);
                    });
            else
            {
                if (negative) pressureEntries.emplace_back(face.velocityRow,nrow,-face.area);
                if (positive) pressureEntries.emplace_back(face.velocityRow,prow,face.area);
            }
        }
        if (interior)
        {
            neighbors[nrow].push_back(prow); neighbors[prow].push_back(nrow);
        }
        else anchored[negative ? nrow : prow] = 1;
    }
    for (int seed = 0; seed < activeCount; ++seed)
    {
        if (visited[seed]) continue;
        std::queue<int> pending;
        pending.push(seed); visited[seed] = 1;
        bool hasAnchor = false;
        while (!pending.empty())
        {
            const int row = pending.front(); pending.pop();
            hasAnchor = hasAnchor || anchored[row] != 0;
            for (int neighbor : neighbors[row]) if (!visited[neighbor])
            {
                visited[neighbor] = 1; pending.push(neighbor);
            }
        }
        require(hasAnchor, "every active fluid component requires a positive-area pressure anchor; no gauge pin");
    }
    result.divergence_.resize(activeCount, static_cast<int>(result.nativeByVelocityRow_.size()));
    result.divergence_.setFromTriplets(divergenceEntries.begin(), divergenceEntries.end());
    result.divergence_.makeCompressed();
    if (result.uniformLeaves()) result.pressureGradientLoad_ = -result.divergence_.transpose();
    else
    {
        result.pressureGradientLoad_.resize(result.nativeByVelocityRow_.size(),activeCount);
        result.pressureGradientLoad_.setFromTriplets(pressureEntries.begin(),pressureEntries.end());
    }
    result.pressureGradientLoad_.makeCompressed();
    if (!result.uniformLeaves())
        preserveInteriorConstantMode(result.pressureGradientLoad_, result.faces_, result.nativeByVelocityRow_);
    return result;
}

void HydroFittedDomain::validateMeshIdentity(const QuadtreeMesh& mesh) const
{
    geometry_.validateMeshIdentity(mesh);
}

Vector2D HydroFittedDomain::spacing() const
{
    return geometry_.spacing();
}

namespace
{
bool sameBoundary(const HydroFittedDomain::Boundary& a, const HydroFittedDomain::Boundary& b)
{
    return a.type == b.type && a.id == b.id && a.pressure == b.pressure
        && a.incomingDensity == b.incomingDensity && a.incomingViscosity == b.incomingViscosity;
}

}

bool HydroFittedDomain::sameDefinition(const HydroFittedDomain& b) const
{
    const auto& a = *this;
    if (a.axisymmetric() != b.axisymmetric() || a.uniformLeaves()!=b.uniformLeaves() || a.cellRowByIndex() != b.cellRowByIndex()
        || a.cellIndexByRow() != b.cellIndexByRow() || a.velocityRowByNative() != b.velocityRowByNative()
        || a.nativeByVelocityRow() != b.nativeByVelocityRow() || a.faces().size() != b.faces().size()
        || a.cellVolumes() != b.cellVolumes()) return false;
    for (std::size_t i = 0; i < a.faces().size(); ++i)
    {
        const auto& x = a.faces()[i]; const auto& y = b.faces()[i];
        if (x.nativeRow != y.nativeRow || x.edgeLabel != y.edgeLabel || x.horizontal != y.horizontal
            || x.negativeCell != y.negativeCell || x.positiveCell != y.positiveCell || x.velocityRow != y.velocityRow
            || x.center.x != y.center.x || x.center.y != y.center.y || x.area != y.area
            || x.dualVolume != y.dualVolume || x.normalSpacing != y.normalSpacing
            || x.tangentSpacing != y.tangentSpacing || x.outwardSign != y.outwardSign
            || !sameBoundary(x.boundary, y.boundary)) return false;
    }
    return true;
}
